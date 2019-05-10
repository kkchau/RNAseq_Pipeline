import os
import subprocess
import re
import json
from datetime import datetime
from chunkypipes.components import Software, Parameter, Redirect, BasePipeline

FIRST_READS_PAIR = 0


class Pipeline(BasePipeline):
    @staticmethod
    def count_gzipped_lines(filepath):
        zcat = subprocess.Popen(['zcat', filepath], stdout=subprocess.PIPE)
        num_lines = subprocess.check_output(['wc', '-l'], stdin=zcat.stdout)
        return num_lines.strip()

    def description(self):
        return """This is an exact replication of the ENCODE long-rna pipeline."""

    def add_pipeline_args(self, parser):
        parser.add_argument('--reads', required=True, nargs='*',
                            help=('Reads to process with this pipeline. Denote paired-end reads with '
                                  'a colon (Ex. read1.fq:read2.fq). Specify multiple times to '
                                  'align multiple libraries (or pairs).'))
        parser.add_argument('--output', required=True,
                            help='Full path to output directory.')
        parser.add_argument('--lib', default=datetime.now().strftime('%Y-%m-%d-%H-%M-%S'),
                            help=('Name of the library, prepended to output file names. Defaults to '
                                  'a date string (YYYY-MM-DD-hh-mm-ss).'))
        parser.add_argument('--step', default=0,
                            help='Pipeline step to start at.')
        parser.add_argument('--forward-adapter', default='ZZZ',
                            help='Adapter sequence for the forward strand.')
        parser.add_argument('--reverse-adapter', default='ZZZ',
                            help='Adapter sequence for the reverse strand.')
        parser.add_argument('--is-stranded', action='store_true',
                            help='Provide this argument if library is stranded.')
        return parser

    def configure(self):
        return {
            'cutadapt': {
                'path': 'Full path to cutadapt',
                'quality-base': 'Phred quality scale [33|64]'
            },
            'STAR': {
                'path': 'Full path to STAR',
                'genome-dir': 'Directory containing a STAR genome index',
                'threads': 'Number of threads to run STAR'
            },
            'RSEM': {
                'path-calculate-expression': 'Full path to RSEM (rsem-calculate-expression)',
                'path-plot-model': 'Full path to RSEM (rsem-plot-model)',
                'reference-dir': 'Full path to RSEM reference, including prefix (Ex. /path/to/rsem-ref/ref-prefix)',
                'threads': 'Number of threads to run RSEM',
                'memory': 'Memory to use for RSEM in MB'
            },
            'bedgraph_to_bw': {
                'path': 'Full path to BedgraphToBW'
            },
            'samtools': {
                'path': 'Full path to samtools'
            },
            'sort': {
                'memory': 'Memory to use for sorting (Ex. 32G)'
            },
            'bedSort': {
                'path': 'Full path to bedSort'
            }
        }

    def run_pipeline(self, pipeline_args, pipeline_config):
        # Instantiate options
        reads = pipeline_args['reads']
        output_dir = pipeline_args['output']
        logs_dir = os.path.join(output_dir, 'logs')
        lib_prefix = pipeline_args['lib']
        step = pipeline_args['step']
        forward_adapter = pipeline_args['forward_adapter']
        reverse_adapter = pipeline_args['reverse_adapter']
        run_is_stranded = pipeline_args['is_stranded']

        # Determine if run is paired-end from input
        run_is_paired_end = len(reads[FIRST_READS_PAIR].split(':')) > 1

        # Create output, tmp, and logs directories
        subprocess.call(['mkdir', '-p', output_dir,
                         logs_dir, os.path.join(output_dir, 'tmp')])

        # Timing functions for getting running time
        start_time = datetime.now()

        # Gather QC data
        qc_data = {
            'total_raw_reads_counts': [],
            'trimmed_reads_counts': [],
            'num_reads_mapped': '0',
            'running_time_seconds': '',
            'running_time_readable': ''
        }

        # Keep list of items to delete
        staging_delete = [os.path.join(output_dir, 'tmp')]

        # Establish software instances
        cat = Software('cat', '/bin/cat')
        cutadapt = Software('cutadapt', pipeline_config['cutadapt']['path'])
        star = Software('STAR', pipeline_config['STAR']['path'])
        rsem_calculate_expression = Software('RSEM', pipeline_config['RSEM']['path-calculate-expression'])
        rsem_plot_model = Software('RSEM', pipeline_config['RSEM']['path-plot-model'])
        bedGraph_to_bw = Software('bedGraphToBigWig', pipeline_config['bedgraph_to_bw']['path'])
        bed_sort = Software('bedSort', pipeline_config['bedSort']['path'])
        samtools_flagstat = Software('samtools flagstat', pipeline_config['samtools']['path'] + ' flagstat')

        # Step 1: If more than one reads pairs are provided, combine them
        if step <= 1 and len(reads) >= 2:
            if run_is_paired_end:
                # Aggregate read1s and read2s
                read1s, read2s = [], []
                for reads_set in reads:
                    read1, read2 = reads_set.split(':')
                    read1s.append(read1)
                    read2s.append(read2)

                # Combine reads groups
                combined_reads = []
                for name, reads_group in [('read1', read1s), ('read2', read2s)]:
                    combined_read_filename = os.path.join(output_dir, '{}.combined.{}.fastq.gz'.format(lib_prefix, name))
                    combined_reads.append(combined_read_filename)
                    staging_delete.append(combined_read_filename)
                    cat.run(
                        Parameter(*[read for read in reads_group]),
                        Redirect(stream=Redirect.STDOUT, dest=combined_read_filename)
                    )

                # Update reads list
                reads = [':'.join(combined_reads)]
            else:
                # Combine reads
                combined_read_filename = os.path.join(output_dir, '{}.combined.fastq.gz'.format(lib_prefix))
                staging_delete.append(combined_read_filename)
                cat.run(
                    Parameter(*[read for read in reads]),
                    Redirect(stream=Redirect.STDOUT, dest=combined_read_filename)
                )

                # Update reads list
                reads = [combined_read_filename]

        # Step 2: Trim adapters with cutadapt
        if step <= 2:
            reads_set = reads[FIRST_READS_PAIR]
            if run_is_paired_end:
                # Get paired-end reads, construct new filenames
                read1, read2 = reads_set.split(':')

                # QC: Get raw fastq read counts
                qc_data['total_raw_reads_counts'].extend([
                    str(int(self.count_gzipped_lines(read1))/4),
                    str(int(self.count_gzipped_lines(read2))/4)
                ])

                trimmed_read1_filename = os.path.join(output_dir, lib_prefix + '_read1.trimmed.fastq.gz')
                trimmed_read2_filename = os.path.join(output_dir, lib_prefix + '_read2.trimmed.fastq.gz')

                staging_delete.append(trimmed_read1_filename)
                staging_delete.append(trimmed_read2_filename)

                # Run cutadapt
                cutadapt.run(
                    Parameter('--quality-base={}'.format(pipeline_config['cutadapt']['quality-base'])),
                    Parameter('--minimum-length=5'),
                    Parameter('--output={}'.format(trimmed_read1_filename)),
                    Parameter('--paired-output={}'.format(trimmed_read2_filename)),
                    Parameter('-a', forward_adapter),
                    Parameter('-A', reverse_adapter),
                    Parameter('-q', '30'),
                    Parameter(read1),
                    Parameter(read2),
                    Redirect(stream=Redirect.STDOUT, dest=os.path.join(logs_dir, 'cutadapt.summary.log'))
                )

                # QC: Get trimmed fastq read counts
                qc_data['trimmed_reads_counts'].extend([
                    str(int(self.count_gzipped_lines(trimmed_read1_filename))/4),
                    str(int(self.count_gzipped_lines(trimmed_read2_filename))/4)
                ])

                # Update reads list
                reads = ':'.join([trimmed_read1_filename, trimmed_read2_filename])

            else:
                # QC: Get raw fastq read count
                qc_data['total_raw_reads_counts'].append(
                    str(int(self.count_gzipped_lines(
                        os.path.join(output_dir, '{}.combined.fastq.gz'.format(lib_prefix))
                    ))/4)
                )

                # Construct new filename
                trimmed_read_filename = os.path.join(output_dir, lib_prefix + '.trimmed.fastq.gz')

                staging_delete.append(trimmed_read_filename)

                # Run cutadapt
                cutadapt.run(
                    Parameter('--quality-base={}'.format(pipeline_config['cutadapt']['quality-base'])),
                    Parameter('--minimum-length=5'),
                    Parameter('--output={}'.format(trimmed_read_filename)),
                    Parameter('-a', forward_adapter),
                    Parameter('-q', '30'),
                    Parameter(reads[FIRST_READS_PAIR]),
                    Redirect(stream=Redirect.STDOUT, dest=os.path.join(logs_dir, 'cutadapt.summary.log'))
                )

                # QC: Get trimmed fastq read count
                qc_data['trimmed_reads_counts'].append(
                    str(int(self.count_gzipped_lines(trimmed_read_filename))/4)
                )

                # Update reads list
                reads = [trimmed_read_filename]

        # Step 3: Alignment
        if step <= 3:
            # Gets reads for paired-end and single-end
            if run_is_paired_end:
                read1, read2 = reads.split(':')
            else:
                read1 = reads[FIRST_READS_PAIR]
                read2 = ''

            # Set up STAR parameters
            star_outfile_prefix = os.path.join(output_dir,
                                               lib_prefix + ('.' if lib_prefix[-1] != '.' else ''))
            star_common = [
                Parameter('--outFileNamePrefix', star_outfile_prefix),
                Parameter('--genomeDir', pipeline_config['STAR']['genome-dir']),
                Parameter('--readFilesIn', read1, read2),
                Parameter('--readFilesCommand', 'zcat'),
                Parameter('--outFilterType', 'BySJout'),
                Parameter('--outFilterMultimapNmax', '20'),
                Parameter('--alignSJoverhangMin', '8'),
                Parameter('--alignSJDBoverhangMin', '1'),
                Parameter('--outFilterMismatchNmax', '999'),
                Parameter('--alignIntronMin', '20'),
                Parameter('--alignIntronMax', '1000000'),
                Parameter('--alignMatesGapMax', '1000000'),
                Parameter('--outSAMunmapped', 'Within'),
                Parameter('--outSAMattributes', 'NH', 'HI', 'AS', 'NM', 'MD'),
                Parameter('--outFilterMismatchNoverReadLmax', '0.04'),
                Parameter('--sjdbScore', '1')
            ]

            star_run = [
                Parameter('--runThreadN', pipeline_config['STAR']['threads']),
                #Parameter('--genomeLoad', 'LoadAndKeep'),
                #Parameter('--limitBAMsortRAM', '10000000000')
            ]

            star_bam = [
                Parameter('--outSAMtype', 'BAM', 'SortedByCoordinate'),
                Parameter('--quantMode', 'TranscriptomeSAM')
            ]

            star_strand, star_wig = [], []

            # STAR strandedness parameters
            if run_is_stranded:
                star_wig.append(Parameter('--outWigStrand', 'Stranded'))
            else:
                star_strand.append(Parameter('--outSAMstrandField', 'intronMotif'))
                star_wig.append(Parameter('--outWigStrand', 'Unstranded'))

            star_meta = []

            # Run STAR alignment step
            star.run(*(star_common + star_run + star_bam + star_strand + star_meta))

            # Store STAR output files
            star_output_bam = star_outfile_prefix + 'Aligned.sortedByCoord.out.bam'

            # QC: Get samtools flagstat
            samtools_flagstat.run(
                Parameter(star_output_bam),
                Redirect(stream=Redirect.STDOUT, dest=star_output_bam + '.flagstat')
            )

            # QC: Get number of mapped reads from this BAM
            with open(star_output_bam + '.flagstat') as flagstats:
                flagstats_contents = flagstats.read()
                target_line = re.search(r'(\d+) \+ \d+ mapped', flagstats_contents)
                if target_line is not None:
                    qc_data['num_reads_mapped'] = str(int(target_line.group(1))/2)

            # Generate bedGraph
            signal_output_dir = os.path.join(output_dir, 'signal')
            subprocess.call(['mkdir', '-p', signal_output_dir])
            signal_output_prefix = os.path.join(signal_output_dir,
                                                lib_prefix + ('.' if lib_prefix[-1] != '.' else ''))

            # Run STAR for signal generation
            star.run(
                Parameter('--runMode', 'inputAlignmentsFromBAM'),
                Parameter('--inputBAMfile', star_output_bam),
                Parameter('--outWigType', 'bedGraph'),
                Parameter('--outFileNamePrefix', signal_output_prefix),
                Parameter('--outWigReferencesPrefix', 'chr'),
                *star_wig
            )

            # Convert bedGraph to bigWig
            chrNL_txt = os.path.join(output_dir, 'chrNL.txt')
            with open(chrNL_txt, 'w') as chrNL_filehandle:
                subprocess.call(['grep', '^chr',
                                 os.path.join(pipeline_config['STAR']['genome-dir'], 'chrNameLength.txt')
                                 ],
                                stdout=chrNL_filehandle)

            # Generate temporary signal file path
            sig_tmp = os.path.join(output_dir, 'sig.tmp')
            staging_delete.append(sig_tmp)
            if run_is_stranded:
                strand = [None, '-', '+']
                for i_strand in [1, 2]:
                    for i_mult in ['Unique', 'UniqueMultiple']:
                        # Get signal file for this iteration
                        signal_file = '{}Signal.{}.str{}.out.bg'.format(signal_output_prefix, i_mult, str(i_strand))
                        # Write to temporary signal file
                        with open(sig_tmp, 'w') as sig_tmp_filehandle:
                            subprocess.call(['grep', '^chr', signal_file],
                                            stdout=sig_tmp_filehandle)
                        # Sort sig.tmp with bedSort
                        bed_sort.run(
                            Parameter(sig_tmp),
                            Parameter(sig_tmp)
                        )
                        # Run bedGraph to bigWig conversion
                        bedGraph_to_bw.run(
                            Parameter(sig_tmp),
                            Parameter(chrNL_txt),
                            Parameter('{}Signal.{}.strand{}.bw'.format(
                                signal_output_prefix,i_mult, strand[i_strand]
                            ))
                        )
            else:
                for i_mult in ['Unique', 'UniqueMultiple']:
                    # Get signal file for this iteration
                    signal_file = '{}Signal.{}.str1.out.bg'.format(signal_output_prefix, i_mult)
                    # Write to temporary signal file
                    with open(sig_tmp, 'w') as sig_tmp_filehandle:
                        subprocess.call(['grep', '^chr', signal_file],
                                        stdout=sig_tmp_filehandle)
                    # Sort sig.tmp with bedSort
                    bed_sort.run(
                        Parameter(sig_tmp),
                        Parameter(sig_tmp)
                    )
                    # Run bedGraph to bigWig conversion
                    bedGraph_to_bw.run(
                        Parameter(sig_tmp),
                        Parameter(chrNL_txt),
                        Parameter('{}Signal.{}.unstranded.bw'.format(signal_output_prefix, i_mult))
                    )

        # Step 4: Sort transcriptome BAM to ensure order of reads to make RSEM output deterministic
        if step <= 4:
            # Set BAM file paths, mv transcriptome BAM to temporary name
            star_outfile_prefix = os.path.join(output_dir,
                                               lib_prefix + ('.' if lib_prefix[-1] != '.' else ''))
            transcriptome_bam = star_outfile_prefix + 'Aligned.toTranscriptome.out.bam'
            tr_bam = star_outfile_prefix + 'Tr.bam'
            staging_delete.append(tr_bam)
            subprocess.call(['mv', transcriptome_bam, tr_bam])

            # Template command
            merge_cmd = 'cat <({input1}) <({input2}) | {compress} > {output}'
            input1_cmd = '{samtools} view -H {bam}'
            compress_cmd = 'samtools view -@ {threads} -bS -'

            if run_is_paired_end:
                input2_cmd = ('{samtools} view -@ {threads} {bam} | ' +
                              'awk \'{{printf "%s", $0 " "; getline; print}}\' | ' +
                              'sort -S {ram} -T {tmpdir} | ' +
                              'tr \' \' \'\\n\'')
            else:
                input2_cmd = ('{samtools} view -@ {threads} {bam} | ' +
                              'sort -S {ram} -T {tmpdir}')

            print merge_cmd.format(
                input1=input1_cmd.format(
                    samtools=pipeline_config['samtools']['path'],
                    bam=tr_bam
                ),
                input2=input2_cmd.format(
                    samtools=pipeline_config['samtools']['path'],
                    threads=pipeline_config['RSEM']['threads'],
                    bam=tr_bam,
                    ram=pipeline_config['sort']['memory'],
                    tmpdir=os.path.join(output_dir, 'tmp')
                ),
                compress=compress_cmd.format(
                    threads=pipeline_config['RSEM']['threads']
                ),
                output=transcriptome_bam
            )

            subprocess.call(merge_cmd.format(
                input1=input1_cmd.format(
                    samtools=pipeline_config['samtools']['path'],
                    bam=tr_bam
                ),
                input2=input2_cmd.format(
                    samtools=pipeline_config['samtools']['path'],
                    threads=pipeline_config['RSEM']['threads'],
                    bam=tr_bam,
                    ram=pipeline_config['sort']['memory'],
                    tmpdir=os.path.join(output_dir, 'tmp')
                ),
                compress=compress_cmd.format(
                    threads=pipeline_config['RSEM']['threads']
                ),
                output=transcriptome_bam
            ), shell=True, executable='/bin/bash')

            subprocess.call(['rm', tr_bam])

        # Step 5: Run RSEM to get quantification
        if step <= 5:
            star_outfile_prefix = os.path.join(output_dir,
                                               lib_prefix + ('.' if lib_prefix[-1] != '.' else ''))
            transcriptome_bam = star_outfile_prefix + 'Aligned.toTranscriptome.out.bam'

            # Set up RSEM parameters
            rsem_common = [
                Parameter('--bam'),
                Parameter('--estimate-rspd'),
                Parameter('--calc-ci'),
                Parameter('--no-bam-output'),
                Parameter('--seed', '12345')
            ]

            rsem_run = [
                Parameter('-p', pipeline_config['RSEM']['threads']),
                Parameter('--ci-memory', pipeline_config['RSEM']['memory'])
            ]

            rsem_type = []
            if run_is_paired_end:
                rsem_type.append(Parameter('--paired-end'))
            if run_is_stranded:
                rsem_type.append(Parameter('--forward-prob', '0'))

            # Run RSEM quantification step
            rsem_calculate_expression.run(*(rsem_common + rsem_run + rsem_type + [
                Parameter(transcriptome_bam),
                Parameter(pipeline_config['RSEM']['reference-dir']),
                Parameter(os.path.join(output_dir, 'RSEM_Quant')),
                Redirect(Redirect.BOTH, dest=os.path.join(logs_dir, 'Log.rsem'))
            ]))

            # Generate RSEM plot model
            rsem_plot_model.run(
                Parameter(os.path.join(output_dir, 'RSEM_Quant'), os.path.join(output_dir, 'Quant.pdf'))
            )

        # QC: Get time delta
        elapsed_time = datetime.now() - start_time
        qc_data['running_time_seconds'] = str(elapsed_time.seconds)
        qc_data['running_time_readable'] = str(elapsed_time)

        # QC: Output QC data to file
        with open(os.path.join(logs_dir, 'qc_metrics.txt'), 'w') as qc_data_file:
            qc_data_file.write(json.dumps(qc_data, indent=4) + '\n')

        # Delete temporary files
        for delete_file in staging_delete:
            subprocess.call(['rm', '-rf', delete_file])

        print 'Complete'
        print 'Elapsed time: {}'.format(str(elapsed_time))
        print 'Elapsed time seconds: {}'.format(str(elapsed_time.seconds))
