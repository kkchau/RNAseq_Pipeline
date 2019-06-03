"""
__authors__: Kevin Chau, Hanqing Zhao
__date__: 2019-05-21
__description__: Perform RNA-Seq quality control, alignment, and quantification
"""

import os

configfile: "config.yaml"

# Setup pipeline directories
if not os.path.exists("data"):
    os.makedirs("data")
if not os.path.exists("log"):
    os.makedirs("log")

# Sub-Workflows
#subworkflow quantification:
#    workdir: "workflows/quantification"
#    snakefile: "workflows/quantification/Snakefile

#subworkflow align_QC:
#    workdir: "workflows/align_QC"
#    snakefile: "workflows/align_QC/Snakefile"

#subworkflow alignment:
    #workdir: "workflows/alignment"
    #snakefile: "workflows/alignment/Snakefile"

subworkflow quantification:
    workdir: "workflows/quantification"
    snakefile: "workflows/quantification/Snakefile"

rule all:
    input:
        quantification(os.path.join("../../data/RSEM_Quant", config["FASTQS"]["prefix"] + ".RSEM_Quant.rsem")),
        quantification("../../log/Log_out.rsem"),
        #quantification(os.path.join(config["GENOME_DIR"], "RSEM_INDEX/rsem_ref_prep.n2g.idx.fa")),
        #quantification(os.path.join(config["GENOME_DIR"], "RSEM_INDEX/rsem_ref_prep.idx.fa"))
        #alignment(os.path.join("../../data/STAR_ALIGN", config["FASTQS"]["prefix"] + ".STARAligned.sortedByCoord.out.bam")),
        #alignment(os.path.join("../../data/STAR_ALIGN", config["FASTQS"]["prefix"] + ".STARAligned.toTranscriptome.out_sorted.bam"))
#        align_QC("../../data/QC/"+config["FASTQS"]["prefix"]+".multiple_metrics.txt"),
#        align_QC("../../data/QC/"+config["FASTQS"]["prefix"]+".marked_duplicates.bam"),
#        align_QC("../../data/QC/"+config["FASTQS"]["prefix"]+".marked_dup_metrics.txt"),
#        align_QC("../../data/QC/"+config["FASTQS"]["prefix"]+".RNA_Metrics.txt"),
#        align_QC("../../data/QC/"+config["FASTQS"]["prefix"]+".insert_size_metrics.txt"),
#        align_QC("../../data/QC/"+config["FASTQS"]["prefix"]+".insert_size_histogram.pdf"),
#        align_QC("../../data/QC/"+config["FASTQS"]["prefix"]+".CollectalignmentSummaryMetricsoutput.txt"),
#        align_QC("../../data/QC/"+config["FASTQS"]["prefix"]+".GC_bias_metrics.pdf"),
#        align_QC("../../data/QC/"+config["FASTQS"]["prefix"]+".summary_metrics.txt"),
#        align_QC(os.path.join("../../data/QC", config["FASTQS"]["prefix"] + ".counts.txt")),
#        align_QC("../../log/"+config["FASTQS"]["prefix"]+".featureCounts.out"),
#        align_QC("../../data/STAR_ALIGN/"+config["FASTQS"]["prefix"]+".STARAligned.sortedByCoord.out.bam.bai"),
#        align_QC(os.path.join("../../data/QC", config["FASTQS"]["prefix"] + ".testReport/")),
#        align_QC("../../log/"+config["FASTQS"]["prefix"]+".RNA-SeQC.out"),
#        align_QC(os.path.join("../../log", config["FASTQS"]["prefix"] + ".Mt.txt"))
