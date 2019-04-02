# RNA-Seq Pipeline
Snakemake pipeline optimized for usage on SLURM

```bash
conda env create -f conda_rnaseq.yaml
source activate rnaseq
```

---

Directory structure:

```
RNAseq_Pipeline  
|- source/              # Source data files (e.g. FASTQ)  
|- references/          # Reference files (e.g. genome FASTA, STAR index)  
|- workflows/           # Individual snakemake workflows  
    |- index/  
    |- align_and_quant/  
    |- ...  
|- cluster.yaml         # Cluster configuration (SLURM)
|- conda_rnaseq.yaml    # Conda environment
|- README.md  
|- Snakefile            # Main Snakefile
```