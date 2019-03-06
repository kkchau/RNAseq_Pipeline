"""
__authors__: Kevin Chau
__date__: 2019-03-05
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

subworkflow align_and_quant:
    workdir: "workflows/align_and_quant"
    snakefile: "workflows/align_and_quant/Snakefile"
    configfile: "workflows/align_and_quant/config.yaml"