from dataclasses import dataclass, field
from typing import List, Set, Dict, Tuple, Iterable, Optional
from pathlib import Path
import glob
import os
import numpy as np
import pandas as pd

configfile: "config/snakemake/config.yml"
workdir: config['workdir']
refalias = config['reference']['refalias']

# include common variables and helper functions
include: "rules/common.smk"

# preprocessing steps
include: "rules/preprocessing.smk"

# Assembly and QC
include: "rules/assembly.smk"

# Alignment (and realignment)
include: "rules/minimap2.smk"
include: "rules/samtools_utils.smk"
include: "rules/coverage_stats.smk"

# Variant calling
include: "rules/sniffles.smk"

# Preliminary graph variant calling
include: "rules/minigraph-cactus.smk"

ruleorder: minimap2_to_hg38_scaffolded > minimap2
ruleorder: sniffles_mosaic_scaffolded > sniffles_mosaic
ruleorder: sniffles_standard_scaffolded > sniffles_standard

rule all:
    input:
        # self-alignment: assembly + QC, variant calls through qc_all stage
        expand("output/assembly/flagger/{specimen}/prediction_summary_final.tsv", specimen = specimens),
        expand("output/assembly/hifiasm/{specimen}/quast/hg38_scaffolded/report.html", specimen = specimens),
        # self-alignment: coverage against the scaffolded personal assembly,
        # variant calls through the qc_all stage
        expand("output/alignment/{ref}_scaffolded/minimap2/standard/coverage_stats/{specimen}.coverage.html", ref = refalias, specimen = specimens),
        expand("output/alignment/{ref}_scaffolded/minimap2/standard/variants/sniffles_mosaic/{specimen}.qc_all.vcf.gz", ref = refalias, specimen = specimens)