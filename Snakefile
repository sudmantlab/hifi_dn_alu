from dataclasses import dataclass, field
from typing import List, Set, Dict, Tuple, Iterable, Optional
from pathlib import Path
import glob
import os
import numpy as np
import pandas as pd

configfile: "config/snakemake/hg38.config.yml"
workdir: config['workdir']
localrules: make_bam_symlink, make_fastq_symlink

### common variables to be accessed in other rules/helper functions ###
sample_table = pd.read_table(config['sample_table'], index_col=False, dtype=str)
specimens = sample_table['specimen'].unique()
chrs = ['chr' + str(n) for n in np.arange(1, 22).tolist()+['X', 'Y']]

# include helper functions
include: "rules/common.smk"

# preprocessing steps
include: "rules/preprocessing.smk"

# Assembly and QC
include: "rules/assembly.smk"
include: "rules/get_assembly_stats.smk"

# Alignment (and realignment)
include: "rules/minimap2.smk"
include: "rules/duplomap.smk"
include: "rules/samtools_utils.smk"
include: "rules/coverage_stats.smk"

# Variant calling
include: "rules/sniffles.smk"

# Preliminary graph variant calling
include: "rules/minigraph-cactus.smk"

ruleorder: minimap2_to_hg38_scaffolded > minimap2
ruleorder: minimap2_to_T2T_scaffolded > minimap2
ruleorder: sniffles_mosaic_scaffolded > sniffles_mosaic
ruleorder: sniffles_standard_scaffolded > sniffles_standard
ruleorder: sniffles_mosaic_scaffolded > sniffles_mosaic
ruleorder: sniffles_standard_scaffolded > sniffles_standard
ruleorder: preprocess_scaffolded_variants > preprocess_variants

# Annotation
include: "rules/process_mosaic.smk"

rule all:
    input:
        # self-alignment: assembly + QC, variant calls through qc_all stage
        expand("output/assembly/flagger/{specimen}/prediction_summary_final.tsv", specimen = specimens),
        expand("output/assembly/hifiasm/{specimen}/quast/hg38_scaffolded/report.html", specimen = specimens),
        # hg38 alignment: reference coverage, variant calls through qc_all stage
        expand("output/alignment/hg38/minimap2/standard/coverage_stats/{specimen}.coverage.html", specimen = specimens),
        expand('output/alignment/hg38/minimap2/standard/variants/sniffles_mosaic/{specimen}.qc_all.vcf.gz', specimen = specimens)