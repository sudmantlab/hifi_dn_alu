# Direct identification of de novo mobile element insertions from single molecule sequencing of human sperm
WIP: Paths might be wonky while testing that the structure generalizes in new setups. Watch your step :)

# Getting started

If not installed, install `mamba` (version ≥1.1.0). 
* `conda` (version ≥ 22.11.1) is also acceptable, but `mamba` will build and update environments faster.

Create a new environment with the minimal requirements by running:
```bash
mamba create -n {env_name} -f requirements.txt -c bioconda
```

Alternatively, install the minimal requirements directly into your base directory.

# Pipeline

## Configuration

Clone this repository and set it as your working directory.
```bash
git clone https://github.com/sudmantlab/hifi_dn_alu.git
cd hifi_dn_alu
```

Most of the pipeline is configured using the config file at `config/snakemake/hg38.config.yml` and a `samples.tsv` file describing your samples (see `config/snakemake/example_samples.tsv` for the expected format).

The pipeline expects the following directories to be present in the working directory:
* `data`: A directory containing uBAM (unaligned BAM) or `fastq.gz` files in a nested subdirectory structure (ex: `data/PacBio-HiFi/{specimen}/{lane}/{smrtcell}.fastq.gz`). This structure is necessary to provide multiple input files per specimen (see `config/snakemake/example_samples.tsv`).
    * If providing `.fastq.gz` files: Make sure that they contain the read quality and methylation tags in the header for each read.
    * If providing uBAM files: The `uBAMtoFastq` helper rule creates `.fastq.gz` files, preserving read quality and methylation tags from the uBAM.
    * Most sequencing facilities provide both lane and SMRTcell IDs, but you can provide placeholders if one or either are unknown (e.g. `data/PacBio-HiFi/your_specimen_name/lane_placeholder/actual_filename.fastq.gz`).
* `references`: A directory containing a reference FASTA and annotation files for the desired reference genome.
    * By default, this pipeline expects the [HGSVC no-ALT reference](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/HGSVC2/technical/reference/20200513_hg38_NoALT/), but can be modified to accept any reference that has a RepeatMasker annotation file (required for Sniffles2).
    * `refalias` will determine the name of the output directory (`output/{refalias}`).

## Dependencies
This pipeline is designed to require minimal dependency management, but you can update most dependencies from the `envs/environment.yml` file. `snakemake` hashes the contents of this file and checks for updates to `envs/environment.yml` at runtime and creates new environments as necessary, ensuring that your results are always up-to-date. 

## Running the pipeline

This pipeline was built and tested on the [Savio cluster at UC Berkeley](https://docs-research-it.berkeley.edu/services/high-performance-computing/user-guide/hardware-config/) using savio4_htc (56 cores @ Intel Xeon Gold 6330 @ 2.0 GHz & 512 GB RAM). You may need to update thread/memory allocation according to your resources.

```

You can test your setup using the dry run (`-np`) setting.
```bash
snakemake --use-conda --cores {cores} -np
```

Recommended: If working on a SLURM cluster, you can use `snakemake` to submit batch jobs to a specified partition (see `config/snakemake/example_slurm.yml` and [the docs](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html)). This allows `snakemake` to automatically re-queue and update memory allocation for jobs that fail with OOM errors. This is particularly useful if calibrating memory allocation or troubleshooting OOM errors on processes where memory allocation is determined by input size (typically assembly and alignment).

```
snakemake --use-conda --profile config/snakemake/example_slurm.yml -np
```