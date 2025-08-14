def get_temp_bams_per_sample(wildcards, sample_table = config['sample_table']) -> List[str]:
    """
    A helper function for creating a list of expected per-smrtcell bam outputs per sample to be merged into a single bam.
    Adjusted to optionally include {hap} in the basepath if wildcards.hap exists.
    """
    table = pd.read_table(sample_table, index_col=False, dtype=str)
    samples = table[table["specimen"] == str(wildcards.specimen)]
    samples = samples.to_records(index=False)
    
    # Check if '{hap}' is in the basepath and format accordingly
    if hasattr(wildcards, 'hap'):
        basepath = "output/alignment/{refalias}/{mapper}/standard/mapped/temp/{specimen}.{hap}/{lane}/{specimen}_{smrtcell}.filt.sorted.bam"
        input_samples = [basepath.format(refalias=wildcards.refalias, mapper=wildcards.mapper, specimen=s[0], hap=wildcards.hap, lane=s[2], smrtcell=s[3]) for s in samples]
    else:
        basepath = "output/alignment/{refalias}/{mapper}/standard/mapped/temp/{specimen}/{lane}/{specimen}_{smrtcell}.filt.sorted.bam"
        input_samples = [basepath.format(refalias=wildcards.refalias, mapper=wildcards.mapper, specimen=s[0], lane=s[2], smrtcell=s[3]) for s in samples]
    
    if len(input_samples) == 0:
        raise Exception("No samples found for specimen {}. Check samples.tsv and try again!".format(wildcards.specimen))
    else:
        return input_samples

def get_fastqs_per_sample(wildcards, sample_table = config['sample_table'], basepath = config['fastq_basepath']) -> List:
    """
    A helper function for creating a list of fastq.gz paths for each sample ahead of hifiasm assembly.
    """
    table = pd.read_table(sample_table, index_col=False, dtype=str)
    samples = table[table["specimen"] == str(wildcards.specimen)]
    samples = samples.to_records(index=False)
    input_samples = [basepath.format(specimen=s[0], group = s[1], lane=s[2], smrtcell = s[3]) for s in samples]
    if len(input_samples) == 0:
        raise Exception("No samples found for specimen {}. Check samples.tsv and try again!".format(wildcards.specimen))
    else:
        return input_samples

rule samtools_sort:
    group: "mapping"
    # TODO: Just change all files that have the .filt.sorted suffixes and just Not Do That
    input: "output/alignment/{refalias}/{mapper}/standard/mapped/temp/{specimen}/{lane}/{smrtcell}.filt.bam"
    output: 
        temp("output/alignment/{refalias}/{mapper}/standard/mapped/temp/{specimen}/{lane}/{specimen}_{smrtcell}.filt.sorted.bam")
    wildcard_constraints:
        specimen = '[A-Za-z0-9]+'
    threads: 10
    conda: "../envs/environment.yml"
    shell: "samtools sort -@ {threads} --output-fmt='BAM' -o {output} {input}"

rule collate_bams:
    group: "mapping"
    input: get_temp_bams_per_sample
    wildcard_constraints:
        specimen = '[A-Za-z0-9]+'
    output: "output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam"
    threads: 10
    conda: "../envs/environment.yml"
    shell: "samtools merge -r -@ {threads} --output-fmt='BAM' {output} {input}"

rule index_bam:
    group: "mapping"
    input: "output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam"
    output: "output/alignment/{refalias}/{mapper}/standard/mapped/{specimen}.sorted.merged.bam.bai"
    wildcard_constraints:
        specimen = '[A-Za-z0-9]+'
    threads: 10
    conda: "../envs/environment.yml"
    shell: 
        """
        samtools index -b {input} -@ {threads}
        """
