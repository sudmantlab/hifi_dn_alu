rule minimap2:
    input:
        hifi = "data/PacBio-HiFi/{specimen}/{lane}/{smrtcell}.fastq.gz"
    output: 
        temp("output/alignment/{refalias}/minimap2/standard/mapped/temp/{specimen}/{lane}/{smrtcell}.filt.bam")
    params:
        refgenome = config['reference']['fasta'],
        readgroup = config['minimap2']['readgroup'],
        minQ = config['samtools']['minQ']
    conda: "../envs/environment.yml"
    threads: 14
    shell: 
        """
        minimap2 --version && minimap2 {params.refgenome} {input.hifi} -t {threads} -ax map-hifi -Y -y -L --eqx --cs --MD -R '{params.readgroup}' | samtools view -q {params.minQ} -bT {params.refgenome} -o {output}
        """

rule minimap2_to_hg38_scaffolded:
    input:
        hifi = "data/PacBio-HiFi/{specimen}/{lane}/{smrtcell}.fastq.gz",
        fa = "output/assembly/hifiasm/{specimen}/hg38_scaffolded/{specimen}.diploid.fasta"
    output:
        temp("output/alignment/hg38_scaffolded/minimap2/standard/mapped/temp/{specimen}/{lane}/{smrtcell}.filt.bam")
    params:
        readgroup = config['minimap2']['readgroup'],
        # minQ = config['samtools']['minQ'] # skip minQ for multimapping
    conda: "../envs/environment.yml"
    threads: 14
    shell:
        """
        minimap2 --version && minimap2 -t {threads} -ax map-hifi -Y -y -L --eqx --cs -I8g --MD -R '{params.readgroup}' {input.fa} {input.hifi} | samtools view -b > {output}
        """