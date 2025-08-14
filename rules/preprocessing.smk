rule uBAMtoFastq:
    # Helper rule that takes an unaligned BAM file and transforms it into a fastq.
    # This transformation retains the read quality (rq) + MM & ML (methylation information) tags for each read.
    input:
        "data/PacBio-HiFi/{specimen}/{lane}/{smrtcell}.ccs.bam"
    output:
        "data/PacBio-HiFi/{specimen}/{lane}/{smrtcell}.fastq.gz"
    conda: "../envs/environment.yml"
    threads: 10
    shell:
        """
        samtools fastq -@ {threads} -c 6 -T MM,ML {input} -0 {output}
        """