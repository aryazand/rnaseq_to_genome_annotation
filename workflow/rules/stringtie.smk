rule stringtie:
    input:
        bam=fastq_process_align.get_bam_2,
    output:
        "results/stringtie/{sample}.gtf",
    conda:
        "../envs/stringtie.yml"
    log:
        "logs/stringtie/{sample}.log",
    params:
        ref_gff=config["annotation"]["stringtie"]["ref_anno"],
        strandedness=lambda wildcards: (
            "--rf"
            if config["annotation"]["stringtie"]["strandedness"] == "reverse"
            else (
                "--fr"
                if config["annotation"]["stringtie"]["strandedness"] == "forward"
                else ""
            )
        ),
    shell:
        "stringtie {input.bam} -G {params.ref_gff} {params.strandedness} -o {output}"


rule stringtie_merge:
    input:
        gtfs=expand("results/stringtie/{sample}.gtf", sample = fastq_process_align.samples.index),
    output:
        gtflist=temp("stringtie_gtflist.txt"),
        gtf="results/stringtie/stringtie_merged.gtf",
    conda:
        "../envs/stringtie.yml"
    log:
        "logs/stringtie/stringtie_merge.log",
    params:
        ref_gff=config["annotation"]["stringtie"]["ref_anno"],
    shell:
        """
       for sample in {input.gtfs}; do echo $sample >> {output.gtflist};done
       stringtie -p --merge {output.gtflist} -G {params.ref_gff} -o {output.gtf}
       """
