rule stringtie:
    input:
        bam="results/samtools/filtered/{sample}.bam",
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
        "stringtie {input.bam} -G {input.ref_gff} {params.strandedness} -o {output}"


rule stringtie_merge:
    input:
        gtfs=lambda wildcards: expand(
            "results/stringtie/{sample}.gtf", sample=wildcards.sample
        ),
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
       stringtie -p {resources.cpus_per_task} --merge {output.gtflist} -G {input.ref_gff} -o {output.gtf}
       """
