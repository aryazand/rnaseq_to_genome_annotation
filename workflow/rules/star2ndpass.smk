rule star_2ndpass:
    input:
        fq1=lambda wildcards: get_processed_fastq(wildcards.sample, regex="read1"),
        fq2=lambda wildcards: (
            get_processed_fastq(wildcards.sample, regex="read2")
            if is_paired_end()
            else []
        ),
        idx=rules.fastq_process_align_star_index.output,
    output:
        aln="results/star/2ndpass/{sample}/mapped.bam",
        log_final="results/star/2ndpass/{sample}/Log.final.out",
    log:
        "results/star/align/{sample}/mapped_2pndpass.log",
    message:
        "make star alignment"
    params:
        extra=lambda wildcards: config["mapping"]["star"]["extra"]
        + "--sjdbFileChrStartEnd "
        + " ".join(
            expand("results/star/align/{sample}/SJ.out.tab", sample=wildcards.sample)
        ),
    threads: 8
    wrapper:
        "v7.2.0/bio/star/align"


rule samtools_sort_2ndpass:
    input:
        "results/star/2ndpass/{sample}/mapped.bam",
    output:
        "results/samtools/sort_2ndpass/{sample}.bam",
    log:
        "results/samtools/sort_2ndpass/{sample}.log",
    message:
        "re-sort reads after mapping 2nd pass"
    params:
        extra=config["mapping"]["samtools_sort"]["extra"],
    threads: 2
    wrapper:
        "v7.0.0/bio/samtools/sort"


rule samtools_index_2ndpass:
    input:
        "results/samtools/sort_2ndpass/{sample}.bam",
    output:
        "results/samtools/sort_2ndpass/{sample}.bai",
    log:
        "results/star/align/{sample}/index_2ndpass.log",
    message:
        "index reads"
    params:
        extra=config["mapping"]["samtools_index"]["extra"],
    threads: 2
    wrapper:
        "v7.0.0/bio/samtools/index"


rule samtools_filter_2ndpass:
    input:
        "results/samtools/sort_2ndpass/{sample}.bam",
    output:
        bam="results/samtools/filtered/{sample}.bam",
        idx="results/samtools/filtered/{sample}.bai",
    log:
        "results/samtools/filtered/{sample}.log",
    params:
        extra="",  # optional params string
        region=config["mapping"]["samtools_view"]["region"],  # optional region string
    threads: 2
    wrapper:
        "v8.1.1/bio/samtools/view"
