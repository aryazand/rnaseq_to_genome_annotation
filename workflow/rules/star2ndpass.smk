rule star_2ndpass:
    input:
        fq1=lambda wildcards: fastq_process_align.get_processed_fastq(wildcards.sample, regex="read1"),
        fq2=lambda wildcards: (
            fastq_process_align.get_processed_fastq(wildcards.sample, regex="read2")
            if fastq_process_align.is_paired_end()
            else []
        ),
        idx=rules.fastq_process_align_star_index.output,
        aln="results/star/align/{sample}/mapped.bam",
        sj="results/star/align/{sample}/SJ.out.tab"
    output:
        aln="results/star/align/{sample}/mapped_2ndpass.bam",
        log_final="results/star/align/{sample}/Log_2ndpass.final.out",
    log:
        "results/star/align/{sample}/mapped_2ndpass.log",
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
