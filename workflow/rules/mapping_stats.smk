use rule samtools_sort from fastq_process_align as fastq_process_align_samtools_sort with:
    input:
        "results/star/align/{sample}/mapped_2ndpass.bam",


rule samtools_filter:
    input:
        "results/star/aligned/{sample}/mapped_2ndpass.bam",
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
