rule get_chr_sizes:
    input:
        "results/get_genome/genome.fasta.fai",
    output:
        "results/get_genome/genome.chrom.sizes",
    log:
        "results/get_genome/chrom_sizes.log",
    shell:
        """
        cut -f1,2 {input} > {output}
        """


rule stringtie_gtf_to_genePred:
    input:
        gtf="results/stringtie/stringtie_merged.gtf",
    output:
        genePred="results/stringtie/stringtie_merged.GenePred",
    params:
        extra=config["ucsc_trackhub"]["process_stringtie"]["gtf_to_genePred"],
    conda:
        "../envs/ucsc_tools.yaml"
    log:
        "results/stringtie/bigbed.log",
    shell:
        """
        gtfToGenePred {params.extra} {input.gtf} {output.genePred}
        """


rule stringtie_GenePred_to_bgpInput:
    input:
        GenePred="results/stringtie/stringtie_merged.GenePred",
    output:
        bgpInput="results/stringtie/stringtie_merged.bgpInput",
    params:
        extra=config["ucsc_trackhub"]["process_stringtie"]["GenePred_to_bgpInput"],
    conda:
        "../envs/ucsc_tools.yaml"
    log:
        "results/stringtie/GenePred_to_bgpInput.log",
    shell:
        """
         genePredToBigGenePred {params.extra} {input.GenePred} stdout | sort -k1,1 -k2,2n > {output.bgpInput}
        """


rule stringtie_bgpInput_to_bigGenePred:
    input:
        bgpInput="results/stringtie/stringtie_merged.bgpInput",
        chrom_sizes="results/get_genome/genome.chrom.sizes",
    output:
        bigGenePred="results/stringtie/stringtie_merged.bb",
    params:
        extra=config["ucsc_trackhub"]["process_stringtie"]["bgpInput_to_bigGenePred"],
    conda:
        "../envs/ucsc_tools.yaml"
    log:
        "results/stringtie/bgpInput_to_bigGenePred.log",
    shell:
        """
        bedToBigBed {params.extra} {input.bgpInput} {input.chrom_sizes} {output.bigGenePred}
        """


rule gff3_to_GenePred:
    input:
        gff="results/get_genome/genome.gff",
    output:
        genePred="results/get_genome/genome.GenePred",
    params:
        extra=config["ucsc_trackhub"]["process_genome_annotation"]["gff3_to_GenePred"],
    conda:
        "../envs/ucsc_tools.yaml"
    log:
        "results/get_genome/gff3_to_GenePred.log",
    shell:
        """
        gff3ToGenePred {params.extra} {input.gff} {output.genePred}
        """


rule GenePred_to_bgpInput:
    input:
        GenePred="results/get_genome/genome.GenePred",
    output:
        bgpInput="results/get_genome/genome.bgpInput",
    params:
        extra=config["ucsc_trackhub"]["process_genome_annotation"][
            "GenePred_to_bgpInput"
        ],
    conda:
        "../envs/ucsc_tools.yaml"
    log:
        "results/get_genome/GenePred_to_bgpInput.log",
    shell:
        """
         genePredToBigGenePred {params.extra} {input.GenePred} stdout | sort -k1,1 -k2,2n > {output.bgpInput}
        """


rule bgpInput_to_bigGenePred:
    input:
        bgpInput="results/get_genome/genome.bgpInput",
        chrom_sizes="results/get_genome/genome.chrom.sizes",
    output:
        bigGenePred="results/get_genome/genome.bb",
    params:
        extra=config["ucsc_trackhub"]["process_genome_annotation"][
            "bgpInput_to_bigGenePred"
        ],
    conda:
        "../envs/ucsc_tools.yaml"
    log:
        "results/get_genome/bgpInput_to_bigGenePred.log",
    shell:
        """
        bedToBigBed {params.extra} {input.bgpInput} {input.chrom_sizes} {output.bigGenePred}
        """


rule faToTwoBit:
    input:
        "results/get_genome/genome.fasta",
    output:
        "results/get_genome/genome.2bit",
    params:
        extra=config["ucsc_trackhub"]["process_fasta"]["faToTwoBit"],
    log:
        "results/get_genome/fa_to_2bit.log",
    wrapper:
        "v7.1.0/bio/ucsc/faToTwoBit"


rule ucsc_trackhub:
    input:
        genome_2bit="results/get_genome/genome.2bit",
        genome_genePred="results/get_genome/genome.bb",
        stringtie="results/stringtie/stringtie_merged.bb",
        plus_bw=lambda wildcards: expand(
            "results/deeptools/coverage/{sample}_forward.bw",
            sample=fastq_process_align.samples.index,
        ),
        minus_bw=lambda wildcards: expand(
            "results/deeptools/coverage/{sample}_reverse.bw",
            sample=fastq_process_align.samples.index,
        ),
    output:
        dir=directory(trackhub_dir),
        trackdb=expand(
            os.path.join(trackhub_dir, "{org}", "trackDb.txt"),
            org=config["ucsc_trackhub"]["genomes"],
        ),
        hubtxt=os.path.join(
            trackhub_dir, config["ucsc_trackhub"]["hub_name"] + ".hub.txt"
        ),
        genomestxt=os.path.join(
            trackhub_dir, config["ucsc_trackhub"]["hub_name"] + ".genomes.txt"
        ),
    conda:
        "../envs/trackhub.yaml"
    log:
        os.path.join(trackhub_dir, "trackhub.log"),
    script:
        "../scripts/trackhub.py"
