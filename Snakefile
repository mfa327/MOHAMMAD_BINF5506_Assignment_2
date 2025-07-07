configfile: "config.yaml"

# 1
rule create_dirs:
    output:
        "results/snakemake/.dirs_created"
    shell:
        "mkdir -p results/raw/SRR1972739 results/aligned results/variants results/annotated results/snpEff/data/reference_db && touch {output}"

# 2
rule download_sra:
    input:
        "results/snakemake/.dirs_created"
    output:
        "results/raw/SRR1972739/SRR1972739.sra"
    shell:
        """
        prefetch SRR1972739 -O results/raw
        """
#3
rule sra_to_fastq:
    input:
        "results/raw/SRR1972739/SRR1972739.sra"
    output:
        "results/raw/SRR1972739/SRR1972739.fastq"
    shell:
        """
        fastq-dump --split-files --gzip --outdir results/raw/SRR1972739 results/raw/SRR1972739/SRR1972739.sra
        """
#4
rule align_reads:
    input:
        ref="results/raw/reference.fasta",
        reads="results/raw/SRR1972739/SRR1972739.fastq"
    output:
        "results/aligned/SRR1972739.sam"
    shell:
        """
        mkdir -p results/aligned
        bwa index {input.ref}
        bwa mem {input.ref} {input.reads} > {output}
        """
# 5
rule sam_to_bam:
    input:
        "results/aligned/SRR1972739.sam"
    output:
        "results/aligned/SRR1972739.bam"
    shell:
        """
        samtools view -bS {input} > {output}
        """
# 6
rule sort_bam:
    input:
        "results/aligned/SRR1972739.bam"
    output:
        "results/aligned/SRR1972739.sorted.bam"
    shell:
        """
        samtools sort {input} -o {output}
        """
# 7
rule call_variants:
    input:
        bam="results/aligned/SRR1972739.sorted.bam",
        ref="results/raw/reference.fasta"
    output:
        vcf="results/variants/raw_variants.vcf"
    shell:
        """
        bcftools mpileup -f {input.ref} {input.bam} | \
        bcftools call -mv -Ov -o {output.vcf}
        """
# 8
rule build_snpeff_db:
    input:
        "results/raw/reference.fasta"
    output:
        touch("results/snpEff/data/reference_db/.db_built")
    shell:
        """
        mkdir -p results/snpEff/data/reference_db
        cp {input} results/snpEff/data/reference_db/sequences.fa
        cd results/snpEff
        snpEff build -genbank -v reference_db
        """
# 9
rule annotate_variants:
    input:
        vcf="results/aligned/variants.raw.vcf",
        db_built="results/snpEff/data/reference_db/.db_built"
    output:
        "results/variants/annotated.vcf"
    shell:
        """
        snpEff -c snpEff.config reference_db {input.vcf} > {output}
        """
# 10
rule upload_s3:
    input:
        "results/variants/annotated.vcf"
    output:
        "results/snakemake/.s3_uploaded"
    shell:
        """
        aws s3 cp {input} s3://mohammad-assignment-2/annotated.vcf
        touch {output}
        """

# 11
rule all:
    input:
        "results/snakemake/.s3_uploaded"
