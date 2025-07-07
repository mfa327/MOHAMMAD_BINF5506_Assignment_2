configfile: "config.yaml"

SRA = "SRR1972739"
REF = "results/raw/reference.fasta"

# Rule: Final output
rule all:
    input:
        "results/snakemake/.s3_uploaded"

# Rule 1: Create required directories
rule create_dirs:
    output:
        "results/snakemake/.dirs_created"
    shell:
        """
        mkdir -p results/raw/{SRA} results/aligned results/variants results/annotated results/snpEff/data/reference_db results/snakemake
        touch {output}
        """

# Rule 2: Download SRA file
rule download_sra:
    input:
        "results/snakemake/.dirs_created"
    output:
        "results/raw/SRR1972739/SRR1972739.sra"
    conda:
        "envs/sra_tools.yaml"
    shell:
        """
        prefetch SRR1972739 -O results/raw
        """

# Rule 3: Convert SRA to FASTQ
rule sra_to_fastq:
    input:
        "results/raw/SRR1972739/SRR1972739.sra"
    output:
        "results/raw/SRR1972739/SRR1972739_1.fastq.gz",
        "results/raw/SRR1972739/SRR1972739_2.fastq.gz"
    conda:
        "envs/sra_tools.yaml"
    shell:
        """
        fasterq-dump {input} -O results/raw/SRR1972739 --split-files
        gzip results/raw/SRR1972739/SRR1972739_1.fastq
        gzip results/raw/SRR1972739/SRR1972739_2.fastq
        """

# Rule 4: Align reads to reference
rule align_reads:
    input:
        ref="results/raw/reference.fasta",
        r1="results/raw/SRR1972739/SRR1972739_1.fastq.gz",
        r2="results/raw/SRR1972739/SRR1972739_2.fastq.gz"

    output:
        f"results/aligned/{SRA}.sam"
    conda:
        "envs/bwa.yaml"
    shell:
        """
        mkdir -p results/aligned
        bwa index {input.ref}
        bwa mem {input.ref} {input.r1} {input.r2} > {output}
        """
# Rule 5: Convert SAM to BAM
rule sam_to_bam:
    input:
        f"results/aligned/{SRA}.sam"
    output:
        f"results/aligned/{SRA}.bam"
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools view -bS {input} > {output}
        """

# Rule 6: Sort BAM file
rule sort_bam:
    input:
        f"results/aligned/{SRA}.bam"
    output:
        f"results/aligned/{SRA}.sorted.bam"
    conda:
        "envs/samtools.yaml"
    shell:
        """
        samtools sort {input} -o {output}
        """

# Rule 7: Call variants
rule call_variants:
    input:
        bam=f"results/aligned/{SRA}.sorted.bam",
        ref=REF
    output:
        vcf="results/variants/raw_variants.vcf"
    conda:
        "envs/bcftools.yaml"
    shell:
        """
        bcftools mpileup -f {input.ref} {input.bam} | \
        bcftools call -mv -Ov -o {output.vcf}
        """

# Rule 8: Build SnpEff database
rule build_snpeff_db:
    input:
        REF
    output:
        touch("results/snpEff/data/reference_db/.db_built")
    conda:
        "envs/snpeff.yaml"
    shell:
        """
        cp {input} results/snpEff/data/reference_db/sequences.fa
        cd results/snpEff
        snpEff build -noCheckProtein -genbank -v reference_db
        """
# Rule 9: Annotate variants
rule annotate_variants:
    input:
        vcf="results/variants/raw_variants.vcf",
        db_built="results/snpEff/data/reference_db/.db_built"
    output:
        "results/variants/annotated.vcf"
    conda:
        "envs/snpeff.yaml"
    shell:
        """
        snpEff -c snpEff.config reference_db {input.vcf} > {output}
        """

# Rule 10: Upload to S3
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
