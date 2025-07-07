configfile: "config.yaml"

# Rule 1: Create directories
rule create_dirs:
    output: "results/snakemake/.dirs_created"
    shell: """
        mkdir -p results/raw/SRR1972739 results/aligned results/variants results/annotated results/snpEff/data/reference_db
        touch {output}
    """

# Rule 2: Download fastq
rule download_fastq_from_ena:
    output:
        "results/raw/SRR1972739/SRR1972739_1.fastq.gz",
        "results/raw/SRR1972739/SRR1972739_2.fastq.gz"
    shell:
        """
        mkdir -p results/raw/SRR1972739
        wget -O {output[0]} ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR197/SRR1972739/SRR1972739_1.fastq.gz
        wget -O {output[1]} ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR197/SRR1972739/SRR1972739_2.fastq.gz
        """


# Rule 3: Convert SRA to FASTQ
rule sra_to_fastq:
    input: "results/raw/SRR1972739/SRR1972739.sra"
    output: "results/raw/SRR1972739/SRR1972739.fastq.gz"
    shell: "fastq-dump --split-files --gzip --outdir results/raw/SRR1972739 {input}"
# Rule 4: Align reads
rule align_reads:
    input:
        ref="results/raw/reference.fasta",
        reads="results/raw/SRR1972739/SRR1972739_1.fastq.gz"

    output: "results/aligned/SRR1972739.sam"
    shell: """
        bwa index {input.ref}
        bwa mem {input.ref} {input.reads} > {output}
    """

# Rule 5: Convert SAM to BAM
rule sam_to_bam:
    input: "results/aligned/SRR1972739.sam"
    output: "results/aligned/SRR1972739.bam"
    shell: "samtools view -bS {input} > {output}"

# Rule 6: Sort BAM
rule sort_bam:
    input: "results/aligned/SRR1972739.bam"
    output: "results/aligned/SRR1972739.sorted.bam"
    shell: "samtools sort {input} -o {output}"

# Rule 7: Variant Calling
rule call_variants:
    input:
        bam="results/aligned/SRR1972739.sorted.bam",
        ref="results/raw/reference.fasta"
    output: vcf="results/variants/raw_variants.vcf"
    shell: """
        bcftools mpileup -f {input.ref} {input.bam} | \
        bcftools call -mv -Ov -o {output.vcf}
    """

# Rule 8: Build SnpEff DB
rule build_snpeff_db:
    input: "results/raw/reference.fasta"
    output: touch("results/snpEff/data/reference_db/.db_built")
    shell: """
        mkdir -p results/snpEff/data/reference_db
        cp {input} results/snpEff/data/reference_db/sequences.fa
        cd results/snpEff && snpEff build -genbank -v reference_db
    """

# Rule 9: Annotate Variants
rule annotate_variants:
    input:
        vcf="results/variants/raw_variants.vcf",
        db="results/snpEff/data/reference_db/.db_built"
    output: "results/variants/annotated.vcf"
    shell: "snpEff -c snpEff.config reference_db {input.vcf} > {output}"

# Rule 10: Upload to S3
rule upload_s3:
    input: "results/variants/annotated.vcf"
    output: "results/snakemake/.s3_uploaded"
    shell: """
        aws s3 cp {input} s3://your-bucket-name/annotated.vcf
        touch {output}
    """

# Rule 11: Final target
rule all:
    input: "results/snakemake/.s3_uploaded"
