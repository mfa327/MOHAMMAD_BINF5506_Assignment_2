configfile: "config.yaml"

rule all:
    input:
        "results/s3_upload.done"

rule create_dirs:
    output:
        "results/snakemake/.dirs_created"
    shell:
        "mkdir -p results/raw/SRR1972739 results/aligned results/variants results/annotated results/snpEff/data/reference_db && touch {output}"

rule download_sra:
    input:
        "results/snakemake/.dirs_created"
    output:
        "results/raw/SRR1972739/SRR1972739.sra"
    shell:
        """
        prefetch SRR1972739 -O results/raw
        """

rule build_snpeff_db:
    input:
        "results/raw/reference.fasta"
    output:
        "results/snpEff/data/reference_db/.db_built"
    shell:
        """
        cp {input} results/snpEff/data/reference_db/sequences.fa
        cd results/snpEff
        snpEff build -genbank -noCheckProtein -v reference_db
        touch {output}
        """

rule align_reads:
    input:
        ref="results/raw/reference.fasta",
        fastq="results/raw/SRR1972739/SRR1972739_1.fastq.gz"
    output:
        "results/aligned/SRR1972739.sam"
    conda:
        "envs/bwa.yaml"
    shell:
        """
        bwa index {input.ref}
        bwa mem {input.ref} {input.fastq} > {output}
        """

rule sam_to_bam:
    input:
        "results/aligned/SRR1972739.sam"
    output:
        "results/aligned/SRR1972739.bam"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools view -bS {input} > {output}"

rule sort_bam:
    input:
        "results/aligned/SRR1972739.bam"
    output:
        "results/aligned/SRR1972739.sorted.bam"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools sort {input} -o {output}"

rule call_variants:
    input:
        bam="results/aligned/SRR1972739.sorted.bam",
        ref="results/raw/reference.fasta"
    output:
        "results/variants/variants.vcf"
    conda:
        "envs/bcftools.yaml"
    shell:
        """
        samtools mpileup -uf {input.ref} {input.bam} | bcftools call -mv -Ov -o {output}
        """

rule annotate_variants:
    input:
        vcf="results/variants/variants.vcf",
        db="results/snpEff/data/reference_db/.db_built"
    output:
        "results/annotated/variants_annotated.vcf"
    conda:
        "envs/snpeff.yaml"
    shell:
        """
        snpEff reference_db {input.vcf} > {output}
        """

rule upload_s3:
    input:
        "results/annotated/variants_annotated.vcf"
    output:
        "results/s3_upload.done"
    shell:
        """
        aws s3 cp {input} s3://your-bucket-name/variants_annotated.vcf
        touch {output}
        """
