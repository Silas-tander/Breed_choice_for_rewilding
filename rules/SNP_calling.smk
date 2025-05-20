rule merge_mustang_reads:
    input:
        r1=lambda wildcards: [f"{raw_file}_R1_001.fastq.gz" for raw_file in config["individuals"][wildcards.samples]],
        r2=lambda wildcards: [f"{raw_file}_R2_001.fastq.gz" for raw_file in config["individuals"][wildcards.samples]]
    output:
        r1="results/testing/merged_reads/{samples}_R1.fastq.gz",
        r2="results/testing/merged_reads/{samples}_R2.fastq.gz"
    threads: 2
    resources:
        mem_mb=16000,
        runtime=60
    shell:
        """
        cat {input.r1} > {output.r1}
        cat {input.r2} > {output.r2}
        """

rule bwa_mem:
    input:
        R1 = "results/testing/merged_reads/{sample}_R1.fastq.gz",
        R2 = "results/testing/merged_reads/{sample}_R2.fastq.gz",
        # Index can be a list of (all) files created by bwa, or one of them
        idx=multiext("data/ref_genome/ncbi_dataset/data/GCF_002863925.1/GCF_002863925.1_EquCab3.0_genomic.fna", "", ".amb", ".ann", ".bwt", ".pac", ".sa"),
    output:
        bam = "results/testing/mapped/{sample}.bam",
    threads: 6
    resources:
        mem_mb = 16000,
        runtime = 4*720
    conda:
        "../envs/bwa.yml"
    script:
        "../scripts/python/bwa_mem.py"

rule get_coverage:
    input:
        bam = "results/testing/mapped/{sample}.bam"
    output:
        stats = "results/testing/mapped_stats/{sample}.txt"
    conda:
        "../envs/bwa.yml"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 60
    shell:
        "samtools coverage -m {input.bam} > {output.stats}"

rule mark_duplicates:
    input:
        bam = "results/testing/mapped/{sample}.bam"
    params:
        mets = "results/testing/mapped_dups/{sample}.dups.bam.txt"
    output:
        dups = "results/testing/mapped_dups/{sample}.dups.bam"
    threads: 6
    resources:
        mem_mb = 8000,
        runtime = 320
    conda:
        "../envs/picard.yml"
    shell:
        '''
        picard MarkDuplicates -I {input.bam} -O {output.dups} -M {params.mets}
        '''

rule index_bam_dups:
    input:
        "results/testing/mapped_dups/{sample}.dups.bam"
    output:
        "results/testing/mapped_dups/{sample}.dups.bam.bai"
    threads: 2
    resources:
        mem_mb = 3000,
        runtime = 60
    conda:
        "../envs/bwa.yml"
    shell:
        """
        samtools index {input} > {output}
        """

rule ref_genome_chrom_filter:
    input:
        "data/ref_genome/ncbi_dataset/data/GCF_002863925.1/GCF_002863925.1_EquCab3.0_genomic.fna"
    output:
        "data/ref_genome/EquCab3.0_chrom.fa"
    params:
        chroms=" ".join(expected_chromosomes)
    conda:
        "../envs/bwa.yml"
    threads: 2
    resources:
        mem_mb = 3000,
        runtime = 10
    shell:
        '''
        samtools faidx {input} {params.chroms} > {output}
        '''

rule chrom_ref_index:
    input:
        "data/ref_genome/EquCab3.0_chrom.fa"
    output:
        "data/ref_genome/EquCab3.0_chrom.fa.fai"
    conda:
        "../envs/bwa.yml"
    threads: 2
    resources:
        mem_mb = 3000,
        runtime = 10
    shell:
        "samtools faidx {input} > {output}"

# create dict index for reference genome
rule index_reference_genome_dict_gatk:
    input:
        "data/ref_genome/ncbi_dataset/data/GCF_002863925.1/GCF_002863925.1_EquCab3.0_genomic.fna"
    output:
        "data/ref_genome/ncbi_dataset/data/GCF_002863925.1/GCF_002863925.1_EquCab3.0_genomic.dict"
    conda:
        "../envs/gatk.yaml"
    threads: 2
    resources:
        mem_mb = 3000,
        runtime = 10
    shell:
        "gatk CreateSequenceDictionary -R {input} -O {output}"

rule genmap_index:
    input:
        "data/ref_genome/EquCab3.0_chrom.fa"
    output:
        directory("data/ref_genome/genmap_idx/")
    conda:
        "../envs/genmap.yml"
    threads: 25
    resources:
        mem_mb = 50000,
        runtime = 720
    conda:
        "../envs/genmap.yml"
    shell:
        "genmap index -F {input} -I {output}"

rule mapability_mask:
    input:
        "data/ref_genome/genmap_idx/"
    output:
        directory("data/ref_genome/genmap_mask/")
    threads: 10
    resources:
        mem_mb = 80000,
        runtime = 720
    conda:
        "../envs/genmap.yml"
    shell:
        "genmap map -K 100 -E 2 -I {input} -O {output} -bg"

rule merge_bed:
    input:
        bed = "data/ref_genome/genmap_mask.bed"
    output:
        m_bed = "data/ref_genome/genmap_mask.merge.bed"
    threads: 2
    resources:
        mem_mb = 8000,
        runtime = 120
    conda:
        "../envs/bedtools.yml"
    shell:
        "bedtools merge -i {input.bed} > {output.m_bed}"

rule maskfasta:
    input:
        ref = "data/ref_genome/EquCab3.0_chrom.fa"
    params:
        bed = "data/ref_genome/genmap_mask.bed"
    output:
        mask_ref = "data/ref_genome/EquCab3.0_chrom.mask.fa"
    threads: 2
    resources:
        mem_mb = 8000,
        runtime = 120
    conda:
        "../envs/bedtools.yml"
    shell:
        '''
        bedtools merge -i {params.bed} | bedtools maskfasta -fi {input.ref} -bed - -fo {output.mask_ref}
        '''