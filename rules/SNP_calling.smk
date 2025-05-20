# # Merge reads for each individual
# rule merge_reads:
#     input:
#         r1=lambda wildcards: [f"/home/silastander/Coregonus/Silas/Horses/data/mustang_fq/24121Bri_N24194/{accession}_R1_001.fastq.gz" for accession in config[wildcards.sample]],
#         r2=lambda wildcards: [f"/home/silastander/Coregonus/Silas/Horses/data/mustang_fq/24121Bri_N24194/{accession}_R2_001.fastq.gz" for accession in config[wildcards.sample]]
#     output:
#         r1="results/testing/merged_reads/{sample}_R1.fastq.gz",
#         r2="results/testing/merged_reads/{sample}_R2.fastq.gz"
#     threads: 4
#     resources:
#         mem_mb = 8000,
#         runtime = 60
#     shell:
#         """
#         cat {input.r1} > {output.r1}
#         cat {input.r2} > {output.r2}
#         """

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


rule trim_galore:
    input:
        R1 = "results/merged_reads/{sample}_1.fastq.gz",
        R2 = "results/merged_reads/{sample}_2.fastq.gz"
    output:
        trimmed1 = "results/trimmed/{sample}_1_val_1.fq.gz",
        trimmed2 = "results/trimmed/{sample}_2_val_2.fq.gz",
        report1 = temp("results/trimmed/{sample}_1.fastq.gz_trimming_report.txt"),
        report2 = temp("results/trimmed/{sample}_2.fastq.gz_trimming_report.txt")
    conda:
       os.path.join(workflow.basedir, "envs/trim_galore.yaml")
    threads: 2
    resources:
        mem_mb = 2000,
        runtime = 420
    shell:
        "trim_galore --paired -q 20 {input.R1} {input.R2} -o results/trimmed/"

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

# rule bwa_mem_pipe:
#     input:
#         uBAM = "data/uBAM_marked/{test_accession}.marked.bam",
#         # Index can be a list of (all) files created by bwa, or one of them
#         idx=multiext("../Horse_ref_genome/EquCab3.0.fa", "", ".amb", ".ann", ".bwt", ".pac", ".sa"),
#     output:
#         bam = "results/mapped/{test_accession}.ubam.bam",
#     threads: 6
#     resources:
#         mem_mb = 12000,
#         runtime = 60
#     conda:
#         "../envs/picard.yml"
#     script:
#         "../scripts/python/bwa_pipe.py"

# rule gather_bams:
#     input:
#         bams = expand("results/mapped/{test_accession}",
#             test_accession = config["test_accession"])
#     output:
#         gath_bam = "results/gath_mapped/{test_samples}.bam"
#     threads: 4
#     resources:
#         mem_mb = 4000,
#         runtime = 60
#     conda:
#         "../envs/picard.yml"
#     script:
#         "../scripts/python/gatherbams.py"

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


expected_chromosomes = ['NC_009144.3','NC_009145.3','NC_009146.3','NC_009147.3',
                        'NC_009148.3','NC_009149.3','NC_009150.3','NC_009151.3',
                        'NC_009152.3','NC_009153.3','NC_009154.3','NC_009155.3',
                        'NC_009156.3','NC_009157.3','NC_009158.3','NC_009159.3',
                        'NC_009160.3','NC_009161.3','NC_009162.3','NC_009163.3',
                        'NC_009164.3','NC_009165.3','NC_009166.3','NC_009167.3',
                        'NC_009168.3','NC_009169.3','NC_009170.3','NC_009171.3',
                        'NC_009172.3','NC_009173.3','NC_009174.3','NC_009175.3']

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

## This one does not work, I did it in the terminal, however the bed file is necessary for next step
# rule create_bed:
#     input:
#         fai = "data/ref_genome/EquCab3.0_chrom.fa.fai"
#     output:
#         bed = "data/ref_genome/EquCab3.0_chrom.bed"
#     threads: 2
#     resources:
#         mem_mb = 3000,
#         runtime = 10
#     shell:
#         '''
#         awk -v FS="\t" -v OFS="\t" '{print $1 FS "0" FS ($2-1)}' {input.fai} > {output.bed}
#         '''

rule bed_to_interval:
    input:
        bed = "data/ref_genome/EquCab3.0_chrom.bed"
    params:
        gen_dict = "data/ref_genome/EquCab3.0_chrom.dict"
    output:
        int_list = "data/ref_genome/EquCab3.0_chrom.interval_list"
    conda:
        "../envs/picard.yml"
    threads: 4
    resources:
        mem_mb = 8000,
        runtime = 60
    shell:
        '''
        picard BedToIntervalList -I {input.bed} -O {output.int_list} -SD {params.gen_dict}
        '''

rule bcftools_call:
    input:
        bam="results/mapped_dups/{sample}.dups.bam",
        ref="data/ref_genome/EquCab3.0_chrom.fa"
    output:
        bc_call = "results/variants/bcftools_call/{sample}_chrom.bcf"
    threads: 6
    resources:
        mem_mb = 20000,
        runtime = 720
    conda:
        "../envs/bcftools.yml"
    shell:
        '''
        bcftools mpileup -f {input.ref} {input.bam} | bcftools call -mv -Ob -o {output.bc_call}
        '''

# # Splitintervals rule from gatk to parralelize the variant calling
# rule splitintervals:
#     input:
#         ref="data/ref_genome/EquCab3.0_chrom.fa"
#     params:
#         scatter = 32,
#         chroms = "data/ref_genome/EquCab3.0_chrom.interval_list"
#     output:
#         bed = directory("data/ref_genome/splitintervals/")
#         # bed = multiext("../ref_genome/splitinterval/genome", ".00.bed", ".01.bed", ".02.bed", ".03.bed")
#     conda:
#         "../envs/gatk.yaml"
#     threads: 4
#     resources:
#         mem_mb = 8000,
#         runtime = 60
#     script:
#         "../scripts/python/splitintervals.py"

rule haplotype_caller:
    input:
        bam="results/testing/mapped_dups/{sample}.dups.bam",
        bai="results/testing/mapped_dups/{sample}.dups.bam.bai",
        ref="data/ref_genome/EquCab3.0_chrom.fa",
        intervals="data/ref_genome/splitintervals/{interval}",
        dict = "data/ref_genome/EquCab3.0_chrom.dict"
    output:
        vcf="results/testing/gatk_variants/{sample}/{interval}.g.vcf"
    params:
        extra="",  # optional
        java_opts="",  # optional
        ploidy = 2
    threads: 2
    resources:
        mem_mb = 4000,
        runtime = 120
    conda:
        "../envs/gatk.yaml"
    script:
        "../scripts/python/haplotype_caller_intervals.py"

rule gather_split_gvcf:
    input:
        gvcf_ints = expand("results/testing/gatk_variants/{{sample}}/{interval}.g.vcf",
            sample = config["samples"],
            interval = config["intervals"])
    output:
        "results/testing/gatk_variants/gathered_variants/{sample}.g.vcf"
    threads: 4
    resources:
        mem_mb = 20000,
        runtime = 720
    conda:
        "../envs/gatk.yaml"
    script:
        "../scripts/python/gather_gvcf.py"

rule sortvcf:
    input:
        "results/testing/gatk_variants/gathered_variants/{sample}.g.vcf"
    output:
        "results/testing/gatk_variants/sorted_variants/{sample}.sorted.g.vcf"
    conda:
        "../envs/gatk.yaml"
    threads: 1
    shell:
        "gatk SortVcf -I {input} -O {output}"

rule genomics_db_import:
    input:
        gvcf=expand("results/testing/gatk_variants/sorted_variants/{sample}.sorted.g.vcf",
            sample = config["samples"]),
        intervals="data/ref_genome/splitintervals/{interval}",
    output:
        db=directory("results/testing/variants/genomics_dbs/db_{interval}"),
    threads: 8
    resources:
        mem_mb = 20000,
        runtime = 720
    conda:
        "../envs/gatk.yaml"
    script:
        "../scripts/python/genomicsDBimport.py"


rule genotype_db:
    input:
        db = "results/testing/variants/genomics_dbs/db_{interval}",
        ref = "data/ref_genome/EquCab3.0_chrom.fa",
        intervals = "data/ref_genome/splitintervals/{interval}"
    output:
        vcf = "results/testing/variants/genomics_vcf/{interval}.vcf"
    threads: 6
    resources:
        mem_mb = 16000,
        runtime = 1440
    conda:
        "../envs/gatk.yaml"
    script:
        "../scripts/python/genotype_dbs.py"
