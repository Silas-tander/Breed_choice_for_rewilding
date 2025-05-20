# Create a dictionary that maps autosomes to ECA names
autosome_to_eca = dict(zip(config["autosomes_only"], config["ECA"]))

# Function to dynamically resolve the recombination map
def get_recomb_map(wildcards):
    eca_chrom = autosome_to_eca[wildcards.autosome]  # Get ECA equivalent
    return f"data/recombination_maps/Averaged_{eca_chrom}_map.txt"

# this has to be refined, how do I run with self calculated allele frequencies?

rule ROH_bcftools:
    input:
        vcf="results/gwf/gvcf_filt/{autosome}.GTmiss.rename.gvcf",
        raw_recomb_map=get_recomb_map,  # Uses function to resolve ECA chromosome
        # AF="results/ROH/bcftools/AF/{autosome}.freqs.tab.gz"
    output:
        roh="results/ROH/bcftools/AF/{autosome}.AF.roh.txt",
        fixed_recomb_map="results/recombination_maps/IMPUTE2/{autosome}_fixed.txt"
    threads: 3
    resources:
        mem_mb=3*8000,
        runtime=2*620
    params:
        samples_file="metadata/samples.tsv",  # Optional: file with sample names
        transition_AZ=0.0001,  # HW to AZ transition probability (-a)
        transition_HW=0.0001,  # AZ to HW transition probability (-H)
        min_markers=10000,  # Minimum number of markers in a ROH region (-G)
    conda:
        "../envs/bcftools.yml"
    shell:
        '''
        # Convert recombination map format to match bcftools roh expected input
        awk 'NR==1 {{print "position COMBINED_rate(cM/Mb) Genetic_Map(cM)"}}
            NR>1 {{print $1, $2, $3}}' {input.raw_recomb_map} > {output.fixed_recomb_map}

        # Run bcftools roh with the corrected recombination map
        bcftools view {input.vcf} | \
        bcftools roh \
            -m {output.fixed_recomb_map} \
            -M 100 \
            -G 30 \
            -e "GT-" \
            -o {output.roh}
        '''

rule plot_genotypes:
    input:
        vcf = "results/gwf/gvcf_filt/{autosome}.GTmiss_wReps.gvcf"
    params:
        chr = "{autosome}",
        outdir = "results/genotype_heatmaps/"
    output:
        heatmap = "results/genotype_heatmaps/{autosome}_genotype_heatmap.png"
    threads: 2
    resources:
        mem_mb = 80000,
        runtime = 120
    conda:
        "../envs/slurm_v2.yaml"
    shell:
        "python scripts/python/plot_variant_heatmap.py {params.chr} {params.outdir} {input.vcf}"

#EquCab3.0.99.genome : Equus_caballus
#EquCab3.0.99.reference : ftp://ftp.ensembl.org/pub/release-99/gtf/
#EquCab3.0.99.retrieval_date : 2020-01-28

###
# it does this automatically so don't bother
###

rule snpEFF_annotation:
    input:
        calls="/home/silastander/Coregonus/Silas/Horses/results/gwf/gvcf_filt/{autosome}.GTmiss_wReps.rename.gvcf", # (vcf, bcf, or vcf.gz)
    params:
        db="EquCab3.0.99", # path to reference db downloaded with the snpeff download wrapper
        snpeff_path = "scripts/snpEff/",
        result_path = "results/snpeff/",
        chr_path = "../../results/snpeff/{autosome}.ann.gvcf",
        stats="{autosome}.html",  # summary statistics (in HTML), optional
        csvstats="{autosome}.csv" # summary statistics in CSV, optional
    output:
        calls="/home/silastander/Coregonus/Silas/Horses/results/snpeff/{autosome}.ann.gvcf",   # annotated calls (vcf, bcf, or vcf.gz)
    # log:
    #     "logs/snpeff/{autosome}.log"
    resources:
        mem_mb=8000,
        runtime=120
    conda:
        "../envs/snpeff.yml"
    shell:
        '''
        # if [ ! -d "{params.result_path}" ]; then
        #     mkdir -p "{params.result_path}"
        cd {params.snpeff_path}
        bcftools view {input.calls} | java -Xmx8g -jar snpEff.jar -v {params.db} > {params.chr_path}
        '''

rule snpEFF_annotation_single:
    input:
        calls="/home/silastander/Coregonus/Silas/Horses/results/gwf/gVCF_single/{autosome}/{sample}.bcf", # (vcf, bcf, or vcf.gz)
    params:
        db="EquCab3.0.99", # path to reference db downloaded with the snpeff download wrapper
        snpeff_path = "scripts/snpEff/",
        result_path = "results/snpeff/single/{sample}/",
        chr_path = "../../results/snpeff/single/{sample}/{autosome}.ann.gvcf",
        stats="{autosome}.html",  # summary statistics (in HTML), optional
        csvstats="{autosome}.csv" # summary statistics in CSV, optional
    output:
        calls="/home/silastander/Coregonus/Silas/Horses/results/snpeff/single/{sample}/{autosome}.ann.gvcf",   # annotated calls (vcf, bcf, or vcf.gz)
    # log:
    #     "logs/snpeff/{autosome}.log"
    resources:
        mem_mb=8000,
        runtime=120
    conda:
        "../envs/snpeff.yml"
    shell:
        '''
        # if [ ! -d "{params.result_path}" ]; then
        #     mkdir -p "{params.result_path}"
        cd {params.snpeff_path}
        bcftools view {input.calls} | java -Xmx8g -jar snpEff.jar -v {params.db} > {params.chr_path}
        '''

rule genotype_genetic_load:
    input:
        vcf = "results/snpeff/{autosome}.ann.gvcf"
    params:
        chr = "{autosome}",
        outdir = "results/genotype_diagnostics/"
    output:
        genload = "results/genotype_diagnostics/{autosome}_genetic_load.csv"
    threads: 2
    resources:
        mem_mb = 20000,
        runtime = 620
    conda:
        "../envs/slurm_v2.yaml"
    shell:
        "python scripts/python/genetic_load.py {params.chr} {params.outdir} {input.vcf}"

# merge for different samples, concat for same samples concat should be the answer here
rule gather_split_vcf:
    input:
        vcf_ints=lambda wildcards: expand(
            "results/snpeff/single/{sample}/{autosome}.ann.gvcf",
            autosome=config["autosomes"],
            sample=wildcards.sample
        )
    output:
        "results/snpeff/gathered_vcf/{sample}.ann.vcf"
    threads: 2
    resources:
        mem_mb=12000,
        runtime=720 // 2
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        bcftools concat {input.vcf_ints} -O z -o {output}
        """

# merge for different samples, concat for same samples concat should be the answer here
rule gather_split_vcf_GTmiss:
    input:
        vcf_ints=lambda wildcards: expand(
            "results/gwf/gvcf_filt/{autosome}.GTmiss.gvcf",
            autosome=config["autosomes"]
        )
    output:
        "results/gwf/gathered_vcf/all_GTmiss.gvcf"
    threads: 2
    resources:
        mem_mb=12000,
        runtime=720 // 2
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        bcftools concat {input.vcf_ints} -O z -o {output}
        """

rule genotype_genetic_load_ratios:
    input:
        vcf = "results/snpeff/{autosome}.ann.gvcf"
    params:
        chr = "{autosome}",
        outdir = "results/genotype_diagnostics/"
    output:
        genload = "results/genotype_diagnostics/{autosome}_genetic_load_ratios.csv"
    threads: 2
    resources:
        mem_mb = 20000,
        runtime = 620
    conda:
        "../envs/slurm_v2.yaml"
    shell:
        "python scripts/python/genetic_load.py {params.chr} {params.outdir} {input.vcf}"

rule deleterious_occupancy:
    input:
        vcf = "results/snpeff/{autosome}.ann.gvcf"
    params:
        chr = "{autosome}",
        outdir = "results/genotype_diagnostics/"
    output:
        genload = "results/genotype_diagnostics/{autosome}_genetic_load_by_site.csv"
    threads: 2
    resources:
        mem_mb = 20000,
        runtime = 620
    conda:
        "../envs/slurm_v2.yaml"
    shell:
        "python scripts/python/deleterious_occupancy.py {params.chr} {params.outdir} {input.vcf}"

rule annotation_count:
    input:
        vcf = "results/snpeff/{autosome}.ann.gvcf"
    params:
        chr = "{autosome}",
        outdir = "results/genotype_diagnostics/"
    output:
        genload = "results/genotype_diagnostics/{autosome}.counts.csv"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 420
    conda:
        "../envs/slurm_v2.yaml"
    shell:
        "python scripts/python/annotation_counter.py {params.chr} {params.outdir} {input.vcf}"

rule ROH_genetic_load:
    input:
        vcf = "results/snpeff/{autosome}.ann.gvcf"
    params:
        chr = "{autosome}",
        outdir = "results/genotype_diagnostics/"
    output:
        genload = "results/genotype_diagnostics/{autosome}_roh_variant_counts_10miss.csv"
    threads: 2
    resources:
        mem_mb = 20000,
        runtime = 620
    conda:
        "../envs/slurm_v2.yaml"
    shell:
        "python scripts/python/ROH_genload.py {params.chr} {params.outdir} {input.vcf}"

rule test_annotations:
    input:
        vcf = "results/snpeff/{autosome}.ann.gvcf"
    params:
        chr = "{autosome}",
        outdir = "results/genotype_annotations/"
    output:
        genload = "results/genotype_annotations/{autosome}_ref_sites.csv"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 120
    conda:
        "../envs/slurm_v2.yaml"
    shell:
        "python scripts/python/test_annotations.py {params.chr} {params.outdir} {input.vcf}"

rule generate_bed_files:
    input:
        vcf="results/snpeff/{autosome}.ann.gvcf"
    output:
        touch("results/bed/{autosome}_done.txt")
    params:
        outdir="results/bed",
        chrom = "{autosome}"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 120
    conda:
        "../envs/slurm_v2.yaml"
    shell:
        """
        python scripts/python/genload_coords.py {params.chrom} {params.outdir} {input.vcf} \
        touch {output}
        """

rule ADMIXTURE:
    params:
        bim = "results/gwf/gathered_vcf/plink/LD_pruned.bim",
        bed = "results/gwf/gathered_vcf/plink/LD_pruned.bed",
        itr = "{n}",
        diro = "results/admixture_results/{n}/"
    output:
        "results/admixture_results/{n}/log.out",
    threads: 4
    resources:
        mem_mb = 8000,
        runtime = 2*720
    conda:
        "../envs/admixture.yaml"
    shell:
        """
        mkdir -p {params.diro}
        # cd {params.diro}
        admixture --cv=10 -B1000 {params.bed} {params.itr} -j4 > {output}
        mv LD_pruned.{params.itr}.Q {params.diro}
        mv LD_pruned.{params.itr}.P {params.diro}
        mv LD_pruned.{params.itr}.Q_bias {params.diro}
        mv LD_pruned.{params.itr}.Q_se {params.diro}
        """