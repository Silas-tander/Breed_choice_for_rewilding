# GONE Ne, ROH (bcftools roh), SnpEff annotation and genetic load, ADMIXTURE, ROH overlap (IoU).

# GONE runs inside scripts/GONE/Linux, so inputs are copied there and outputs moved out
rule GOne:
    params:
        mapa = "data/gone/{interval}.map",
        ped = "data/gone/{interval}.ped",
        prefix = "{interval}",
        input_params = "scripts/GONE/Linux/INPUT_PARAMETERS_FILE"
    output:
        general = "scripts/GONE/chr_results_{interval}/OUTPUT_{interval}",
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 120
    shell:
        '''
        if [ ! -f scripts/GONE/Linux/$(basename {params.mapa}) ]; then
            cp {params.mapa} scripts/GONE/Linux/
        fi
        if [ ! -f scripts/GONE/Linux/$(basename {params.ped}) ]; then
            cp {params.ped} scripts/GONE/Linux/
        fi
        cd scripts/GONE/
        mkdir -p chr_results_{params.prefix}/
        cd Linux/
        ./script_GONE.sh {params.prefix}
        mv OUTPUT_{params.prefix} ../results_{params.prefix}/
        mv Output_Ne_{params.prefix} ../results_{params.prefix}/
        mv Output_d2_{params.prefix} ../results_{params.prefix}/
        mv outfileHWD ../results_{params.prefix}/
        mv seedfile ../results_{params.prefix}/
        mv timefile ../results_{params.prefix}/
        rm -r TEMPORARY_FILES/
        rm {params.prefix}.map
        rm {params.prefix}.ped
        cd - # Returns to the previous directory
        '''

rule GOne_all:
    params:
        mapa = "results/breed_choice/variants/gathered_vcf/plink/{population}.map",
        ped = "results/breed_choice/variants/gathered_vcf/plink/{population}.ped",
        prefix = "{population}",
        input_params = "scripts/GONE/Linux/INPUT_PARAMETERS_FILE"
    output:
        general = "scripts/GONE/results_{population}/OUTPUT_{population}",
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 120
    shell:
        '''
        if [ ! -f scripts/GONE/Linux/$(basename {params.mapa}) ]; then
            cp {params.mapa} scripts/GONE/Linux/
        fi
        if [ ! -f scripts/GONE/Linux/$(basename {params.ped}) ]; then
            cp {params.ped} scripts/GONE/Linux/
        fi
        cd scripts/GONE/
        mkdir -p results_{params.prefix}/
        cd Linux/
        ./script_GONE.sh {params.prefix}
        mv OUTPUT_{params.prefix} ../results_{params.prefix}/
        mv Output_Ne_{params.prefix} ../results_{params.prefix}/
        mv Output_d2_{params.prefix} ../results_{params.prefix}/
        mv outfileHWD ../results_{params.prefix}/
        mv seedfile ../results_{params.prefix}/
        mv timefile ../results_{params.prefix}/
        rm -r TEMPORARY_FILES/
        rm {params.prefix}.map
        rm {params.prefix}.ped
        cd - # Returns to the previous directory
        '''

rule GOne_bootstrap:
    params:
        mapa = "results/breed_choice/variants/gathered_vcf/plink/{population}.map",
        ped = "results/breed_choice/variants/gathered_vcf/plink/{population}.ped",
        prefix = "{population}",
        input_params = "scripts/GONE/Linux/INPUT_PARAMETERS_FILE",
        result_path = "results/breed_choice/GOne_results",
        itr = "{n}"
    output:
        general = "results/breed_choice/GOne_results/{n}/results_{population}/OUTPUT_{population}",
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 120
    shell:
        '''
        cp -R scripts/GONE/Linux/ {params.result_path}/{params.itr}/
        if [ ! -f {params.result_path}/{params.itr}/Linux/$(basename {params.mapa}) ]; then
            cp {params.mapa} {params.result_path}/{params.itr}/Linux/
        fi
        if [ ! -f {params.result_path}/{params.itr}/Linux/$(basename {params.ped}) ]; then
            cp {params.ped} {params.result_path}/{params.itr}/Linux/
        fi
        cd {params.result_path}/{params.itr}
        mkdir -p results_{params.prefix}/
        cd Linux/
        ./script_GONE.sh {params.prefix}
        mv OUTPUT_{params.prefix} ../results_{params.prefix}/
        mv Output_Ne_{params.prefix} ../results_{params.prefix}/
        mv Output_d2_{params.prefix} ../results_{params.prefix}/
        mv outfileHWD ../results_{params.prefix}/
        mv seedfile ../results_{params.prefix}/
        mv timefile ../results_{params.prefix}/
        rm -r TEMPORARY_FILES/
        rm {params.prefix}.map
        rm {params.prefix}.ped
        cd ..
        rm -r Linux/
        '''


autosome_to_eca = dict(zip(config["autosomes_only"], config["ECA"]))

# Averaged recombination map per autosome (chrN -> ECA name)
def get_recomb_map(wildcards):
    eca_chrom = autosome_to_eca[wildcards.autosome]
    return f"data/recombination_maps/Averaged_{eca_chrom}_map.txt"


# Per-autosome rename for ROH: all samples, one recombination map per chromosome
rule rename_chr_GTmiss_roh_breed_choice:
    input:
        vcf = "results/breed_choice/gvcf_filt/{autosome}.GTmiss.gvcf",
        chr_dir = "results/breed_choice/chr"  # ordering only: the directory() output is the DAG edge, not chr.names.txt
    params:
        names = "results/breed_choice/chr/chr.names.txt"
    output:
        vcf_re = "results/breed_choice/gvcf_filt/{autosome}.GTmiss.rename.gvcf",
        vcf_csi = "results/breed_choice/gvcf_filt/{autosome}.GTmiss.rename.gvcf.csi"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 120
    conda:
        "../../envs/bcftools.yml"
    log:
        "logs/breed_choice_roh/rename_chr_GTmiss_roh_breed_choice.{autosome}.log"
    shell:
        '''
        exec > {log} 2>&1
        set -euo pipefail
        bcftools annotate --rename-chrs {params.names} {input.vcf} -Oz -o {output.vcf_re}
        bcftools index --csi {output.vcf_re}
        '''

# Cohort allele frequencies (AC/AN) for bcftools roh --AF-file
rule af_for_roh:
    input:
        vcf="results/breed_choice/gvcf_filt/{autosome}.GTmiss.rename.gvcf"
    output:
        af_tab_gz="results/breed_choice/AF/{autosome}.freqs.tab.gz",
        af_tbi="results/breed_choice/AF/{autosome}.freqs.tab.gz.tbi"
    threads: 4
    resources:
        mem_mb=2*8000,
        runtime=2*180
    conda:
        "../../envs/bcftools2.yml"
    log:
        "logs/ROH_AF/{autosome}.log"
    shell:
        r'''
        set -euo pipefail
        mkdir -p results/breed_choice/AF logs/ROH_AF

        # Ensure AC/AN exist, keep biallelic SNPs, compute AF=AC/AN.
        # Output exactly: CHROM<TAB>POS<TAB>REF,ALT<TAB>AF
        bcftools +fill-tags {input.vcf} -Ou -- -t AC,AN \
        | bcftools view -m2 -M2 -v snps -Ou \
        | bcftools query -f '%CHROM\t%POS\t%REF,%ALT\t%INFO/AC\t%INFO/AN\n' \
        | awk '($5 > 0) {{ printf "%s\t%s\t%s\t%.6f\n", $1, $2, $3, $4/$5 }}' \
        | bgzip -c > {output.af_tab_gz} 2> {log}

        tabix -s 1 -b 2 -e 2 {output.af_tab_gz} 2>> {log}
        '''

# -G 30: genotypes only (PLs ignored), as recommended in the bcftools roh docs
rule ROH_bcftools:
    input:
        vcf="results/breed_choice/gvcf_filt/{autosome}.GTmiss.rename.gvcf",
        vcf_csi="results/breed_choice/gvcf_filt/{autosome}.GTmiss.rename.gvcf.csi",
        raw_recomb_map=get_recomb_map,
        af_tab_gz="results/breed_choice/AF/{autosome}.freqs.tab.gz",
        af_tbi="results/breed_choice/AF/{autosome}.freqs.tab.gz.tbi"
    output:
        roh="results/breed_choice/ROH/bcftools/AF/{autosome}.AF.roh.txt",
        fixed_recomb_map="results/breed_choice/recombination_maps/IMPUTE2/{autosome}_fixed.txt"
    threads: 3
    resources:
        mem_mb=3*8000,
        runtime=2*620
    conda:
        "../../envs/bcftools2.yml"
    log:
        "logs/ROH_bcftools/{autosome}.log"
    shell:
        r'''
        set -euo pipefail
        mkdir -p results/breed_choice/ROH/bcftools/AF logs/ROH_bcftools results/breed_choice/recombination_maps/IMPUTE2

        # Fix recombination map header/cols for bcftools roh
        awk 'NR==1 {{print "position COMBINED_rate(cM/Mb) Genetic_Map(cM)"}}
             NR>1 {{print $1, $2, $3}}' {input.raw_recomb_map} > {output.fixed_recomb_map}

        # Run bcftools roh with indexed VCF and indexed AF file
        bcftools roh \
            --AF-file {input.af_tab_gz} \
            -m {output.fixed_recomb_map} \
            -M 100 \
            -G 30 \
            -o {output.roh} \
            {input.vcf} &> {log}
        '''

rule combine_roh:
    input:
        roh_files = lambda wildcards: expand(
            "results/breed_choice/ROH/bcftools/AF/{autosome}.AF.roh.txt",
            autosome=config["autosomes_only"],
        )
    params:
        roh_dir = "results/breed_choice/ROH/bcftools/AF"
    output:
        csv = "results/breed_choice/ROH/bcftools/AF_all/all_samples.ROH_df.csv"
    log:
        "logs/combine_roh.log"
    threads: 2
    resources:
        mem_mb = 80000,
        runtime = 600
    conda:
        "../../envs/slurm_v2.yaml"
    shell:
        """
        Rscript scripts/r/process_bcftoolsROH_df_snakemake.r {params.roh_dir} {output.csv} > {log} 2>&1
        """


rule calculate_snp_intervals:
    input:
        vcf = "results/breed_choice/gvcf_filt/{autosome}.GTmiss_wReps.gvcf"
    params:
        outdir = "results/breed_choice/ROH/GTmiss_wReps/",
        filename = "{autosome}"
    output:
        outy = "results/breed_choice/ROH/GTmiss_wReps/{autosome}"
    threads: 1
    resources:
        mem_mb = 40000,
        runtime = 120
    conda:
        "../../envs/slurm_v2.yaml"
    shell:
        '''
        python scripts/python/interval_calculation.py {input.vcf} {params.outdir} {params.filename}
        '''

rule plot_genotypes:
    input:
        vcf = "results/breed_choice/gvcf_filt/{autosome}.GTmiss_wReps.gvcf"
    params:
        chr = "{autosome}",
        outdir = "results/breed_choice/genotype_heatmaps/"
    output:
        heatmap = "results/breed_choice/genotype_heatmaps/{autosome}_genotype_heatmap.png"
    threads: 2
    resources:
        mem_mb = 80000,
        runtime = 120
    conda:
        "../../envs/slurm_v2.yaml"
    shell:
        "python scripts/python/plot_variant_heatmap.py {params.chr} {params.outdir} {input.vcf}"


# SnpEff's EquCab3.0.99 database uses numeric chromosome names
rule rename_chr_wReps_breed_choice:
    input:
        vcf = "results/breed_choice/gvcf_filt/{autosome}.GTmiss_wReps.gvcf",
        chr_dir = "results/breed_choice/chr"  # ordering only: needs stage_chr_names_breed_choice
    params:
        names = "results/breed_choice/chr/chr.names.txt"
    output:
        vcf_re = "results/breed_choice/gvcf_filt/{autosome}.GTmiss_wReps.rename.gvcf"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 120
    conda:
        "../../envs/bcftools.yml"
    log:
        "logs/breed_choice_snpeff/rename_chr_wReps_breed_choice.{autosome}.log"
    shell:
        '''
        exec > {log} 2>&1
        set -euo pipefail
        bcftools annotate --rename-chrs {params.names} {input.vcf} -Oz -o {output.vcf_re}
        '''

rule snpEFF_annotation:
    input:
        calls="results/breed_choice/gvcf_filt/{autosome}.GTmiss_wReps.rename.gvcf",
    params:
        db="EquCab3.0.99",  # SnpEff database
        snpeff_path = "scripts/snpEff/",
        result_path = "results/breed_choice/snpeff/",
        chr_path = "../../results/breed_choice/snpeff/{autosome}.ann.gvcf",
        stats="{autosome}.html",
        csvstats="{autosome}.csv"
    output:
        calls="results/breed_choice/snpeff/{autosome}.ann.gvcf",
    resources:
        mem_mb=8000,
        runtime=120
    conda:
        "../../envs/snpeff.yml"
    shell:
        '''
        # if [ ! -d "{params.result_path}" ]; then
        #     mkdir -p "{params.result_path}"
        cd {params.snpeff_path}
        bcftools view ../../{input.calls} | java -Xmx8g -jar snpEff.jar -v {params.db} > {params.chr_path}
        '''

rule snpEFF_annotation_single:
    input:
        calls="data/gvcf_single/{autosome}/{sample}.bcf",
    params:
        db="EquCab3.0.99",
        snpeff_path = "scripts/snpEff/",
        result_path = "results/breed_choice/snpeff/single/{sample}/",
        chr_path = "../../results/breed_choice/snpeff/single/{sample}/{autosome}.ann.gvcf",
        stats="{autosome}.html",
        csvstats="{autosome}.csv"
    output:
        calls="results/breed_choice/snpeff/single/{sample}/{autosome}.ann.gvcf",
    resources:
        mem_mb=8000,
        runtime=120
    conda:
        "../../envs/snpeff.yml"
    shell:
        '''
        # if [ ! -d "{params.result_path}" ]; then
        #     mkdir -p "{params.result_path}"
        cd {params.snpeff_path}
        bcftools view {input.calls} | java -Xmx8g -jar snpEff.jar -v {params.db} > {params.chr_path}
        '''


rule genotype_genetic_load:
    input:
        vcf = "results/breed_choice/snpeff/{autosome}.ann.gvcf"
    params:
        chr = "{autosome}",
        outdir = "results/breed_choice/genotype_diagnostics/"
    output:
        genload = "results/breed_choice/genotype_diagnostics/{autosome}_genetic_load.csv"
    threads: 2
    resources:
        mem_mb = 20000,
        runtime = 620
    conda:
        "../../envs/slurm_v2.yaml"
    shell:
        "python scripts/python/genetic_load.py {params.chr} {params.outdir} {input.vcf}"

rule gather_split_vcf:
    input:
        vcf_ints=lambda wildcards: expand(
            "results/breed_choice/snpeff/single/{sample}/{autosome}.ann.gvcf",
            autosome=config["autosomes"],
            sample=wildcards.sample
        )
    output:
        "results/breed_choice/snpeff/gathered_vcf/{sample}.ann.vcf"
    threads: 2
    resources:
        mem_mb=12000,
        runtime=720 // 2
    conda:
        "../../envs/bcftools.yml"
    shell:
        """
        bcftools concat {input.vcf_ints} -O z -o {output}
        """

rule gather_split_vcf_GTmiss:
    input:
        vcf_ints=lambda wildcards: expand(
            "results/breed_choice/gvcf_filt/{autosome}.GTmiss.gvcf",
            autosome=config["autosomes_only"]
        )
    output:
        "results/breed_choice/gathered_vcf/all_GTmiss.gvcf"
    threads: 2
    resources:
        mem_mb=12000,
        runtime=720 // 2
    conda:
        "../../envs/bcftools.yml"
    shell:
        """
        bcftools concat {input.vcf_ints} -O z -o {output}
        """

rule exclude_2nd_degree_related:
    input:
        vcf = "results/breed_choice/gathered_vcf/all_GTmiss.gvcf",
        samples = "metadata/domestic_notinbred.txt"
    output:
        vcf = "results/breed_choice/gathered_vcf/all_bcf_exclude2nd.gvcf",
    threads: 2
    resources:
        mem_mb=12000,
        runtime=720 // 2
    conda:
        "../../envs/bcftools.yml"
    shell:
        """
        bcftools view -S {input.samples} -Ov -o {output.vcf} {input.vcf} --force-samples
        """

rule deleterious_occupancy:
    input:
        vcf = "results/breed_choice/snpeff/{autosome}.ann.gvcf"
    params:
        chr = "{autosome}",
        outdir = "results/breed_choice/genotype_diagnostics/"
    output:
        genload = "results/breed_choice/genotype_diagnostics/{autosome}_genetic_load_by_site.csv"
    threads: 2
    resources:
        mem_mb = 20000,
        runtime = 620
    conda:
        "../../envs/slurm_v2.yaml"
    shell:
        "python scripts/python/deleterious_occupancy.py {params.chr} {params.outdir} {input.vcf}"

# Genetic load without sites whose alternate-allele occupancy exceeds the threshold at the
# given impact levels (near-fixed / reference-bias sites): mild = 0.9, strict = 0.5
rule genotype_genetic_load_mild_filtering:
    input:
        vcf = "results/breed_choice/snpeff/{autosome}.ann.gvcf",
        occupancy = "results/breed_choice/genotype_diagnostics/{autosome}_genetic_load_by_site.csv"
    params:
        chr = "{autosome}",
        outdir = "results/breed_choice/genotype_diagnostics/",
        threshold = "0.9",
        impact_levels = "HIGH,MODERATE"
    output:
        genload = "results/breed_choice/genotype_diagnostics/{autosome}_genetic_load_filtered_0.9_HIGH+MODERATE.csv"
    threads: 2
    resources:
        mem_mb = 20000,
        runtime = 620
    conda:
        "../../envs/slurm_v2.yaml"
    shell:
        "python scripts/python/genetic_load_filtered.py {params.chr} {params.outdir} {input.vcf} {input.occupancy} {params.threshold} {params.impact_levels}"

rule genotype_genetic_load_strict_filtering:
    input:
        vcf = "results/breed_choice/snpeff/{autosome}.ann.gvcf",
        occupancy = "results/breed_choice/genotype_diagnostics/{autosome}_genetic_load_by_site.csv"
    params:
        chr = "{autosome}",
        outdir = "results/breed_choice/genotype_diagnostics/",
        threshold = "0.5",
        impact_levels = "HIGH,MODERATE"
    output:
        genload = "results/breed_choice/genotype_diagnostics/{autosome}_genetic_load_filtered_0.5_HIGH+MODERATE.csv"
    threads: 2
    resources:
        mem_mb = 20000,
        runtime = 620
    conda:
        "../../envs/slurm_v2.yaml"
    shell:
        "python scripts/python/genetic_load_filtered.py {params.chr} {params.outdir} {input.vcf} {input.occupancy} {params.threshold} {params.impact_levels}"

rule shared_polymorphism:
    input:
        vcf = "results/breed_choice/snpeff/{autosome}.ann.gvcf"
    params:
        chr = "{autosome}",
        outdir = "results/breed_choice/genotype_matrices/"
    output:
        shared_pol = "results/breed_choice/genotype_matrices/{autosome}_heterozygous_shared_genotypes.csv",
        hom = "results/breed_choice/genotype_matrices/{autosome}_homozygous_alt_shared_genotypes.csv",
        alt = "results/breed_choice/genotype_matrices/{autosome}_homozygous_ref_shared_genotypes.csv",
        miss = "results/breed_choice/genotype_matrices/{autosome}_missing_shared_genotypes.csv"
    threads: 2
    resources:
        mem_mb = 20000,
        runtime = 620
    conda:
        "../../envs/slurm_v2.yaml"
    shell:
        "python scripts/python/shared_polymorphism.py {params.chr} {params.outdir} {input.vcf}"

rule annotation_count:
    input:
        vcf = "results/breed_choice/snpeff/{autosome}.ann.gvcf"
    params:
        chr = "{autosome}",
        outdir = "results/breed_choice/genotype_diagnostics/"
    output:
        genload = "results/breed_choice/genotype_diagnostics/{autosome}.counts.csv"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 420
    conda:
        "../../envs/slurm_v2.yaml"
    shell:
        "python scripts/python/annotation_counter.py {params.chr} {params.outdir} {input.vcf}"

rule ROH_genetic_load:
    input:
        vcf = "results/breed_choice/snpeff/{autosome}.ann.gvcf"
    params:
        chr = "{autosome}",
        outdir = "results/breed_choice/genotype_diagnostics/"
    output:
        genload = "results/breed_choice/genotype_diagnostics/{autosome}_roh_variant_counts_10miss.csv"
    threads: 2
    resources:
        mem_mb = 20000,
        runtime = 620
    conda:
        "../../envs/slurm_v2.yaml"
    shell:
        "python scripts/python/ROH_genload.py {params.chr} {params.outdir} {input.vcf}"

rule test_annotations:
    input:
        vcf = "results/breed_choice/snpeff/{autosome}.ann.gvcf"
    params:
        chr = "{autosome}",
        outdir = "results/breed_choice/genotype_annotations/"
    output:
        genload = "results/breed_choice/genotype_annotations/{autosome}_ref_sites.csv"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 120
    conda:
        "../../envs/slurm_v2.yaml"
    shell:
        "python scripts/python/test_annotations.py {params.chr} {params.outdir} {input.vcf}"

rule generate_bed_files:
    input:
        vcf="results/breed_choice/snpeff/{autosome}.ann.gvcf"
    output:
        touch("results/breed_choice/bed/{autosome}_done.txt")
    params:
        outdir="results/breed_choice/bed",
        chrom = "{autosome}"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 120
    conda:
        "../../envs/slurm_v2.yaml"
    shell:
        """
        python scripts/python/genload_coords.py {params.chrom} {params.outdir} {input.vcf} \
        touch {output}
        """

# ADMIXTURE K = {n} (--cv=10, 1000 bootstraps) on the LD-pruned keep-list set
rule ADMIXTURE:
    input:
        bim = f"{REL_ANALYSIS_DIR}/plink/all_bcf_exclude2nd_dom2.pruned.bim",
        bed = f"{REL_ANALYSIS_DIR}/plink/all_bcf_exclude2nd_dom2.pruned.bed",
        fam = f"{REL_ANALYSIS_DIR}/plink/all_bcf_exclude2nd_dom2.pruned.fam"
    params:
        itr = "{n}",
        diro = f"{REL_ANALYSIS_DIR}/admixture_results/{{n}}/"
    output:
        f"{REL_ANALYSIS_DIR}/admixture_results/{{n}}/log.out",
    threads: 4
    resources:
        mem_mb = 8000,
        runtime = 2*720
    conda:
        "../../envs/admixture.yaml"
    shell:
        """
        mkdir -p {params.diro}
        admixture --cv=10 -B1000 {input.bed} {params.itr} -j4 > {output}
        mv all_bcf_exclude2nd_dom2.pruned.{params.itr}.Q {params.diro}
        mv all_bcf_exclude2nd_dom2.pruned.{params.itr}.P {params.diro}
        mv all_bcf_exclude2nd_dom2.pruned.{params.itr}.Q_bias {params.diro}
        mv all_bcf_exclude2nd_dom2.pruned.{params.itr}.Q_se {params.diro}
        """

# Same on the MAF >= 0.05 subset
rule ADMIXTURE_MAF:
    input:
        bim = f"{REL_ANALYSIS_DIR}/plink/all_bcf_exclude2nd_dom2.pruned.maf.bim",
        bed = f"{REL_ANALYSIS_DIR}/plink/all_bcf_exclude2nd_dom2.pruned.maf.bed",
        fam = f"{REL_ANALYSIS_DIR}/plink/all_bcf_exclude2nd_dom2.pruned.maf.fam"
    params:
        itr = "{n}",
        diro = f"{REL_ANALYSIS_DIR}/admixture_results_maf/{{n}}/"
    output:
        f"{REL_ANALYSIS_DIR}/admixture_results_maf/{{n}}/log.out",
    threads: 4
    resources:
        mem_mb = 8000,
        runtime = 2*720
    conda:
        "../../envs/admixture.yaml"
    shell:
        """
        mkdir -p {params.diro}
        admixture --cv=10 -B1000 {input.bed} {params.itr} -j4 > {output}
        mv all_bcf_exclude2nd_dom2.pruned.maf.{params.itr}.Q {params.diro}
        mv all_bcf_exclude2nd_dom2.pruned.maf.{params.itr}.P {params.diro}
        mv all_bcf_exclude2nd_dom2.pruned.maf.{params.itr}.Q_bias {params.diro}
        mv all_bcf_exclude2nd_dom2.pruned.maf.{params.itr}.Q_se {params.diro}
        """

# ROH lengths in cM from the recombination maps
rule roh_to_cm:
    input:
        roh_csv="results/breed_choice/ROH/bcftools/AF_all/all_samples.ROH_df.csv"
    output:
        out_csv="results/breed_choice/ROH/bcftools/AF_all/chr_FROHints_bcf_cM_AF.csv"
    params:
        recomb_dir="data/recombination_maps"
    threads: 2
    resources:
        mem_mb = 40000,
        runtime = 720
    conda:
        "../../envs/slurm_v3.yml"
    shell:
        """
        Rscript scripts/r/bcftools_ROH_cM_snakemake.r {input.roh_csv} {params.recomb_dir} {output.out_csv}
        """

# Pairwise ROH overlap (intersection over union) between individuals
rule IOU_matrix:
    input:
        roh_csv = "results/breed_choice/ROH/bcftools/AF_all/chr_FROHints_bcf_cM_AF.csv",
        vcf = "results/breed_choice/gathered_vcf/all_GTmiss.gvcf"
    params:
        metadata = "metadata/samples.csv",
        chr = "all",
        outdir = "results/breed_choice/ROH/IOU/",
        cM_cutoff = 5  # minimum ROH length (cM)
    output:
        iou_matrix = "results/breed_choice/ROH/IOU/all_iou_similarity_matrix_5.0.csv"
    threads: 2
    resources:
        mem_mb = 200000,
        runtime = 4*720
    conda:
        "../../envs/slurm_v3.yml"
    log:
        "logs/IOU_matrix_5.0cM.log"
    shell:
        """
        python scripts/python/IOU_matrix.py \
            {input.roh_csv} \
            {params.metadata} \
            {input.vcf} \
            {params.outdir} \
            {params.chr} \
            {params.cM_cutoff} \
            > {log} 2>&1
        """

rule join_AF_files:
    input:
        af_files=expand("results/breed_choice/AF/{autosome}.freqs.tab.gz",
                        autosome=config["autosomes_only"])
    output:
        joined_af="results/breed_choice/AF/all_autosomes.freqs.tab.gz",
        joined_af_tbi="results/breed_choice/AF/all_autosomes.freqs.tab.gz.tbi"
    threads: 4
    resources:
        mem_mb=12000,
        runtime=720
    conda:
        "../../envs/bcftools2.yml"
    shell:
        r"""
        set -euo pipefail
        # concat, drop any comment lines, sort by chrom+pos, bgzip & index
        zcat -f {input.af_files} \
          | awk 'BEGIN{{FS=OFS="\t"}} !/^#/{{print $0}}' \
          | LC_ALL=C sort -T "$(dirname {output.joined_af})" -k1,1V -k2,2n \
          | bgzip -@ {threads} -c > {output.joined_af}

        tabix -s 1 -b 2 -e 2 {output.joined_af}
        """

# Cohort AF for ngsrelate_run (not in rule all: the keep list uses data/ngsrelate/ngsrelate.res)
rule annotate_cohort_vcf_with_AF:
    input:
        vcf="results/breed_choice/gathered_vcf/all_GTmiss.gvcf",
        af="results/breed_choice/AF/all_autosomes.freqs.tab.gz",
        af_tbi="results/breed_choice/AF/all_autosomes.freqs.tab.gz.tbi"
    output:
        vcf="results/breed_choice/gathered_vcf/all_GTmiss.withAF.vcf.gz",
        tbi="results/breed_choice/gathered_vcf/all_GTmiss.withAF.vcf.gz.tbi"
    threads: 4
    resources:
        mem_mb=12000,
        runtime=720
    conda:
        "../../envs/bcftools2.yml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.vcf})"

        # 1) Strip any existing INFO/AF (works on plain .vcf/.gvcf)
        tmp_in=$(mktemp --suffix=.vcf)
        bcftools annotate -x INFO/AF -Ov -o "$tmp_in" {input.vcf}

        # 2) Add AF from the external table; inject header if AF not already defined
        if bcftools view -h "$tmp_in" | grep -q '^##INFO=<ID=AF'; then
            bcftools annotate \
              -a {input.af} \
              -c CHROM,POS,REF,ALT,INFO/AF \
              -Oz -o {output.vcf} "$tmp_in"
        else
            bcftools annotate \
              -a {input.af} \
              -c CHROM,POS,REF,ALT,INFO/AF \
              -h <(echo '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency from cohort (AC/AN)">') \
              -Oz -o {output.vcf} "$tmp_in"
        fi

        # 3) Tabix index (.tbi)
        bcftools index -t -f {output.vcf}

        # Cleanup
        rm -f "$tmp_in"
        """

OUTDIR = "results/breed_choice/ngsrelate"

rule ngsrelate_run:
    input:
        vcf = "results/breed_choice/gathered_vcf/all_GTmiss.withAF.vcf.gz"
    output:
        res = f"results/breed_choice/ngsrelate/ngsrelate.res"
    threads: 8
    resources:
        mem_mb=24000,
        runtime=720
    conda:
        "../../envs/ngsrelate.yml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{OUTDIR}"

        # Ensure index exists (harmless if already present)
        bcftools index -f {input.vcf} >/dev/null 2>&1 || true

        # Run ngsRelate using INFO/AF from the VCF
        ngsRelate -h {input.vcf} -A AF -T GT -c 1 -O {output.res}
        """

# Per-sample ROH segment QC for one chromosome (config: roh_debug)
rule roh_debug:
    input:
        vcf = "results/breed_choice/gvcf_filt/{chrom}.GTmiss.gvcf",
        roh = "results/breed_choice/ROH/bcftools/AF_all/chr_FROHints_bcf_cM_AF.csv"
    output:
        seg = "results/breed_choice/roh_debug/{chrom}_ROH_segments.csv",
        summary = "results/breed_choice/roh_debug/{chrom}_ROH_summary.csv"
    params:
        outdir = "results/breed_choice/roh_debug",
        include_missing = lambda wildcards: str(config.get("roh_debug", {}).get("include_missing", True)).lower(),
        include_het = lambda wildcards: str(config.get("roh_debug", {}).get("include_het", True)).lower(),
        per_sample = lambda wildcards: str(config.get("roh_debug", {}).get("per_sample", True)).lower(),
        max_rows_per_fig = lambda wildcards: int(config.get("roh_debug", {}).get("max_rows_per_fig", 30)),
        per_sample_files = lambda wildcards: str(config.get("roh_debug", {}).get("per_sample_files", True)).lower()
    threads: 4
    resources:
        mem_mb = 16000,
        runtime = 720
    conda:
        "../../envs/slurm_v3.yml"
    shell:
        r"""
        set -euo pipefail

        python scripts/python/ROH_seg_smk.py \
          --vcf {input.vcf} \
          --chrom {wildcards.chrom} \
          --roh {input.roh} \
          --outdir {params.outdir} \
          --include-missing {params.include_missing} \
          --include-het {params.include_het} \
          --per-sample {params.per_sample} \
          --max-rows-per-fig {params.max_rows_per_fig} \
          --per-sample-files {params.per_sample_files}

        # sanity checks for declared outputs
        test -s {output.seg}
        test -s {output.summary}
        """
