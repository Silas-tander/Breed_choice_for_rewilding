# breed_choice: variant filtering, population structure, relatedness, ROH and genetic load of 146 horses.
# Run: snakemake --use-conda -j <cores>. Inputs from outside the workflow live under data/ and metadata/.

import os

configfile: "config/config.yaml"

n_bootstraps = 11  # ADMIXTURE K = 1..n_bootstraps

PREFIX = "all_bcf_exclude2nd_dom2"

REL_ANALYSIS_DIR = "results/breed_choice/unrelated_theta0.125"  # PCA / ADMIXTURE / KING on the ngsRelate theta <= 0.125 keep list

GENO_VALUES = ["0.99"]  # --geno values of the PCA sensitivity runs

PLINK_DIR  = f"{REL_ANALYSIS_DIR}/plink"
PRUNE_PREF = f"{PLINK_DIR}/{PREFIX}"

rule all:
    input:
        # Chromosome-name lookup
        "results/breed_choice/chr",

        # Filtered and LD-pruned VCFs
        expand(
            "results/breed_choice/gvcf_filt/{autosome}.GTmiss.gvcf",
            autosome=config["autosomes_only"],
        ),
        expand(
            "results/breed_choice/gvcf_LD/{autosome}_GTmiss.LDprune_bcf.gvcf",
            autosome=config["autosomes_only"],
        ),

        # PCA (LD-pruned; MAF >= 0.05)
        f"{PRUNE_PREF}.pruned.pca.eigenvec",
        f"{PRUNE_PREF}.pruned.pca.eigenval",

        f"{PRUNE_PREF}.pruned.pca.maf.eigenvec",
        f"{PRUNE_PREF}.pruned.pca.maf.eigenval",

        # Missingness QC
        "results/breed_choice/genomics_stats/missingness/all_GTmiss.imiss",
        "results/breed_choice/genomics_stats/missingness/all_GTmiss.lmiss",
        "results/breed_choice/genomics_stats/missingness/all_reps.imiss",
        "results/breed_choice/genomics_stats/missingness/all_reps.lmiss",
        "results/breed_choice/genomics_stats/missingness/all_GTmiss_thin50k.traw",

        # ADMIXTURE (LD-pruned; MAF >= 0.05)
        expand(
            f"{REL_ANALYSIS_DIR}/admixture_results/{{n}}/log.out",
            n=[i + 1 for i in range(n_bootstraps)],
        ),
        expand(
            f"{REL_ANALYSIS_DIR}/admixture_results_maf/{{n}}/log.out",
            n=[i + 1 for i in range(n_bootstraps)],
        ),

        # KING on the keep list
        f"{REL_ANALYSIS_DIR}/king/king_related.kin",
        f"{REL_ANALYSIS_DIR}/king/king_kinship.kin",
        f"{REL_ANALYSIS_DIR}/king/king_ibs.seg",

        # ROH in cM
        "results/breed_choice/ROH/bcftools/AF_all/chr_FROHints_bcf_cM_AF.csv",

        # Genetic load: all sites, and without near-fixed sites (occupancy > 0.9 / 0.5)
        expand(
            "results/breed_choice/genotype_diagnostics/{autosome}_genetic_load.csv",
            autosome=config["autosomes_only"],
        ),
        expand(
            "results/breed_choice/genotype_diagnostics/{autosome}_genetic_load_filtered_0.9_HIGH+MODERATE.csv",
            autosome=config["autosomes_only"],
        ),
        expand(
            "results/breed_choice/genotype_diagnostics/{autosome}_genetic_load_filtered_0.5_HIGH+MODERATE.csv",
            autosome=config["autosomes_only"],
        ),

        # Genotype counts, metadata table, heterozygosity denominator
        expand(
            "results/breed_choice/genotype_diagnostics/{autosome}_genotype_counts.csv",
            autosome=config["autosomes_only"],
        ),
        "results/breed_choice/metadata/metadata_breed_choice_updated.csv",
        "results/breed_choice/genomics_stats/autosome_uncallable_bp.tsv",

        # Binned depth and SNP density
        "results/breed_choice/depth_binned/all_autosomes.depth_binned.tsv",
        "results/breed_choice/depth_binned/all_autosomes.snp_density_binned.tsv",

        # ngsRelate keep list
        f"{REL_ANALYSIS_DIR}/relatedness/excluded_samples_breed_choice.csv",
        f"{REL_ANALYSIS_DIR}/relatedness/related_pairs_breed_choice.csv",

        # Alternative KING-based keep lists
        "results/breed_choice/relatedness_king/unrelated_samples_king_breed_choice.txt",
        "results/breed_choice/relatedness_king/related_pairs_king.csv",
        "results/breed_choice/relatedness_king/king_mind0.20/king_ibs.seg",
        "results/breed_choice/relatedness_king/king_mind0.20/king_kinship.kin0",
        "results/breed_choice/relatedness_king/king_mind0.20/unrelated_samples_king_mind_breed_choice.txt",

        # Pairwise individual Fst
        "results/breed_choice/fst_individuals/pairwise_hudson_fst_roh_horses.csv",

        # Alternative keep list, all eigenvalues, --geno sensitivity, filtration counts
        "results/breed_choice/unrelated_ngsKING0.0884/plink/all_bcf_exclude2nd_dom2.pruned.pca.eigenvec",
        "results/breed_choice/unrelated_ngsKING0.0884/plink/all_bcf_exclude2nd_dom2.pruned.pca.maf.eigenvec",
        "results/breed_choice/unrelated_theta0.125/plink/all_bcf_exclude2nd_dom2.pruned.pca.all.eigenval",
        "results/breed_choice/unrelated_ngsKING0.0884/plink/all_bcf_exclude2nd_dom2.pruned.pca.all.eigenval",
        expand(
            "results/breed_choice/unrelated_theta0.125/plink_geno{geno}/all_bcf_exclude2nd_dom2.pruned.pca.{x}",
            geno=GENO_VALUES, x=["eigenvec", "maf.eigenvec", "all.eigenval"],
        ),
        "results/breed_choice/genomics_stats/filtration_counts.tsv",
        expand(
            "results/breed_choice/unrelated_ngsKING0.0884/admixture_results/{n}/log.out",
            n=[i + 1 for i in range(n_bootstraps)],
        ),
        expand(
            "results/breed_choice/unrelated_ngsKING0.0884/admixture_results_maf/{n}/log.out",
            n=[i + 1 for i in range(n_bootstraps)],
        ),

include: "rules/rules_breed_choice/gvcf_filter_clean.smk"
include: "rules/rules_breed_choice/popgen.smk"
include: "rules/rules_breed_choice/plink_clean_old.smk"
include: "rules/rules_breed_choice/missingness_breed_choice.smk"
include: "rules/rules_breed_choice/fastq_formatting.smk"
include: "rules/rules_breed_choice/insert_size_breed_choice.smk"
include: "rules/rules_breed_choice/metadata_breed_choice.smk"
include: "rules/rules_breed_choice/depth_binned_breed_choice.smk"
include: "rules/rules_breed_choice/relatedness_exclude_breed_choice.smk"
include: "rules/rules_breed_choice/king_all_breed_choice.smk"
include: "rules/rules_breed_choice/king_exclude_breed_choice.smk"
include: "rules/rules_breed_choice/fst_individuals_breed_choice.smk"
include: "rules/rules_breed_choice/popstruct_relset_breed_choice.smk"
include: "rules/rules_breed_choice/het_denominator_breed_choice.smk"
include: "rules/rules_breed_choice/popstruct_geno_breed_choice.smk"
include: "rules/rules_breed_choice/filtration_counts_breed_choice.smk"
