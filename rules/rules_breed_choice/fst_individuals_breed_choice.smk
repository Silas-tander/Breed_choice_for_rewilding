# Pairwise Hudson Fst between individual horses (2 haplotypes each) on the LD-pruned VCFs.

FST_IND_DIR = "results/breed_choice/fst_individuals"

# Horses of the ROH relatedness heatmap
FST_IND_SAMPLES = [
    "DART55", "IS079", "WP006", "FORT3344", "FORT3347", "MUST1096",
    "modoc3", "modoc4", "modoc5", "STCZ3328", "STCZ3332",
    "mols1", "mols2", "mols3", "mols4", "mols5",
    "UoN005", "UoN006", "UoN008", "UoN010", "UoN013",
    "NGSHORSE003", "NGSHORSE016", "TG1111D2630", "TG1111D2631",
    "Yak1", "Yak2", "Yak3", "Yak4", "Yak5", "Yak6", "Yak8", "Yak9",
]


rule fst_individuals_roh_breed_choice:
    input:
        vcfs = expand(
            "results/breed_choice/gvcf_LD/{autosome}_GTmiss.LDprune_bcf.gvcf",
            autosome=config["autosomes_only"],
        )
    output:
        long = f"{FST_IND_DIR}/pairwise_hudson_fst_roh_horses.csv",
        matrix = f"{FST_IND_DIR}/pairwise_hudson_fst_roh_horses_matrix.csv"
    params:
        samples = FST_IND_SAMPLES
    log:
        "results/breed_choice/logs/fst_individuals/fst_individuals_roh.log"
    resources:
        mem_mb = 16000,
        runtime = 120
    conda:
        "../../envs/allel.yml"
    script:
        "../../scripts/python/pairwise_individual_hudson_fst.py"
