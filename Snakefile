# Snakefile
import os
configfile: "config_files/complete_config.yaml"

rule all:
    input:
        expand("results/gwf/gvcf_filt/{autosome}.GTmiss.gvcf", autosome = config["autosomes"]),
        expand("results/gwf/gvcf_filt/{autosome}_GTmiss.LDprune.gvcf.recode.vcf", autosome = config["autosomes"]),
        expand("results/gwf/gvcf_LD/{autosome}_GTmiss.LDprune_bcf.gvcf", autosome = config["autosomes"]),


include: "rules/gvcf_filter.smk"
include: "rules/popgen.smk"
include: "rules/plink.smk"
include: "rules/SNP_calling.smk"