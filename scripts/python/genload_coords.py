"""BED files of variant positions per SnpEff impact and effect (synonymous, intergenic).
Usage: python genload_coords.py <chr> <outdir> <vcf>
"""

import allel
import os
import sys
import pandas as pd
import numpy as np

chr = sys.argv[1]
outdir = sys.argv[2]
vcf_path = sys.argv[3]

os.makedirs(outdir, exist_ok=True)

callset = allel.read_vcf(vcf_path, fields=['variants/POS', 'variants/ANN', 'variants/LOF'])

positions = callset['variants/POS']
annotations = callset['variants/ANN']
LOF = callset['variants/LOF']

impact_levels = ['LOW', 'MODIFIER', 'HIGH', 'MODERATE']
effects_of_interest = ["synonymous_variant", "intergenic_region"]

def parse_annotation(annotation):
    fields = annotation.split('|')
    effect = fields[1]
    impact = fields[2]
    return effect, impact

impact_bed = {impact: [] for impact in impact_levels}
effects_bed = {effect: [] for effect in effects_of_interest}

for variant_idx, variant_ann in enumerate(annotations):
    if variant_ann:
        effect, impact = parse_annotation(variant_ann)
        pos = positions[variant_idx]

        if impact in impact_levels:
            impact_bed[impact].append([chr, pos, pos])

        if effect in effects_of_interest:
            effects_bed[effect].append([chr, pos, pos])

for impact, bed_data in impact_bed.items():
    impact_bed_path = os.path.join(outdir, f"{chr}_{impact}.bed")
    pd.DataFrame(bed_data, columns=["chr", "start", "end"]).to_csv(
        impact_bed_path, sep="\t", header=False, index=False
    )

for effect, bed_data in effects_bed.items():
    effect_bed_path = os.path.join(outdir, f"{chr}_{effect}.bed")
    pd.DataFrame(bed_data, columns=["chr", "start", "end"]).to_csv(
        effect_bed_path, sep="\t", header=False, index=False
    )

print(f"BED files saved in {outdir}")
