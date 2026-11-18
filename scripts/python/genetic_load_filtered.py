"""genetic_load.py without sites whose alternate-allele occupancy (hom_alt + het/2) / (ref + hom_alt + het/2)
exceeds <threshold> at one of <impacts>.
Usage: python genetic_load_filtered.py <chr> <outdir> <vcf> <by_site.csv> <threshold> <impacts, e.g. HIGH,MODERATE>
"""

import allel
import os
import sys
import pandas as pd
import numpy as np
from collections import Counter, defaultdict
from itertools import combinations

chr = sys.argv[1]
outdir = sys.argv[2]
vcf_path = sys.argv[3]
occupancy_csv = sys.argv[4]
threshold_str = sys.argv[5]
threshold = float(threshold_str)
impact_levels_str = sys.argv[6]

if not (0.0 <= threshold <= 1.0):
    sys.exit(f"threshold must be between 0 and 1, got {threshold}")

VALID_IMPACTS = {"LOW", "MODIFIER", "HIGH", "MODERATE", "LOF"}
filter_impacts = [tok.strip() for tok in impact_levels_str.split(",") if tok.strip()]
invalid = set(filter_impacts) - VALID_IMPACTS
if invalid:
    sys.exit(f"invalid impact level(s) {sorted(invalid)}, must be a subset of {sorted(VALID_IMPACTS)}")
if not filter_impacts:
    sys.exit("at least one impact level must be given")

os.makedirs(outdir, exist_ok=True)

occ_df = pd.read_csv(occupancy_csv)
for col in ("homozygous_alternate", "heterozygous", "reference"):
    if col not in occ_df.columns:
        occ_df[col] = 0

hom_alt = occ_df["homozygous_alternate"].astype(float)
het = occ_df["heterozygous"].astype(float)
ref = occ_df["reference"].astype(float)
denom = ref + hom_alt + het / 2
occ_df["deleterious_occupancy"] = np.where(denom > 0, (hom_alt + het / 2) / denom, np.nan)

# A POS can have several Impact rows; only rows in filter_impacts count
excluded_positions = set(
    occ_df.loc[
        occ_df["Impact"].isin(filter_impacts) & (occ_df["deleterious_occupancy"] > threshold),
        "POS",
    ].tolist()
)

callset = allel.read_vcf(vcf_path, fields=['calldata/GT', 'variants/ANN', 'samples', 'variants/LOF', 'variants/POS'])

samples = callset['samples']
genotypes = callset['calldata/GT']
annotations = callset['variants/ANN']
LOF = callset['variants/LOF']
pos = callset['variants/POS']

impact_levels = ['LOW', 'MODIFIER', 'HIGH', 'MODERATE']
effects_of_interest = ["SYNONYMOUS_CODING", "INTERGENIC"]

sample_counts = {
    sample: {
        'heterozygous': Counter(),
        'homozygous': Counter(),
    } for sample in samples
}

def parse_annotation(annotation):
    fields = annotation.split('|')
    effect = fields[1]
    impact = fields[2]
    return effect, impact

impact_levels = ['LOW', 'MODIFIER', 'HIGH', 'MODERATE']
effects_of_interest = ["synonymous_variant", "intergenic_region"]

n_excluded = 0
for variant_idx, variant_ann in enumerate(annotations):
    if pos[variant_idx] in excluded_positions:
        n_excluded += 1
        continue
    if variant_ann:
        effect, impact = parse_annotation(variant_ann)
        for sample_idx, genotype in enumerate(genotypes[variant_idx]):
            sample = samples[sample_idx]
            if np.array_equal(genotype, [0, 1]) or np.array_equal(genotype, [1, 0]):
                sample_counts[sample]['heterozygous'][impact] += 1
                if LOF[variant_idx]:
                    sample_counts[sample]['heterozygous']['LOF'] += 1
                if effect in effects_of_interest:
                    sample_counts[sample]['heterozygous'][effect] += 1

            elif np.array_equal(genotype, [1, 1]):
                sample_counts[sample]['homozygous'][impact] += 1
                if LOF[variant_idx]:
                    sample_counts[sample]['homozygous']['LOF'] += 1
                if effect in effects_of_interest:
                    sample_counts[sample]['homozygous'][effect] += 1

print(f"{chr}: excluded {n_excluded} of {len(annotations)} sites (deleterious_occupancy > {threshold} at impact levels {filter_impacts})", file=sys.stderr)

heterozygous_df = pd.DataFrame.from_dict(
    {sample: counts['heterozygous'] for sample, counts in sample_counts.items()},
    orient='index').fillna(0).astype(int)

homozygous_df = pd.DataFrame.from_dict(
    {sample: counts['homozygous'] for sample, counts in sample_counts.items()},
    orient='index').fillna(0).astype(int)

combined_df = pd.concat({'Heterozygous': heterozygous_df, 'Homozygous': homozygous_df}, axis=1)

impact_tag = "+".join(filter_impacts)
output_file = os.path.join(outdir, f'{chr}_genetic_load_filtered_{threshold_str}_{impact_tag}.csv')

combined_df.to_csv(output_file, index=True)
