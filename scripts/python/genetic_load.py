"""Per-sample het and hom-alt counts by SnpEff impact (+ LOF, synonymous, intergenic) of one chromosome VCF.
Usage: python genetic_load.py <chr> <outdir> <vcf>  ->  <outdir>/<chr>_genetic_load.csv
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

os.makedirs(outdir, exist_ok=True)

callset = allel.read_vcf(vcf_path, fields=['calldata/GT', 'variants/ANN', 'samples', 'variants/LOF'])

samples = callset['samples']
genotypes = callset['calldata/GT']
annotations = callset['variants/ANN']
LOF = callset['variants/LOF']

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

for variant_idx, variant_ann in enumerate(annotations):
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

heterozygous_df = pd.DataFrame.from_dict(
    {sample: counts['heterozygous'] for sample, counts in sample_counts.items()},
    orient='index').fillna(0).astype(int)

homozygous_df = pd.DataFrame.from_dict(
    {sample: counts['homozygous'] for sample, counts in sample_counts.items()},
    orient='index').fillna(0).astype(int)

combined_df = pd.concat({'Heterozygous': heterozygous_df, 'Homozygous': homozygous_df}, axis=1)

output_file = os.path.join(outdir, f'{chr}_genetic_load.csv')

combined_df.to_csv(output_file, index=True)
