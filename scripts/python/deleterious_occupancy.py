"""Per-site genotype counts (het, hom-alt, hom-ref) by SnpEff impact of one chromosome VCF.
Usage: python deleterious_occupancy.py <chr> <outdir> <vcf>  ->  <outdir>/<chr>_genetic_load_by_site.csv
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

callset = allel.read_vcf(vcf_path, fields=['calldata/GT', 'variants/ANN', 'samples', 'variants/LOF', 'POS'])

samples = callset['samples']
genotypes = callset['calldata/GT']
annotations = callset['variants/ANN']
pos = callset['variants/POS']
LOF = callset['variants/LOF']

impact_levels = ['LOW', 'HIGH', 'MODERATE', 'MODIFIER']

def parse_annotation(annotation):
    fields = annotation.split('|')
    effect = fields[1]
    impact = fields[2]
    return effect, impact


site_counts = defaultdict(lambda: defaultdict(Counter))

for variant_idx, variant_ann in enumerate(annotations):
    # only the first ANN entry is used
    effect, impact = parse_annotation(variant_ann)

    if impact in impact_levels:
        site_pos = pos[variant_idx]

        for sample_idx, genotype in enumerate(genotypes[variant_idx]):
            if np.array_equal(genotype, [0, 1]) or np.array_equal(genotype, [1, 0]):
                site_counts[site_pos][impact]['heterozygous'] += 1
                if LOF[variant_idx] != '':
                    site_counts[site_pos]['LOF']['heterozygous'] += 1

            elif np.array_equal(genotype, [1, 1]):
                site_counts[site_pos][impact]['homozygous_alternate'] += 1
                if LOF[variant_idx] != '':
                    site_counts[site_pos]['LOF']['homozygous_alternate'] += 1

            elif np.array_equal(genotype, [0, 0]):
                site_counts[site_pos][impact]['reference'] += 1

rows = []
for site_pos, impacts in site_counts.items():
    for impact, genotypes in impacts.items():
        row = {'POS': site_pos, 'Impact': impact}
        row.update(genotypes)
        rows.append(row)

site_df = pd.DataFrame(rows).fillna(0)

output_file = os.path.join(outdir, f'{chr}_genetic_load_by_site.csv')
site_df.to_csv(output_file, index=False)
