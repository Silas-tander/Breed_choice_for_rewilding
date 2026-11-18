"""Pairwise counts of shared genotypes (missing, hom-ref, het, hom-alt) between samples of one chromosome VCF.
Usage: python shared_polymorphism.py <chr> <outdir> <vcf>
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

n_samples = len(samples)

heterozygous_shared = pd.DataFrame(0, index=samples, columns=samples)
homozygous_alt_shared = pd.DataFrame(0, index=samples, columns=samples)
homozygous_ref_shared = pd.DataFrame(0, index=samples, columns=samples)
missing_shared = pd.DataFrame(0, index=samples, columns=samples)

for variant_idx, variant_gt in enumerate(genotypes):
    for (i, j) in combinations(range(n_samples), 2):
        sample_i = samples[i]
        sample_j = samples[j]
        genotype_i = variant_gt[i]
        genotype_j = variant_gt[j]

        if np.array_equal(genotype_i, [-1, -1]) and np.array_equal(genotype_j, [-1, -1]):
            missing_shared.loc[sample_i, sample_j] += 1
            missing_shared.loc[sample_j, sample_i] += 1
        elif np.array_equal(genotype_i, [0, 0]) and np.array_equal(genotype_j, [0, 0]):
            homozygous_ref_shared.loc[sample_i, sample_j] += 1
            homozygous_ref_shared.loc[sample_j, sample_i] += 1
        elif (np.array_equal(genotype_i, [0, 1]) or np.array_equal(genotype_i, [1, 0])) and \
            (np.array_equal(genotype_j, [0, 1]) or np.array_equal(genotype_j, [1, 0])):
            heterozygous_shared.loc[sample_i, sample_j] += 1
            heterozygous_shared.loc[sample_j, sample_i] += 1
        elif np.array_equal(genotype_i, [1, 1]) and np.array_equal(genotype_j, [1, 1]):
            homozygous_alt_shared.loc[sample_i, sample_j] += 1
            homozygous_alt_shared.loc[sample_j, sample_i] += 1

output_file = os.path.join(outdir, f'{chr}_heterozygous_shared_genotypes.csv')
output_file2 = os.path.join(outdir, f'{chr}_homozygous_alt_shared_genotypes.csv')
output_file3 = os.path.join(outdir, f'{chr}_homozygous_ref_shared_genotypes.csv')
output_file4 = os.path.join(outdir, f'{chr}_missing_shared_genotypes.csv')

heterozygous_shared.to_csv(output_file)
homozygous_alt_shared.to_csv(output_file2)
homozygous_ref_shared.to_csv(output_file3)
missing_shared.to_csv(output_file4)
