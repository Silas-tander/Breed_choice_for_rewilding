"""Genotype heatmap (samples x sites) of one chromosome VCF.
Usage: python plot_variant_heatmap.py <chr> <outdir> <vcf>
"""

import os
import pandas as pd
import sys
import allel
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import math
import matplotlib.cm as cm

chr = sys.argv[1]
vcf_path = sys.argv[3]
callset = allel.read_vcf(vcf_path, fields=['POS', 'samples', 'calldata/GT', 'calldata/DP', 'variants/DP', 'variants/QD'])

outdir = sys.argv[2]

os.makedirs(outdir, exist_ok=True)

base_df = pd.DataFrame({
    'pos': callset['variants/POS']
})

sample_dfs = {}

for sample_index, sample in enumerate(callset['samples']):
    sample_data = base_df.copy()

    genotypes = callset['calldata/GT'][:, sample_index]

    sample_data['genotype'] = genotypes.tolist()

    sample_dfs[sample] = sample_data

GT_map = {
    (-1, -1): -1,
    (0, 0): 0,  
    (0, 1): 1,   
    (1, 1): 2
}

def map_genotype_converter(genotype):
    return GT_map.get(tuple(genotype), -1)

combined_df = pd.DataFrame()
for sample, df in sample_dfs.items():
    df = df.copy()
    df['sample'] = sample
    df['GT'] = df['genotype'].apply(map_genotype_converter)
    combined_df = pd.concat([combined_df, df], ignore_index=True)

df_pivot_all = combined_df.pivot(index='sample', columns='pos', values='GT')

cmap = {
    -1: 'red',
    0: 'grey',
    1: 'blue',
    2: 'black'
}

color_palette = sns.color_palette([cmap.get(val, 'white') for val in sorted(cmap)])

plt.figure(figsize=(16, len(sample_dfs) * 0.5))

sns.heatmap(
    df_pivot_all,
    cmap=color_palette,
    cbar=False,
    linewidths=0.0001,
    linecolor='black',
    yticklabels=True
)

plt.xticks([])

plt.xlabel('Position (pos)')
plt.ylabel('Sample')
plt.title('Heatmap of Genotype Values by Position for All Samples')

output_file = os.path.join(outdir, f'{chr}_genotype_heatmap.png')
plt.savefig(output_file, dpi=300, bbox_inches='tight')
