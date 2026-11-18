"""Distances (and log10) between consecutive heterozygous sites per sample of one VCF, saved to HDF5.
Usage: python interval_calculation.py <vcf> <outdir> <filename>
"""

import os
import pandas as pd
import sys
import allel
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

vcf_path = sys.argv[1]
save_directory = sys.argv[2]
filename = sys.argv[3]

os.makedirs(save_directory, exist_ok=True)

file_path = os.path.join(save_directory, filename)

callset = allel.read_vcf(vcf_path, fields=['POS', 'samples', 'calldata/GT'])

base_df = pd.DataFrame({
    'pos': callset['variants/POS']
})

sample_dfs = {}

for sample_index, sample in enumerate(callset['samples']):
    sample_data = base_df.copy()

    genotypes = callset['calldata/GT'][:, sample_index]

    sample_data['genotype'] = genotypes.tolist()

    sample_dfs[sample] = sample_data

def filter_for_01_genotype(df):
    return df[df['genotype'].apply(lambda x: x == [0, 1])]

filtered_sample_dfs = {sample: filter_for_01_genotype(df) for sample, df in sample_dfs.items()}

def calculate_pos_diff(df):
    df = df.sort_values(by='pos')
    pos_diff = df['pos'].diff().dropna().tolist()

    log10_diffs = np.log10([diff for diff in pos_diff if diff > 0])

    return pos_diff, log10_diffs

combined_diffs = {sample: calculate_pos_diff(df) for sample, df in filtered_sample_dfs.items()}

with pd.HDFStore(file_path) as store:
    for sample, (pos_diff, log10_diff) in combined_diffs.items():
        df = pd.DataFrame({'pos_diff': pos_diff, 'log10_diff': log10_diff})
        store.put(sample, df)
