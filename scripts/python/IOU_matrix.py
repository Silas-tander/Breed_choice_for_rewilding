"""Pairwise ROH overlap between samples: intersection over union of the SNPs in their ROH (>= cM_cutoff).
Usage: python IOU_matrix.py <roh.csv> <metadata.csv> <vcf> <outdir> <chr|all> <cM_cutoff>
"""

import os
import pandas as pd
import sys
import allel
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import re
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import average
from ete3 import Tree
from itertools import combinations
from itertools import combinations_with_replacement
import pyranges as pr
from tqdm import tqdm


roh_path = sys.argv[1]
metadata_path = sys.argv[2]
file_path = sys.argv[3]
outdir = sys.argv[4]
chr = sys.argv[5]
cM_cutoff = float(sys.argv[6])

chr_combined_bcf = pd.read_csv(roh_path)
chr_combined_bcf = chr_combined_bcf[chr_combined_bcf["cM"] > cM_cutoff]

meta = pd.read_csv(metadata_path)

callset = allel.read_vcf(file_path, fields=['POS', 'samples', 'calldata/GT','is_snp'])

samples = callset['samples']
genotypes = callset['calldata/GT']
genotypes_flat = genotypes[..., 0] + genotypes[..., 1]
pos = callset['variants/POS']
is_snp = callset['variants/is_snp']
gt = allel.GenotypeArray(callset['calldata/GT'])[is_snp]
snp_pos = callset['variants/POS'][is_snp]
ac = gt.count_alleles()
alt = gt.to_n_alt()

if chr == 'all':
    roh_chr = chr_combined_bcf.copy()
else:
    roh_chr = chr_combined_bcf[chr_combined_bcf['chrom'] == chr].copy()


roh_chr['snps_list'] = roh_chr.apply(
    lambda row: pos[(pos >= row['from']) & (pos <= row['to'])].tolist(),
    axis=1
)

def get_alt_values_for_sample(row, snp_pos, alt, samples):
    idx = np.where((snp_pos >= row['from']) & (snp_pos <= row['to']))[0]
    if len(idx) == 0:
        return []
    try:
        sample_idx = np.where(samples == row['sample'])[0][0]
    except IndexError:
        return []
    alt_subset = alt[idx, sample_idx]
    return alt_subset.tolist()

roh_chr['alt_values'] = roh_chr.apply(
    lambda row: get_alt_values_for_sample(row, snp_pos, alt, samples),
    axis=1
)

sample_idx = {s: i for i, s in enumerate(samples)}

chromosomes = roh_chr['chrom'].unique()
sample_list = roh_chr['sample'].unique()

intersect_sum = pd.DataFrame(0, index=sample_list, columns=sample_list, dtype=int)
union_sum = pd.DataFrame(0, index=sample_list, columns=sample_list, dtype=int)

tqdm_chromosomes = tqdm(chromosomes, desc="Chromosomes")
tqdm_sample_pairs = lambda chrom: tqdm(combinations_with_replacement(sample_list, 2), desc=f"Pairs in chr{chrom}", leave=False)
tqdm_samples = lambda chrom: tqdm(chrom_df.groupby('sample'), desc=f"Samples in chr{chrom}", leave=False)
tqdm_snps = lambda chrom: tqdm(sample_ranges.items(), desc=f"SNP+ALT in chr{chrom}", leave=False)

for chrom in tqdm_chromosomes:
    chrom_df = roh_chr[roh_chr['chrom'] == chrom]

    sample_ranges = {}
    for sample, df in tqdm_samples(chrom):
        pr_sample = pr.PyRanges(chromosomes=[chrom]*len(df), starts=df['from'], ends=df['to'])
        pr_sample = pr_sample.merge()
        sample_ranges[sample] = pr_sample

    sample_snp_alt = {}
    for sample, prange in tqdm_snps(chrom):
        idx_in_roh = []
        for interval in prange.as_df().itertuples():
            idx = np.where((snp_pos >= interval.Start) & (snp_pos <= interval.End))[0]
            idx_in_roh.extend(idx)
        idx_in_roh = np.unique(idx_in_roh)
        if len(idx_in_roh) == 0:
            continue
        alt_in_roh = alt[idx_in_roh, sample_idx[sample]]
        snp_in_roh = snp_pos[idx_in_roh]
        sample_snp_alt[sample] = dict(zip(snp_in_roh, alt_in_roh))

    for s1, s2 in tqdm_sample_pairs(chrom):
        if s1 not in sample_snp_alt or s2 not in sample_snp_alt:
            continue

        snps1 = set(sample_snp_alt[s1].keys())
        snps2 = set(sample_snp_alt[s2].keys())

        overlap_snps = snps1 & snps2
        union_snps = snps1 | snps2

        if len(union_snps) == 0:
            continue

        match_count = sum(
            sample_snp_alt[s1][pos] == sample_snp_alt[s2][pos]
            for pos in overlap_snps
        )

        intersect_sum.at[s1, s2] += match_count
        union_sum.at[s1, s2] += len(union_snps)

        intersect_sum.at[s2, s1] = intersect_sum.at[s1, s2]
        union_sum.at[s2, s1] = union_sum.at[s1, s2]
iou_matrix = intersect_sum / union_sum

iou_matrix.to_csv(os.path.join(outdir, f"{chr}_iou_similarity_matrix_{cM_cutoff}.csv"))
