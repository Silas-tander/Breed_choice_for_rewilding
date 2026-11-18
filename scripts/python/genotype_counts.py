"""Per-sample genotype counts (missing, hom-ref, het, hom-alt) of one chromosome VCF; other calls
(multi-allelic, half-missing) only count towards total_genotypes.
Usage: python genotype_counts.py <chr> <outdir> <vcf>  ->  <outdir>/<chr>_genotype_counts.csv
"""
import os
import sys

import allel
import numpy as np
import pandas as pd

chrom, outdir, vcf_path = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(outdir, exist_ok=True)

_, samples, _, chunks = allel.iter_vcf_chunks(vcf_path, fields=["calldata/GT"], chunk_length=50_000)

n = len(samples)
total = 0
missing = np.zeros(n, dtype=np.int64)
hom_ref = np.zeros(n, dtype=np.int64)
het = np.zeros(n, dtype=np.int64)
hom_alt = np.zeros(n, dtype=np.int64)

for chunk, _, _, _ in chunks:
    gt = chunk["calldata/GT"]
    a, b = gt[..., 0], gt[..., 1]
    total += gt.shape[0]
    missing += ((a == -1) & (b == -1)).sum(axis=0)
    hom_ref += ((a == 0) & (b == 0)).sum(axis=0)
    het += (((a == 0) & (b == 1)) | ((a == 1) & (b == 0))).sum(axis=0)
    hom_alt += ((a == 1) & (b == 1)).sum(axis=0)

counts_df = pd.DataFrame(
    {
        "total_genotypes": total,
        "homozygous_alt": hom_alt,
        "heterozygous": het,
        "homozygous_ref": hom_ref,
        "missing": missing,
    },
    index=samples,
)
counts_df["homozygous_alt_fraction"] = counts_df["homozygous_alt"] / counts_df["total_genotypes"]
counts_df["heterozygous_fraction"] = counts_df["heterozygous"] / counts_df["total_genotypes"]
counts_df["homozygous_ref_fraction"] = counts_df["homozygous_ref"] / counts_df["total_genotypes"]
counts_df["missing_fraction"] = counts_df["missing"] / counts_df["total_genotypes"]

counts_df.to_csv(os.path.join(outdir, f"{chrom}_genotype_counts.csv"), index=True)
