"""Pairwise Hudson Fst between individuals (2 haplotypes each), summed over the input VCFs:
Fst = sum(num) / sum(den) over sites called in both (scikit-allel); not clamped at 0. Snakemake script.
"""
from itertools import combinations

import allel
import numpy as np
import pandas as pd

samples = list(snakemake.params.samples)
pairs = list(combinations(range(len(samples)), 2))
num = np.zeros(len(pairs))
den = np.zeros(len(pairs))
n_sites = np.zeros(len(pairs), dtype=int)
n_total = 0

for vcf in snakemake.input.vcfs:
    cs = allel.read_vcf(vcf, samples=samples, fields=["samples", "calldata/GT", "variants/ALT"])
    if cs is None:
        continue
    missing = set(samples) - set(cs["samples"])
    if missing:
        raise ValueError(f"{vcf}: samples not found: {sorted(missing)}")
    order = [list(cs["samples"]).index(s) for s in samples]
    gt = allel.GenotypeArray(cs["calldata/GT"])[:, order]
    biallelic = (cs["variants/ALT"][:, 1] == "") if cs["variants/ALT"].ndim > 1 else np.ones(gt.shape[0], bool)
    gt = gt[biallelic]
    n_total += gt.shape[0]
    ac = np.stack([gt[:, [i]].count_alleles(max_allele=1) for i in range(len(samples))], axis=1)  # (variants, samples, 2)
    called = ac.sum(axis=2) == 2
    for k, (i, j) in enumerate(pairs):
        ok = called[:, i] & called[:, j]
        if not ok.any():
            continue
        a, b = allel.hudson_fst(ac[ok, i], ac[ok, j])
        num[k] += np.nansum(a)
        den[k] += np.nansum(b)
        n_sites[k] += ok.sum()

long = pd.DataFrame({
    "sample1": [samples[i] for i, _ in pairs],
    "sample2": [samples[j] for _, j in pairs],
    "fst": np.where(den > 0, num / np.where(den > 0, den, 1), np.nan),
    "n_sites": n_sites,
})
long.to_csv(snakemake.output.long, index=False)

mat = pd.DataFrame(0.0, index=samples, columns=samples)
for r in long.itertuples():
    mat.loc[r.sample1, r.sample2] = mat.loc[r.sample2, r.sample1] = r.fst
mat.to_csv(snakemake.output.matrix)

with open(snakemake.log[0], "w") as fh:
    fh.write(
        f"{len(samples)} samples, {len(pairs)} pairs, {len(snakemake.input.vcfs)} VCFs, "
        f"{n_total} biallelic sites read; sites used per pair: "
        f"min {n_sites.min()}, median {int(np.median(n_sites))}, max {n_sites.max()}; "
        f"Fst range {np.nanmin(long['fst']):.4f} to {np.nanmax(long['fst']):.4f}.\n"
    )
