"""Add per-sample genome statistics to a metadata table (left join on sample_alias):
  heterozygosity = het / autosome length (fai minus --x-contig), missing_frac = missing / total;
  --imiss: missing_frac from a PLINK .imiss instead (+ missing_GTmiss, sites_GTmiss);
  --counts-suffix S: counts-based missing columns renamed missing_S / missing_frac_S;
  --roh: sum_FROH and sum_FROH_cM (sum cM / --genome-cm) of ROH longer than --roh-min-cm.
"""
import argparse
import sys

import pandas as pd

COUNT_COLS = ["total_genotypes", "heterozygous", "homozygous_ref", "homozygous_alt", "missing"]

p = argparse.ArgumentParser()
p.add_argument("--meta", required=True)
p.add_argument("--counts", required=True, nargs="+")
p.add_argument("--fai", required=True)
p.add_argument("--x-contig", default="NC_009175.3")
p.add_argument("--counts-suffix", default=None)
p.add_argument("--imiss", default=None)
p.add_argument("--roh", default=None)
p.add_argument("--roh-min-cm", type=float, default=1.0)
p.add_argument("--genome-cm", type=float, default=2794.81 - 192.04)
p.add_argument("--out", required=True)
args = p.parse_args()

meta = pd.read_csv(args.meta, dtype={"sample_alias": str})

counts = pd.concat(
    [pd.read_csv(f, index_col=0, dtype={0: str}) for f in args.counts]
)
counts.index = counts.index.astype(str)
counts = counts.groupby(level=0)[COUNT_COLS].sum()

fai = pd.read_csv(args.fai, sep="\t", header=None, usecols=[0, 1], names=["contig", "length"])
if args.x_contig not in set(fai["contig"]):
    sys.exit(f"X contig {args.x_contig} not found in {args.fai}")
autosome_len = fai.loc[fai["contig"] != args.x_contig, "length"].sum()

counts["heterozygosity"] = counts["heterozygous"] / autosome_len
counts["missing_frac"] = counts["missing"] / counts["total_genotypes"]
if args.counts_suffix:
    s = args.counts_suffix
    counts = counts.rename(columns={"missing": f"missing_{s}", "missing_frac": f"missing_frac_{s}"})
added = counts

if args.imiss:
    imiss = pd.read_csv(args.imiss, sep=r"\s+", dtype={"IID": str}).set_index("IID")
    miss = pd.DataFrame({"missing_GTmiss": imiss["N_MISS"], "sites_GTmiss": imiss["N_GENO"]})
    miss["missing_frac"] = miss["missing_GTmiss"] / miss["sites_GTmiss"]
    added = added.drop(columns="missing_frac", errors="ignore").join(miss, how="outer")

if args.roh:
    roh = pd.read_csv(args.roh, dtype={"ind": str})
    roh = roh[roh["cM"] > args.roh_min_cm]
    froh = roh.groupby("ind").agg(sum_FROH=("FROH", "sum"), sum_FROH_cM=("cM", "sum"))
    froh["sum_FROH_cM"] = froh["sum_FROH_cM"] / args.genome_cm
    added = added.join(froh, how="outer")

added.index.name = "sample_alias"
meta = meta.drop(columns=[c for c in added.columns if c in meta.columns])
out = meta.merge(added.reset_index(), on="sample_alias", how="left")
out.to_csv(args.out, index=False)

in_meta = set(meta["sample_alias"])
print(f"autosome length (fai minus {args.x_contig}): {autosome_len}", file=sys.stderr)
print(f"metadata rows: {len(meta)}; with genotype counts: {out['total_genotypes'].notna().sum()}", file=sys.stderr)
if args.roh:
    print(f"with FROH: {out['sum_FROH'].notna().sum()}", file=sys.stderr)
not_in_meta = sorted(set(added.index) - in_meta)
if not_in_meta:
    print(f"{len(not_in_meta)} samples in counts/ROH but not in metadata (dropped): {not_in_meta[:20]}", file=sys.stderr)
