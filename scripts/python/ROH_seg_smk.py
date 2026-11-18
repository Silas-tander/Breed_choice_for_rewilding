#!/usr/bin/env python3
"""Genotypes inside bcftools ROH segments of one chromosome: debug plots and segment/summary tables.
Run by rule roh_debug (popgen.smk); options: --help.
"""

import argparse
import json
import os
import math
import numpy as np
import pandas as pd

# non-interactive backend, set before importing pyplot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import allel


def parse_bool(x: str) -> bool:
    return str(x).strip().lower() in {"1", "true", "t", "yes", "y"}

def get_args():
    ap = argparse.ArgumentParser(description="ROH debug plot + tables")
    ap.add_argument("--vcf", required=True, help="VCF/BCF path (ideally per-chromosome)")
    ap.add_argument("--chrom", required=True, help="Chromosome label (e.g., '1' or 'chr1')")
    ap.add_argument("--roh", required=True, help="CSV with columns: sample, chrom, from, to [, cM, ...]")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--include-missing", default="true", type=parse_bool)
    ap.add_argument("--include-het", default="true", type=parse_bool)
    ap.add_argument("--per-sample", default="true", type=parse_bool,
                    help="One figure per sample (else one big figure)")
    ap.add_argument("--max-rows-per-fig", default=30, type=int)
    ap.add_argument("--per-sample-files", default="true", type=parse_bool,
                    help="Also dump per-sample CSVs for segments and summaries")
    ap.add_argument("--region", default=None, help="Optional region string for scikit-allel (e.g., '1', 'chr1', '1:1-100000')")
    return ap.parse_args()


GT_COLORS = {
    -1: 'red',
    0:  'grey',
    1:  'blue',
    2:  'black'
}

def _normalize_chrom(x):
    return str(x).replace("chr", "")

def _compress_runs_idx(pos, gt):
    """Return (start_idx, end_idx) arrays for constant-GT runs over sorted positions."""
    if len(gt) == 0:
        return np.array([], dtype=int), np.array([], dtype=int)
    change = np.r_[0, np.flatnonzero(gt[1:] != gt[:-1]) + 1]
    end = np.r_[change[1:] - 1, gt.size - 1]
    return change, end


def build_combined_df(callset):
    """Vectorized build of long-form (sample, pos, GT) using to_n_alt."""
    gt = allel.GenotypeArray(callset['calldata/GT'])
    nalt = gt.to_n_alt(fill=-1)
    pos = callset['variants/POS']
    samples = callset['samples']

    df = (
        pd.DataFrame(nalt, columns=samples)
          .assign(pos=pos)
          .melt(id_vars='pos', var_name='sample', value_name='GT')
          .sort_values(['sample', 'pos'], kind='mergesort')
          .reset_index(drop=True)
    )
    df['sample'] = df['sample'].astype('category')
    df['GT'] = df['GT'].astype(np.int8)
    return df

def compute_roh_segments_df(combined_df, roh_rows, chrom_label, include_het=True, include_missing=True):
    """
    Build tidy dataframe of segments for the provided ROH rows.
    One row per (sample, ROH, contiguous GT run), clipped to ROH bounds.
    """
    allowed = {0, 2}
    if include_het: allowed.add(1)
    if include_missing: allowed.add(-1)

    records = []
    for _, r in roh_rows.iterrows():
        sname = r['sample']
        r_from, r_to = int(r['from']), int(r['to'])

        gsub = (combined_df[(combined_df['sample'] == sname) &
                            (combined_df['pos'] >= r_from) &
                            (combined_df['pos'] <= r_to)]
                [['pos', 'GT']]
                .sort_values('pos'))
        if gsub.empty:
            records.append({
                'sample': sname, 'chrom': chrom_label, 'roh_from': r_from, 'roh_to': r_to,
                'seg_start_bp': np.nan, 'seg_end_bp': np.nan, 'seg_gt': np.nan,
                'seg_len_bp': 0, 'seg_n_sites': 0, 'seg_frac_bp_of_roh': 0.0,
                'cM': r['cM'] if 'cM' in r and pd.notnull(r['cM']) else np.nan
            })
            continue

        gsub = gsub[gsub['GT'].isin(allowed)]
        if gsub.empty:
            continue

        pos = gsub['pos'].to_numpy()
        gt = gsub['GT'].to_numpy()
        starts, ends = _compress_runs_idx(pos, gt)

        for si, ei in zip(starts, ends):
            s_bp, e_bp = int(pos[si]), int(pos[ei])
            s_bp = max(s_bp, r_from)
            e_bp = min(e_bp, r_to)
            if e_bp < s_bp:
                continue

            seg_gt = int(gt[si])
            seg_len_bp = max(1, e_bp - s_bp)
            seg_n_sites = int(ei - si + 1)
            roh_len_bp = max(1, r_to - r_from)
            frac = seg_len_bp / roh_len_bp

            records.append({
                'sample': sname,
                'chrom': chrom_label,
                'roh_from': r_from,
                'roh_to': r_to,
                'seg_start_bp': s_bp,
                'seg_end_bp': e_bp,
                'seg_gt': seg_gt,
                'seg_len_bp': seg_len_bp,
                'seg_n_sites': seg_n_sites,
                'seg_frac_bp_of_roh': frac,
                'cM': r['cM'] if 'cM' in r and pd.notnull(r['cM']) else np.nan
            })

    seg_df = pd.DataFrame.from_records(records)
    if not seg_df.empty:
        seg_df = seg_df.sort_values(['sample', 'roh_from', 'seg_start_bp'], kind='mergesort')
    return seg_df

def summarize_roh_by_sample(seg_df):
    """
    Collapse segments into per-ROH summaries and (optionally) a compact JSON list of segments per ROH row.
    Returns one row per (sample, ROH) with GT composition and majority GT.
    """
    if seg_df.empty:
        return pd.DataFrame(columns=[
            'sample','chrom','roh_from','roh_to','n_segments','roh_len_bp_est',
            'bp_gt_-1','bp_gt_0','bp_gt_1','bp_gt_2','majority_gt','majority_frac_bp','segments_json'
        ])

    grp = seg_df.groupby(['sample','chrom','roh_from','roh_to'], sort=False)

    def _agg(group):
        total_bp = group['seg_len_bp'].sum()
        comp = group.groupby('seg_gt')['seg_len_bp'].sum()
        bp_m = int(comp.get(-1, 0))
        bp_0 = int(comp.get(0, 0))
        bp_1 = int(comp.get(1, 0))
        bp_2 = int(comp.get(2, 0))
        majority_gt = int(np.array([-1,0,1,2])[np.argmax([bp_m, bp_0, bp_1, bp_2])])
        maj_bp = max(bp_m, bp_0, bp_1, bp_2)
        maj_frac = (maj_bp / total_bp) if total_bp > 0 else np.nan

        return pd.Series({
            'n_segments': int(len(group)),
            'roh_len_bp_est': int(total_bp),
            'bp_gt_-1': bp_m,
            'bp_gt_0': bp_0,
            'bp_gt_1': bp_1,
            'bp_gt_2': bp_2,
            'majority_gt': majority_gt,
            'majority_frac_bp': maj_frac,
            'segments_json': json.dumps(group[['seg_start_bp','seg_end_bp','seg_gt','seg_len_bp','seg_n_sites']].to_dict('records'))
        })

    out = grp.apply(_agg).reset_index()
    out = out.sort_values(['sample','roh_from'], kind='mergesort')
    return out

def save_roh_tables(seg_df, roh_summary_df, outdir, chrom_label, per_sample_files=False):
    os.makedirs(outdir, exist_ok=True)
    seg_path = os.path.join(outdir, f"{chrom_label}_ROH_segments.csv")
    sum_path = os.path.join(outdir, f"{chrom_label}_ROH_summary.csv")
    seg_df.to_csv(seg_path, index=False)
    roh_summary_df.to_csv(sum_path, index=False)

    if per_sample_files and not seg_df.empty:
        for s, df_s in seg_df.groupby('sample'):
            df_s.to_csv(os.path.join(outdir, f"{chrom_label}_{s}_ROH_segments.csv"), index=False)
        for s, df_s in roh_summary_df.groupby('sample'):
            df_s.to_csv(os.path.join(outdir, f"{chrom_label}_{s}_ROH_summary.csv"), index=False)

    print(f"[write] {seg_path}")
    print(f"[write] {sum_path}")


def compress_runs(positions, gts):
    """Return list of (start_pos, end_pos, GT) for constant GT runs."""
    if len(positions) == 0:
        return []
    pos = np.asarray(positions)
    gt = np.asarray(gts)
    change = np.r_[0, np.flatnonzero(gt[1:] != gt[:-1]) + 1]
    end = np.r_[change[1:] - 1, gt.size - 1]
    return list(zip(pos[change], pos[end], gt[change]))

def plot_roh_debug(
    combined_df,
    roh_df,
    chrom,
    outdir,
    include_het=True,
    include_missing=True,
    per_sample=True,
    max_rows_per_fig=30
):
    """
    Create ROH debug plots with genotype segments inside each ROH.
    - per_sample=True: one PNG per sample (recommended if many ROHs).
    - per_sample=False: one big plot for all samples on this chromosome.
    """
    def normalize_chrom(x):
        return str(x).replace("chr", "")
    chrom_norm = normalize_chrom(chrom)

    roh_chr = roh_df.copy()
    roh_chr['chrom_norm'] = roh_chr['chrom'].apply(normalize_chrom)
    roh_chr = roh_chr[roh_chr['chrom_norm'] == chrom_norm]

    if roh_chr.empty:
        print(f"[warn] No ROH rows for chrom {chrom}.")
        return

    samples_in_data = set(combined_df['sample'].unique())
    roh_chr = roh_chr[roh_chr['sample'].isin(samples_in_data)]
    if roh_chr.empty:
        print(f"[warn] No overlapping samples between ROH table and genotype data on chrom {chrom}.")
        return

    allowed_gts = {0, 2}
    if include_het:
        allowed_gts.add(1)
    if include_missing:
        allowed_gts.add(-1)

    def _plot_rows(rows, title, outfile):
        n = len(rows)
        if n == 0:
            return
        fig_h = max(3, min(0.6 * n + 1.5, 18))
        fig, ax = plt.subplots(figsize=(16, fig_h))

        xmin = rows['from'].min()
        xmax = rows['to'].max()
        span = xmax - xmin
        margin = max(1_000, int(0.01 * span)) if span > 0 else 1_000
        ax.set_xlim(xmin - margin, xmax + margin)

        yticklabels = []
        y_positions = np.arange(n)
        band_height = 0.8

        for idx, (_, r) in enumerate(rows.reset_index(drop=True).iterrows()):
            y = y_positions[idx]

            ax.add_patch(Rectangle(
                (r['from'], y - band_height/2),
                r['to'] - r['from'],
                band_height,
                facecolor='#f9e79f',
                edgecolor='#b7950b',
                linewidth=0.8,
                alpha=0.6
            ))

            gsub = (combined_df[(combined_df['sample'] == r['sample']) &
                                (combined_df['pos'] >= r['from']) &
                                (combined_df['pos'] <= r['to'])]
                    [['pos', 'GT']]
                    .sort_values('pos'))
            gsub = gsub[gsub['GT'].isin(allowed_gts)]
            segs = compress_runs(gsub['pos'].to_numpy(), gsub['GT'].to_numpy())

            for s, e, gt in segs:
                color = GT_COLORS.get(int(gt), 'white')
                s2 = max(int(s), int(r['from']))
                e2 = min(int(e), int(r['to']))
                if e2 >= s2:
                    ax.add_patch(Rectangle(
                        (s2, y - (band_height*0.45)),
                        max(1, e2 - s2),
                        band_height*0.9,
                        facecolor=color,
                        edgecolor='none',
                        alpha=1.0
                    ))

            label = f"{r['sample']} | {chrom}:{int(r['from']):,}-{int(r['to']):,}"
            if 'cM' in r and pd.notnull(r['cM']):
                label += f" | {float(r['cM']):.2f} cM"
            yticklabels.append(label)

        ax.set_yticks(y_positions)
        ax.set_yticklabels(yticklabels, fontsize=9)
        ax.set_xlabel("Genomic position (bp)")
        ax.set_title(title, loc='left')

        legend_patches = [
            Rectangle((0,0),1,1, facecolor='#f9e79f', edgecolor='#b7950b', alpha=0.6, label='ROH interval'),
        ]
        for gt_code, name in [(0, 'Hom-Ref (0)'), (2, 'Hom-Alt (2)'), (1, 'Het (1)'), (-1, 'Missing (-1)')]:
            if gt_code in allowed_gts:
                legend_patches.append(Rectangle((0,0),1,1, facecolor=GT_COLORS[gt_code], label=name))
        ax.legend(handles=legend_patches, loc='upper right', frameon=True)

        for spine in ['top','right']:
            ax.spines[spine].set_visible(False)

        plt.tight_layout()
        outpath = os.path.join(outdir, outfile)
        plt.savefig(outpath, dpi=300, bbox_inches='tight')
        plt.close(fig)
        print(f"[write] {outpath}")

    if per_sample:
        for sample, rows in roh_chr.groupby('sample'):
            rows = rows.sort_values(['from', 'to'])
            total = len(rows)
            pages = math.ceil(total / max_rows_per_fig)
            for p in range(pages):
                chunk = rows.iloc[p*max_rows_per_fig:(p+1)*max_rows_per_fig]
                title = f"ROH Genotype Segments | chr{chrom_norm} | {sample}"
                suffix = f"_p{p+1}" if pages > 1 else ""
                outfile = f"{str(chrom)}_{sample}_ROH_debug{suffix}.png"
                _plot_rows(chunk, title, outfile)
    else:
        rows = roh_chr.sort_values(['sample', 'from', 'to'])
        title = f"ROH Genotype Segments | chr{chrom_norm} | All samples"
        outfile = f"{str(chrom)}_ALL_ROH_debug.png"
        _plot_rows(rows, title, outfile)


def main():
    args = get_args()
    os.makedirs(args.outdir, exist_ok=True)

    roh_df = pd.read_csv(args.roh)
    required_cols = {'sample','chrom','from','to'}
    missing = required_cols - set(map(str, roh_df.columns))
    if missing:
        raise ValueError(f"ROH CSV missing required columns: {missing}")

    for col in ['from','to']:
        roh_df[col] = pd.to_numeric(roh_df[col], errors='coerce').astype('Int64')

    fields = ['POS', 'samples', 'calldata/GT']
    callset = allel.read_vcf(args.vcf, fields=fields, region=args.region)

    combined_df = build_combined_df(callset)

    chrom_norm = _normalize_chrom(args.chrom)
    roh_chr = roh_df.copy()
    roh_chr['chrom_norm'] = roh_chr['chrom'].apply(_normalize_chrom)
    roh_chr = roh_chr[roh_chr['chrom_norm'] == chrom_norm]
    roh_chr = roh_chr[roh_chr['sample'].isin(set(combined_df['sample'].unique()))]
    if roh_chr.empty:
        print(f"[warn] No ROHs for chrom {args.chrom} with samples present in VCF. Exiting after building combined_df.")
        return

    seg_df = compute_roh_segments_df(
        combined_df=combined_df,
        roh_rows=roh_chr,
        chrom_label=f"chr{chrom_norm}",
        include_het=args.include_het,
        include_missing=args.include_missing
    )

    roh_summary_df = summarize_roh_by_sample(seg_df)

    save_roh_tables(seg_df, roh_summary_df, args.outdir, chrom_label=f"chr{chrom_norm}", per_sample_files=args.per_sample_files)

    plot_roh_debug(
        combined_df=combined_df,
        roh_df=roh_df,
        chrom=args.chrom,
        outdir=args.outdir,
        include_het=args.include_het,
        include_missing=args.include_missing,
        per_sample=args.per_sample,
        max_rows_per_fig=args.max_rows_per_fig
    )

if __name__ == "__main__":
    main()
