"""SnpEff-annotated variant counts per ROH interval and sample for one chromosome.
Usage: python ROH_genload.py <chr> <outdir> <vcf> <roh.csv>
"""

import allel
import os
import sys
import pandas as pd
import numpy as np

chr = sys.argv[1]
outdir = sys.argv[2]
vcf_path = sys.argv[3]
ROH_path = sys.argv[4]

os.makedirs(outdir, exist_ok=True)

callset = allel.read_vcf(vcf_path, fields=['calldata/GT', 'variants/ANN', 'variants/POS', 'samples'])

samples = callset['samples']
genotypes = callset['calldata/GT']
annotations = callset['variants/ANN']
positions = callset['variants/POS']

roh_df = pd.read_csv(ROH_path)

chr_int = int(chr.lstrip('chr'))

roh_df = roh_df[roh_df['chrom'] == chr_int]

def parse_annotation(annotation):
    fields = annotation.split('|')
    effect = fields[1]
    impact = fields[2]
    return effect, impact

results = []

for sample_idx, sample in enumerate(samples):
    sample_roh = roh_df[roh_df['sample'] == sample]

    for _, interval in sample_roh.iterrows():
        interval_start = interval['from']
        interval_end = interval['to']
        length_bps = interval['lengthBps']
        interval_id = f"{sample}_{interval_start}-{interval_end}"

        interval_variants = (positions >= interval_start) & (positions <= interval_end)

        effect_counts = {}
        impact_counts = {}

        for variant_idx, within_interval in enumerate(interval_variants):
            if within_interval:
                genotype = genotypes[variant_idx, sample_idx]
                annotation = annotations[variant_idx]

                if annotation:
                    effect, impact = parse_annotation(annotation)

                    if np.array_equal(genotype, [1, 1]):
                        count = 2
                    elif np.array_equal(genotype, [1, 0]) or np.array_equal(genotype, [0, 1]):
                        count = 1
                    else:
                        count = 0

                    if count > 0:
                        if effect not in effect_counts:
                            effect_counts[effect] = 0
                        effect_counts[effect] += count

                        if impact not in impact_counts:
                            impact_counts[impact] = 0
                        impact_counts[impact] += count

        results.append({
            'sample': sample,
            'chrom': chr,
            'interval': f"{interval_start}-{interval_end}",
            'lengthBps': length_bps,
            'effect_counts': effect_counts,
            'impact_counts': impact_counts,
        })

output_df = pd.DataFrame(results)

output_file = os.path.join(outdir, f'{chr}_roh_variant_counts_10miss.csv')
output_df.to_csv(output_file, index=False)
