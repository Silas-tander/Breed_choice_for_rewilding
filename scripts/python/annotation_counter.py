"""Unique variant positions per SnpEff annotation of one chromosome VCF.
Usage: python annotation_counter.py <chr> <outdir> <vcf>  ->  <outdir>/<chr>.counts.csv
"""

import allel
import os
import sys
import pandas as pd
from collections import Counter

chr = sys.argv[1]
outdir = sys.argv[2]
vcf_path = sys.argv[3]

os.makedirs(outdir, exist_ok=True)

callset = allel.read_vcf(vcf_path, fields=['variants/ANN', 'variants/POS'])

annotations = callset['variants/ANN']
pos = callset['variants/POS']

def parse_annotation(annotation):
    fields = annotation.split('|')
    effect = fields[1]
    return effect

annotation_counts = Counter()

for variant_idx, variant_ann in enumerate(annotations):
    if variant_ann:
        effect = parse_annotation(variant_ann)
        site_pos = pos[variant_idx]
        annotation_counts[effect] += 1

counts_df = pd.DataFrame(annotation_counts.items(), columns=['Annotation', 'Count'])

output_file = os.path.join(outdir, f'{chr}.counts.csv')
counts_df.to_csv(output_file, index=False)
