"""Debug check: first variant whose first ANN field equals REF.
Usage: python test_annotations.py <chr> <outdir> <vcf>
"""

import allel
import os
import sys
import pandas as pd
import numpy as np

chr = sys.argv[1]
outdir = sys.argv[2]
vcf_path = sys.argv[3]

os.makedirs(outdir, exist_ok=True)

callset = allel.read_vcf(vcf_path, fields=['variants/ANN', 'variants/REF', 'POS'])

pos = callset['variants/POS']
annotations = callset['variants/ANN']
ref = callset['variants/REF']

ref_sites = []

for i in range(len(ref)):
    if annotations[i][0] == ref[i]:
        ref_sites.append(pos[i])
        break

output_df = pd.DataFrame(ref_sites)
output_file = os.path.join(outdir, f'{chr}_ref_sites.csv')
output_df.to_csv(output_file, index=False)
