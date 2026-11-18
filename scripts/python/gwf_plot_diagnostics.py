"""QUAL and DP diagnostic plots of one gVCF.
Usage: python gwf_plot_diagnostics.py <vcf> <outdir>
"""

import os
import pandas as pd
import sys
import allel
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

path = sys.argv[1]
diri = sys.argv[2]

if not os.path.exists(diri):
    os.makedirs(diri)

callset = allel.read_vcf(path, fields=['QUAL', 'is_snp', 'DP'])

df = pd.DataFrame({'snp': callset['variants/is_snp'], 'qual': callset['variants/QUAL'], 'dp': callset['variants/DP']})

df_clean = df.dropna(subset=['qual'])


plt.figure(figsize=(10, 6))
ax = sns.histplot(data=df_clean, x='qual', hue='snp', multiple='stack', bins=8000, kde=False)
ax.set_xlabel('Quality')
ax.set_ylabel('Frequency')
ax.set_title('Histogram of Qual Values Stratified by SNP Status axlim')
ax.set_xlim(0, 50000)
ax.set_ylim(0, 5000)

plot_title = ax.get_title()
filename = plot_title.replace(' ', '_') + '.png'

plt.savefig(f"{diri}/{filename}")
plt.close()


plt.figure(figsize=(10, 6))
ax2 = sns.histplot(data=df_clean, x='qual', hue='snp', multiple='stack', bins=8000, kde=False)
ax2.set_xlabel('Quality')
ax2.set_ylabel('Frequency')
ax2.set_title('Histogram of Qual Values Stratified by SNP Status')
ax2.set_ylim(0, 7000)

plot_title2 = ax2.get_title()
filename2 = plot_title2.replace(' ', '_') + '.png'

plt.savefig(f"{diri}/{filename2}")
plt.close()

df_clean2 = df.dropna()

plt.figure(figsize=(10, 6))
ax3 = sns.scatterplot(data=df_clean2, x='qual', y='dp', hue='snp')
ax3.set_title('quality depth scatterplot')

plot_title3 = ax3.get_title()
filename3 = plot_title3.replace(' ', '_') + '.png'

plt.savefig(f"{diri}/{filename3}")
plt.close()

df_cleany = df.dropna(subset=['qual', 'dp'])
qd = df_cleany['qual'].div(df_cleany['dp'].values)

df_clean3 = pd.DataFrame({
    'qual': df_cleany['qual'],
    'qd': qd,
    'dp': df_cleany['qual'],
    'snp': df_cleany['snp']
})

plt.figure(figsize=(10, 6))
ax4 = sns.histplot(data=df_clean3, x='qd', hue='snp', multiple='stack', bins=8000, kde=False)
ax4.set_xlabel('Quality depth')
ax4.set_ylabel('Frequency')
ax4.set_title('Histogram of Qual depth Values Stratified by SNP Status')
ax4.set_ylim(0, 7000)

plot_title4 = ax4.get_title()
filename4 = plot_title4.replace(' ', '_') + '.png'

plt.savefig(f"{diri}/{filename4}")
plt.close()

plt.figure(figsize=(10, 6))
ax5 = sns.histplot(data=df, x='dp', hue='snp', multiple='stack', bins=8000, kde=False)
ax5.set_xlabel('Depth')
ax5.set_ylabel('Frequency')
ax5.set_title('Histogram of Depth Values Stratified by SNP Status')
ax5.set_xlim(0, 15000)

plot_title5 = ax5.get_title()
filename5 = plot_title5.replace(' ', '_') + '.png'

plt.savefig(f"{diri}/{filename5}")
plt.close()
