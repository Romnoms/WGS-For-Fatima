#!/usr/bin/env python3
"""
Genomic Landscape Figure — Cancer Discovery Style
Inspired by Li et al. Cancer Discovery 2025 (CD-24-1379) Figure 1.
Multi-panel figure showing:
  Panel A: TMB bar chart with variant type breakdown (SNV/Indel/MNV)
  Panel B: Relative mutation spectrum (6-class substitution) per sample
  Panel C: Rainfall plot — inter-mutation distance across the genome
  Panel D: Oncoplot — gene x sample matrix of recurrent HIGH/MODERATE mutations

Data: Filtered somatic variants from sucrose-induced tumors in C57BL/6J mice.
"""

import os
import gzip
import glob
import re
from collections import Counter, defaultdict
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Patch
from matplotlib.gridspec import GridSpec
from matplotlib.colors import ListedColormap, BoundaryNorm
import warnings
warnings.filterwarnings('ignore')

plt.rcParams.update({
    'font.size': 10,
    'font.family': 'sans-serif',
    'axes.titlesize': 12,
    'axes.labelsize': 10,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'figure.facecolor': 'white',
})

DATA_DIR = os.path.expanduser('~/wgs-project/data/filtered_somatic')
OUT_DIR = os.path.expanduser('~/wgs-project/figures_filtered')
os.makedirs(OUT_DIR, exist_ok=True)

SAMPLES = sorted([
    os.path.basename(f).replace('.filtered.vcf.gz', '')
    for f in glob.glob(f'{DATA_DIR}/*.filtered.vcf.gz')
])
SAMPLE_SHORT = {s: s.replace('26034XD-04-', 'S') for s in SAMPLES}
LABELS = [SAMPLE_SHORT[s] for s in SAMPLES]
GENOME_SIZE_MB = 2728.22
COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

SUB_COLORS = {
    'C>A': '#1EBBD7', 'C>G': '#050708', 'C>T': '#E62725',
    'T>A': '#CBCACB', 'T>C': '#A1CF64', 'T>G': '#EDC8C5',
}
TYPE_COLORS = {'SNV': '#3498db', 'Insertion': '#2ecc71', 'Deletion': '#e74c3c', 'MNV': '#f39c12'}
IMPACT_COLORS = {'HIGH': '#d62728', 'MODERATE': '#ff7f0e', 'LOW': '#2ca02c', 'MODIFIER': '#aec7e8'}

# GRCm39 chromosome sizes (bp)
CHROM_SIZES = {
    '1': 195154279, '2': 181755017, '3': 159745316, '4': 156860686,
    '5': 151758149, '6': 149588044, '7': 144995196, '8': 130127694,
    '9': 124359700, '10': 130530862, '11': 121973369, '12': 120092757,
    '13': 120883175, '14': 125139656, '15': 104073951, '16': 98008968,
    '17': 95294699, '18': 90720763, '19': 61420004, 'X': 169476592,
    'Y': 91455967,
}
MAIN_CHROMS = [str(i) for i in range(1, 20)] + ['X', 'Y']

# Cumulative offsets for genome-wide coordinate
cum_offset = {}
running = 0
for c in MAIN_CHROMS:
    cum_offset[c] = running
    running += CHROM_SIZES[c]
GENOME_TOTAL = running


def parse_filtered_vcf(vcf_path):
    variants = []
    csq_fields = None
    with gzip.open(vcf_path, 'rt') as fh:
        for line in fh:
            if line.startswith('##INFO=<ID=CSQ'):
                fmt = re.search(r'Format: (.+?)"', line)
                if fmt:
                    csq_fields = fmt.group(1).split('|')
                continue
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            chrom, pos, ref, alt = parts[0], int(parts[1]), parts[3], parts[4]
            info = parts[7]
            fmt_keys = parts[8].split(':')
            fmt_vals = parts[9].split(':')
            fmt_dict = dict(zip(fmt_keys, fmt_vals))
            af = float(fmt_dict.get('AF', '0'))
            dp = int(fmt_dict.get('DP', '0'))

            if len(ref) == 1 and len(alt) == 1:
                var_type = 'SNV'
            elif len(ref) < len(alt):
                var_type = 'Insertion'
            elif len(ref) > len(alt):
                var_type = 'Deletion'
            else:
                var_type = 'MNV'

            csq_str = ''
            for item in info.split(';'):
                if item.startswith('CSQ='):
                    csq_str = item[4:]
                    break
            gene, consequence, impact = '', '', ''
            if csq_str and csq_fields:
                first_csq = csq_str.split(',')[0]
                csq_vals = first_csq.split('|')
                csq_dict = dict(zip(csq_fields, csq_vals))
                gene = csq_dict.get('SYMBOL', '')
                consequence = csq_dict.get('Consequence', '')
                impact = csq_dict.get('IMPACT', '')

            variants.append({
                'chrom': chrom, 'pos': pos, 'ref': ref, 'alt': alt,
                'af': af, 'dp': dp, 'var_type': var_type,
                'gene': gene, 'consequence': consequence, 'impact': impact
            })
    return pd.DataFrame(variants)


def get_substitution_class(ref, alt):
    if ref in ('C', 'T'):
        return f'{ref}>{alt}'
    else:
        return f'{COMPLEMENT[ref]}>{COMPLEMENT[alt]}'


# =============================================================================
# Parse all filtered VCFs
# =============================================================================
print("Parsing filtered somatic VCFs...")
all_variants = {}
for sample in SAMPLES:
    vcf_path = f'{DATA_DIR}/{sample}.filtered.vcf.gz'
    df = parse_filtered_vcf(vcf_path)
    df['sample'] = sample
    df['short'] = SAMPLE_SHORT[sample]
    all_variants[sample] = df
    print(f"  {SAMPLE_SHORT[sample]}: {len(df)} filtered variants")

combined = pd.concat(all_variants.values(), ignore_index=True)


# =============================================================================
# Build the multi-panel landscape figure
# =============================================================================
print("\nGenerating Cancer Discovery-style genomic landscape figure...")

fig = plt.figure(figsize=(20, 22))
gs = GridSpec(4, 1, figure=fig, height_ratios=[1, 0.8, 1.2, 1.2],
              hspace=0.30)

# ---- Panel A: TMB bar chart with variant type stacked ----
ax_a = fig.add_subplot(gs[0])

# Sort samples by TMB for visual clarity
sample_tmb = {s: len(all_variants[s]) / GENOME_SIZE_MB for s in SAMPLES}
sorted_samples = sorted(SAMPLES, key=lambda s: sample_tmb[s], reverse=True)
sorted_labels = [SAMPLE_SHORT[s] for s in sorted_samples]

# Build stacked data
type_order = ['SNV', 'Insertion', 'Deletion', 'MNV']
bottoms = np.zeros(len(sorted_samples))
x_pos = np.arange(len(sorted_samples))

for vtype in type_order:
    vals = []
    for s in sorted_samples:
        df = all_variants[s]
        count = (df['var_type'] == vtype).sum() / GENOME_SIZE_MB
        vals.append(count)
    vals = np.array(vals)
    ax_a.bar(x_pos, vals, bottom=bottoms, color=TYPE_COLORS[vtype],
             label=vtype, edgecolor='white', linewidth=0.3, width=0.8)
    bottoms += vals

ax_a.set_xticks(x_pos)
ax_a.set_xticklabels(sorted_labels, rotation=0, ha='center', fontsize=9)
ax_a.set_ylabel('Mutations / Mb', fontsize=11)
ax_a.set_title('A.  Tumor Mutational Burden by Variant Type', fontsize=13,
               fontweight='bold', loc='left', pad=10)
ax_a.legend(title='Variant Type', bbox_to_anchor=(1.01, 1), loc='upper left',
            fontsize=9, title_fontsize=10)
ax_a.spines['top'].set_visible(False)
ax_a.spines['right'].set_visible(False)

# Annotate hypermutated
for i, s in enumerate(sorted_samples):
    tmb = sample_tmb[s]
    if tmb > 4:
        ax_a.annotate(f'{tmb:.1f}', (i, tmb + 0.1), ha='center', fontsize=8,
                      fontweight='bold', color='#e74c3c')

# ---- Panel B: Relative mutation spectrum (100% stacked) ----
ax_b = fig.add_subplot(gs[1])

sub_order = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
bottoms = np.zeros(len(sorted_samples))

for sub in sub_order:
    vals = []
    for s in sorted_samples:
        df = all_variants[s]
        snvs = df[df['var_type'] == 'SNV']
        total_snv = len(snvs)
        if total_snv == 0:
            vals.append(0)
            continue
        count = sum(1 for _, r in snvs.iterrows()
                    if get_substitution_class(r['ref'], r['alt']) == sub)
        vals.append(count / total_snv * 100)
    vals = np.array(vals)
    ax_b.bar(x_pos, vals, bottom=bottoms, color=SUB_COLORS[sub],
             label=sub, edgecolor='white', linewidth=0.3, width=0.8)
    bottoms += vals

ax_b.set_xticks(x_pos)
ax_b.set_xticklabels(sorted_labels, rotation=0, ha='center', fontsize=9)
ax_b.set_ylabel('Proportion (%)', fontsize=11)
ax_b.set_ylim(0, 100)
ax_b.set_title('B.  Relative SNV Mutation Spectrum (6-class)', fontsize=13,
               fontweight='bold', loc='left', pad=10)
ax_b.legend(title='Substitution', bbox_to_anchor=(1.01, 1), loc='upper left',
            fontsize=9, title_fontsize=10, ncol=1)
ax_b.spines['top'].set_visible(False)
ax_b.spines['right'].set_visible(False)

# ---- Panel C: Rainfall plot for select samples ----
ax_c = fig.add_subplot(gs[2])

# Pick 4 representative samples
rainfall_samples = []
for s in sorted_samples:
    short = SAMPLE_SHORT[s]
    if short in ('S15', 'S04', 'S16', 'S01'):
        rainfall_samples.append(s)

# Use top 4 by TMB if those aren't all found
if len(rainfall_samples) < 4:
    rainfall_samples = sorted_samples[:4]

# Colors for rainfall by substitution type
rain_colors = {
    'C>A': '#1EBBD7', 'C>G': '#050708', 'C>T': '#E62725',
    'T>A': '#CBCACB', 'T>C': '#A1CF64', 'T>G': '#EDC8C5',
}

# Plot all 4 samples overlaid
for si, sample in enumerate(rainfall_samples):
    df = all_variants[sample]
    snvs = df[(df['var_type'] == 'SNV') & (df['chrom'].isin(MAIN_CHROMS))].copy()
    snvs = snvs.sort_values(['chrom', 'pos'],
                            key=lambda x: x.map(lambda v: MAIN_CHROMS.index(v) if v in MAIN_CHROMS else 999) if x.name == 'chrom' else x)

    # Compute genome-wide position
    snvs['genome_pos'] = snvs.apply(
        lambda r: cum_offset.get(r['chrom'], 0) + r['pos'], axis=1)
    snvs = snvs.sort_values('genome_pos')

    # Inter-mutation distance
    snvs['imd'] = snvs['genome_pos'].diff().fillna(1e6)
    snvs['imd'] = snvs['imd'].clip(lower=1)
    snvs['log_imd'] = np.log10(snvs['imd'])

    # Substitution class for color
    snvs['sub'] = snvs.apply(lambda r: get_substitution_class(r['ref'], r['alt']), axis=1)

    if si == 0:  # Only plot the first (highest TMB) fully, others as background
        for sub in sub_order:
            mask = snvs['sub'] == sub
            ax_c.scatter(snvs.loc[mask, 'genome_pos'] / 1e6,
                         snvs.loc[mask, 'log_imd'],
                         c=rain_colors[sub], s=1.5, alpha=0.6, label=sub,
                         rasterized=True)
    else:
        ax_c.scatter(snvs['genome_pos'] / 1e6, snvs['log_imd'],
                     c='#d0d0d0', s=0.3, alpha=0.15, rasterized=True)

# Add chromosome boundaries
for c in MAIN_CHROMS:
    bnd = cum_offset[c] / 1e6
    ax_c.axvline(x=bnd, color='gray', linewidth=0.3, alpha=0.4)

# Chromosome labels
for c in MAIN_CHROMS:
    mid = (cum_offset[c] + CHROM_SIZES[c] / 2) / 1e6
    label = c if c not in ('18', '19', 'Y') else ''
    if c in ('18', '19', 'Y'):
        label = c
    ax_c.text(mid, -0.3, c, ha='center', va='top', fontsize=6, color='gray')

ax_c.set_ylabel('log10(Inter-mutation distance)', fontsize=11)
ax_c.set_xlabel('Genome position (Mb)', fontsize=11)
ax_c.set_title(f'C.  Rainfall Plot — {SAMPLE_SHORT[rainfall_samples[0]]} '
               f'(colored) vs cohort (gray)', fontsize=13,
               fontweight='bold', loc='left', pad=10)
ax_c.set_ylim(0, 8)
ax_c.set_xlim(0, GENOME_TOTAL / 1e6)

# Add kataegis threshold line
ax_c.axhline(y=3, color='#e74c3c', linestyle='--', alpha=0.5, linewidth=1)
ax_c.text(GENOME_TOTAL / 1e6 * 0.98, 3.15, 'Kataegis threshold (1 kb)',
          ha='right', fontsize=8, color='#e74c3c', alpha=0.7)

ax_c.legend(title='Substitution', bbox_to_anchor=(1.01, 1), loc='upper left',
            fontsize=8, title_fontsize=9, markerscale=5)
ax_c.spines['top'].set_visible(False)
ax_c.spines['right'].set_visible(False)

# ---- Panel D: Oncoplot — gene x sample heatmap ----
ax_d = fig.add_subplot(gs[3])

coding = combined[(combined['impact'].isin(['HIGH', 'MODERATE'])) & (combined['gene'] != '')].copy()

# Filter gene families
skip_prefixes = ('Vmn', 'Zfp', 'Or', 'Gm', 'Olfr')
coding = coding[~coding['gene'].str.startswith(tuple(skip_prefixes))]

gene_sample_sets = coding.groupby('gene')['short'].apply(set)
genes_recurrent = {g: s for g, s in gene_sample_sets.items() if len(s) >= 3}
genes_sorted = sorted(genes_recurrent.items(), key=lambda x: len(x[1]), reverse=True)[:25]
gene_names = [g for g, _ in genes_sorted]

# Build matrix
impact_val = {'HIGH': 3, 'MODERATE': 2}
# Use sorted_labels order for columns
matrix = np.zeros((len(gene_names), len(sorted_labels)))
for i, gene in enumerate(gene_names):
    gene_data = coding[coding['gene'] == gene]
    for _, row in gene_data.iterrows():
        if row['short'] in sorted_labels:
            j = sorted_labels.index(row['short'])
            val = impact_val.get(row['impact'], 0)
            matrix[i, j] = max(matrix[i, j], val)

cmap = ListedColormap(['#f5f5f5', '#ff7f0e', '#d62728'])
bounds = [0, 1.5, 2.5, 3.5]
norm = BoundaryNorm(bounds, cmap.N)

im = ax_d.imshow(matrix, aspect='auto', cmap=cmap, norm=norm, interpolation='none')

ax_d.set_xticks(range(len(sorted_labels)))
ax_d.set_xticklabels(sorted_labels, rotation=0, ha='center', fontsize=9)
ax_d.set_yticks(range(len(gene_names)))
ax_d.set_yticklabels(gene_names, fontsize=9)
ax_d.set_title('D.  Somatic Gene Mutation Landscape\n'
               '     HIGH/MODERATE impact, recurrent in >=3 samples',
               fontsize=12, fontweight='bold', loc='left', pad=8)

# Recurrence annotation on right
ax_d2 = ax_d.twinx()
ax_d2.set_ylim(ax_d.get_ylim())
ax_d2.set_yticks(range(len(gene_names)))
recurrence = [len(genes_recurrent[g]) for g in gene_names]
ax_d2.set_yticklabels([f'{r}/16' for r in recurrence], fontsize=8)
ax_d2.set_ylabel('Samples', fontsize=10)

# Grid lines
for i in range(len(gene_names) + 1):
    ax_d.axhline(y=i - 0.5, color='white', linewidth=1.5)
for j in range(len(sorted_labels) + 1):
    ax_d.axvline(x=j - 0.5, color='white', linewidth=1.5)

# Legend
legend_elements = [
    Patch(facecolor='#f5f5f5', edgecolor='gray', label='Wild-type'),
    Patch(facecolor='#ff7f0e', label='MODERATE'),
    Patch(facecolor='#d62728', label='HIGH'),
]
ax_d.legend(handles=legend_elements, loc='lower right', fontsize=9, title='Impact')

# Main title
fig.suptitle('Genomic Landscape of Somatic Mutations in Sucrose-Induced Mouse Tumors',
             fontsize=16, fontweight='bold', y=0.995)

plt.savefig(f'{OUT_DIR}/17_genomic_landscape.png')
plt.close()
print("  -> 17_genomic_landscape.png")


# =============================================================================
# Bonus: Individual sample rainfall plots (2x2 grid)
# =============================================================================
print("Generating individual rainfall plots...")

fig, axes = plt.subplots(2, 2, figsize=(18, 10))
fig.suptitle('Rainfall Plots — Inter-Mutation Distance Across the Genome',
             fontsize=15, fontweight='bold')

plot_samples_rain = []
for short_name in ['S15', 'S04', 'S16', 'S01']:
    for s in SAMPLES:
        if SAMPLE_SHORT[s] == short_name:
            plot_samples_rain.append(s)
            break

# Fallback
while len(plot_samples_rain) < 4:
    for s in sorted_samples:
        if s not in plot_samples_rain:
            plot_samples_rain.append(s)
            break

for ax, sample in zip(axes.flat, plot_samples_rain):
    df = all_variants[sample]
    snvs = df[(df['var_type'] == 'SNV') & (df['chrom'].isin(MAIN_CHROMS))].copy()
    snvs['genome_pos'] = snvs.apply(
        lambda r: cum_offset.get(r['chrom'], 0) + r['pos'], axis=1)
    snvs = snvs.sort_values('genome_pos')
    snvs['imd'] = snvs['genome_pos'].diff().fillna(1e6).clip(lower=1)
    snvs['log_imd'] = np.log10(snvs['imd'])
    snvs['sub'] = snvs.apply(lambda r: get_substitution_class(r['ref'], r['alt']), axis=1)

    for sub in sub_order:
        mask = snvs['sub'] == sub
        ax.scatter(snvs.loc[mask, 'genome_pos'] / 1e6,
                   snvs.loc[mask, 'log_imd'],
                   c=rain_colors[sub], s=2, alpha=0.5, label=sub,
                   rasterized=True)

    for c in MAIN_CHROMS:
        ax.axvline(x=cum_offset[c] / 1e6, color='gray', linewidth=0.3, alpha=0.3)

    ax.axhline(y=3, color='#e74c3c', linestyle='--', alpha=0.4, linewidth=0.8)
    ax.set_ylim(0, 8)
    ax.set_xlim(0, GENOME_TOTAL / 1e6)
    short = SAMPLE_SHORT[sample]
    tmb = len(df) / GENOME_SIZE_MB
    ax.set_title(f'{short} (TMB: {tmb:.1f} mut/Mb, n={len(snvs):,} SNVs)',
                 fontweight='bold', fontsize=11)
    ax.set_ylabel('log10(IMD)')
    ax.set_xlabel('Genome position (Mb)')

# Single legend
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncol=6, fontsize=9,
           markerscale=4, title='Substitution Type',
           bbox_to_anchor=(0.5, -0.02))

plt.tight_layout(rect=[0, 0.03, 1, 0.96])
plt.savefig(f'{OUT_DIR}/18_rainfall_plots.png')
plt.close()
print("  -> 18_rainfall_plots.png")

print("\nLandscape figures complete!")
