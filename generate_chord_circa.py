#!/usr/bin/env python3
"""
Chord Diagrams & Circa-Style Plots — WGS Somatic Variants
Sucrose-Induced Tumors in C57BL/6J Mice (GRCm39)

Generates:
  1. Chord diagram: chromosome-to-chromosome shared gene mutations
  2. Chord diagram: sample-to-sample co-mutated genes
  3. Chord diagram: pathway-to-pathway co-mutation
  4. Circa-style plot: genome-wide scatter with multi-track data rings
  5. Circa-style plot: per-sample comparison (2x2)
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
from matplotlib.patches import Patch, FancyArrowPatch
from matplotlib.gridspec import GridSpec
import warnings
warnings.filterwarnings('ignore')

from pycirclize import Circos

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

SUB_COLORS = {
    'C>A': '#1EBBD7', 'C>G': '#050708', 'C>T': '#E62725',
    'T>A': '#CBCACB', 'T>C': '#A1CF64', 'T>G': '#EDC8C5',
}
SUB_ORDER = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']

IMPACT_COLORS = {'HIGH': '#d62728', 'MODERATE': '#ff7f0e', 'LOW': '#2ca02c', 'MODIFIER': '#aec7e8'}

# Chromosome colors (alternating for visual clarity)
CHROM_COLORS = {}
palette = ['#4e79a7', '#f28e2b', '#e15759', '#76b7b2', '#59a14f',
           '#edc948', '#b07aa1', '#ff9da7', '#9c755f', '#bab0ac',
           '#4e79a7', '#f28e2b', '#e15759', '#76b7b2', '#59a14f',
           '#edc948', '#b07aa1', '#ff9da7', '#9c755f', '#bab0ac', '#aaa']
for i, c in enumerate(MAIN_CHROMS):
    CHROM_COLORS[c] = palette[i]

WINDOW_SIZE = 2_000_000


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


def get_sub_class(ref, alt):
    if ref in ('C', 'T'):
        return f'{ref}>{alt}'
    return f'{COMPLEMENT[ref]}>{COMPLEMENT[alt]}'


# =============================================================================
# Parse all VCFs
# =============================================================================
print("Parsing filtered somatic VCFs...")
all_variants = {}
for sample in SAMPLES:
    vcf_path = f'{DATA_DIR}/{sample}.filtered.vcf.gz'
    df = parse_filtered_vcf(vcf_path)
    df['sample'] = sample
    df['short'] = SAMPLE_SHORT[sample]
    all_variants[sample] = df
    print(f"  {SAMPLE_SHORT[sample]}: {len(df)} variants")

combined = pd.concat(all_variants.values(), ignore_index=True)
coding = combined[(combined['impact'].isin(['HIGH', 'MODERATE'])) & (combined['gene'] != '')].copy()

# Filter gene families
skip_prefixes = ('Vmn', 'Zfp', 'Or', 'Gm', 'Olfr')
coding = coding[~coding['gene'].str.startswith(tuple(skip_prefixes))]


# =============================================================================
# 1. CHORD DIAGRAM: Chromosome-to-Chromosome Shared Gene Mutations
# =============================================================================
print("\nGenerating chromosome-to-chromosome chord diagram...")

# Build adjacency: for each gene mutated on chrom_i, if same sample has a
# coding mutation on chrom_j, count the link
chrom_gene_sample = coding.groupby(['chrom', 'gene', 'short']).size().reset_index()
chrom_gene_sample.columns = ['chrom', 'gene', 'short', 'count']

# For each sample, find all chromosomes with coding mutations
sample_chroms = coding.groupby('short')['chrom'].apply(set)

# Build adjacency matrix: number of samples with coding mutations on BOTH chromosomes
chrom_list = [c for c in MAIN_CHROMS if c in coding['chrom'].unique()]
n_chrom = len(chrom_list)
adj_chrom = np.zeros((n_chrom, n_chrom))

for short_name, chroms in sample_chroms.items():
    chroms_in = [c for c in chroms if c in chrom_list]
    for i, c1 in enumerate(chrom_list):
        for j, c2 in enumerate(chrom_list):
            if i < j and c1 in chroms_in and c2 in chroms_in:
                adj_chrom[i, j] += 1
                adj_chrom[j, i] += 1

# Use pycirclize for chord diagram
chrom_labels = [f'chr{c}' for c in chrom_list]
chrom_colors_list = [CHROM_COLORS.get(c, '#999') for c in chrom_list]

# Scale sector sizes by number of coding mutations per chromosome
chrom_mut_counts = coding[coding['chrom'].isin(chrom_list)].groupby('chrom').size()
sector_sizes = {f'chr{c}': max(int(chrom_mut_counts.get(c, 1)), 1) for c in chrom_list}

circos = Circos(sector_sizes, space=5)

for sector in circos.sectors:
    chrom = sector.name.replace('chr', '')
    color = CHROM_COLORS.get(chrom, '#999')
    sector.axis(fc=color, ec='white', lw=1.5, alpha=0.8)
    sector.text(sector.name, r=125, fontsize=6, fontweight='bold')

# Add chord links for strong connections (>= 8 shared samples)
for i in range(n_chrom):
    for j in range(i + 1, n_chrom):
        val = adj_chrom[i, j]
        if val >= 8:
            alpha = min(val / 16, 0.7)
            lw = max(val / 4, 0.5)
            color1 = CHROM_COLORS.get(chrom_list[i], '#999')
            circos.link(
                (chrom_labels[i], 0, sector_sizes[chrom_labels[i]] * 0.8),
                (chrom_labels[j], 0, sector_sizes[chrom_labels[j]] * 0.8),
                color=color1, alpha=alpha * 0.6,
            )

fig = circos.plotfig()
fig.text(0.5, 1.06,
         'Chromosome-to-Chromosome Co-Mutation Links',
         ha='center', va='top', fontsize=14, fontweight='bold')
fig.text(0.5, 1.02,
         'Chords connect chromosomes sharing coding mutations in the same sample\n'
         '(HIGH/MODERATE impact, shown for pairs with shared mutations in 8+ samples)',
         ha='center', va='top', fontsize=9, color='gray')
fig.text(0.5, -0.02,
         'Sector size proportional to coding mutation count per chromosome',
         ha='center', va='top', fontsize=8, color='gray')

plt.savefig(f'{OUT_DIR}/22_chord_chromosome.png')
plt.close()
print("  -> 22_chord_chromosome.png")


# =============================================================================
# 2. CHORD DIAGRAM: Sample-to-Sample Shared Gene Mutations
# =============================================================================
print("Generating sample-to-sample chord diagram...")

# Adjacency: number of genes mutated in both samples
sample_genes = coding.groupby('short')['gene'].apply(set)
sample_list = sorted(sample_genes.index)
n_samples = len(sample_list)
adj_sample = np.zeros((n_samples, n_samples))

for i in range(n_samples):
    for j in range(i + 1, n_samples):
        shared = len(sample_genes[sample_list[i]] & sample_genes[sample_list[j]])
        adj_sample[i, j] = shared
        adj_sample[j, i] = shared

# Sector sizes by total coding mutations
sample_mut_counts = coding.groupby('short').size()
sector_sizes_s = {s: max(int(sample_mut_counts.get(s, 1)), 1) for s in sample_list}

# Sample colors
sample_palette = plt.cm.tab20(np.linspace(0, 1, n_samples))
sample_color_map = {s: matplotlib.colors.to_hex(sample_palette[i]) for i, s in enumerate(sample_list)}

circos2 = Circos(sector_sizes_s, space=6)

for sector in circos2.sectors:
    color = sample_color_map.get(sector.name, '#999')
    sector.axis(fc=color, ec='white', lw=1.5, alpha=0.85)
    sector.text(sector.name, r=125, fontsize=8, fontweight='bold')

# Add chords for top connections (>= 5 shared genes)
threshold = sorted(adj_sample[np.triu_indices(n_samples, k=1)], reverse=True)
if len(threshold) > 20:
    chord_thresh = max(threshold[20], 3)
else:
    chord_thresh = 3

for i in range(n_samples):
    for j in range(i + 1, n_samples):
        val = adj_sample[i, j]
        if val >= chord_thresh:
            alpha = min(val / adj_sample.max() * 0.8, 0.6)
            color = sample_color_map[sample_list[i]]
            circos2.link(
                (sample_list[i], 0, sector_sizes_s[sample_list[i]] * 0.8),
                (sample_list[j], 0, sector_sizes_s[sample_list[j]] * 0.8),
                color=color, alpha=alpha,
            )

fig = circos2.plotfig()
fig.text(0.5, 1.06,
         'Sample-to-Sample Co-Mutated Genes',
         ha='center', va='top', fontsize=14, fontweight='bold')
fig.text(0.5, 1.02,
         f'Chords connect samples sharing coding gene mutations\n'
         f'(HIGH/MODERATE impact, threshold: {int(chord_thresh)}+ shared genes)',
         ha='center', va='top', fontsize=9, color='gray')
fig.text(0.5, -0.02,
         'Sector size proportional to number of coding mutations per sample',
         ha='center', va='top', fontsize=8, color='gray')

plt.savefig(f'{OUT_DIR}/23_chord_samples.png')
plt.close()
print("  -> 23_chord_samples.png")


# =============================================================================
# 3. CHORD DIAGRAM: Pathway-to-Pathway Co-Mutation
# =============================================================================
print("Generating pathway co-mutation chord diagram...")

pathways = {
    'Chromatin\nRemodeling': ['Chd5', 'Kdm5a', 'Smarca4', 'Arid1a'],
    'Protein\nFolding': ['Pfdn2', 'Hsp90aa1', 'Cct2'],
    'Cell\nProliferation': ['Cdv3', 'Pimreg', 'Cdk4', 'Ccnd1'],
    'RNA\nSplicing': ['Prpf40a', 'Polr3d', 'Sf3b1', 'U2af1'],
    'DNA Damage\nResponse': ['Setx', 'Casp2', 'Brca1', 'Rad51'],
    'Metabolic\nRegulation': ['Helz2', 'Fasn', 'Acaca', 'Slc2a1'],
    'Immune\nSignaling': ['Il12rb2', 'Csf1r', 'Cd274', 'Jak2'],
    'EGFR/RTK\nSignaling': ['Rhbdf2', 'Egfr', 'Erbb2', 'Fgfr1'],
}

pathway_names = list(pathways.keys())
pathway_clean = [p.replace('\n', ' ') for p in pathway_names]
n_pw = len(pathway_names)

# Find which samples are affected by each pathway
pw_samples = {}
for pw_name, genes in pathways.items():
    affected = set()
    for gene in genes:
        gene_data = coding[coding['gene'] == gene]
        affected.update(gene_data['short'].unique())
    pw_samples[pw_name] = affected

# Adjacency: number of samples with mutations in BOTH pathways
adj_pw = np.zeros((n_pw, n_pw))
for i in range(n_pw):
    for j in range(i + 1, n_pw):
        shared = len(pw_samples[pathway_names[i]] & pw_samples[pathway_names[j]])
        adj_pw[i, j] = shared
        adj_pw[j, i] = shared

pw_colors = ['#8e44ad', '#2980b9', '#27ae60', '#e67e22',
             '#c0392b', '#16a085', '#d35400', '#7f8c8d']
pw_color_map = {pw_clean: pw_colors[i] for i, pw_clean in enumerate(pathway_clean)}

# Sector sizes
pw_sector_sizes = {pw_clean: max(len(pw_samples[pw_name]), 1)
                   for pw_name, pw_clean in zip(pathway_names, pathway_clean)}

circos3 = Circos(pw_sector_sizes, space=8)

for sector in circos3.sectors:
    color = pw_color_map.get(sector.name, '#999')
    sector.axis(fc=color, ec='white', lw=2, alpha=0.85)
    sector.text(sector.name, r=130, fontsize=7, fontweight='bold')

# Add all chords (pathways are few, show all connections)
for i in range(n_pw):
    for j in range(i + 1, n_pw):
        val = adj_pw[i, j]
        if val >= 1:
            alpha = min(val / 16, 0.6)
            color = pw_colors[i]
            circos3.link(
                (pathway_clean[i], 0, pw_sector_sizes[pathway_clean[i]] * 0.8),
                (pathway_clean[j], 0, pw_sector_sizes[pathway_clean[j]] * 0.8),
                color=color, alpha=max(alpha, 0.15),
            )

fig = circos3.plotfig()
fig.text(0.5, 1.06,
         'Pathway Co-Mutation Network',
         ha='center', va='top', fontsize=14, fontweight='bold')
fig.text(0.5, 1.02,
         'Chords connect pathways with shared affected samples\n'
         'Sector size = number of samples with mutations in that pathway',
         ha='center', va='top', fontsize=9, color='gray')

plt.savefig(f'{OUT_DIR}/24_chord_pathways.png')
plt.close()
print("  -> 24_chord_pathways.png")


# =============================================================================
# 4. CIRCA-STYLE PLOT: Genome-wide scatter with multi-track rings
# =============================================================================
print("\nGenerating circa-style genome-wide plot...")

fig, ax = plt.subplots(figsize=(18, 10))

# Draw chromosome blocks along x-axis
chrom_centers = {}
for i, c in enumerate(MAIN_CHROMS):
    start = cum_offset[c] / 1e6
    end = (cum_offset[c] + CHROM_SIZES[c]) / 1e6
    mid = (start + end) / 2
    chrom_centers[c] = mid
    color = CHROM_COLORS[c]
    # Alternating background
    if i % 2 == 0:
        ax.axvspan(start, end, alpha=0.06, color='gray')

# Plot variants as scatter points
# Y-axis = VAF, color = substitution class, size = impact
snvs = combined[(combined['var_type'] == 'SNV') & (combined['chrom'].isin(MAIN_CHROMS))].copy()
snvs['genome_pos'] = snvs.apply(lambda r: (cum_offset.get(r['chrom'], 0) + r['pos']) / 1e6, axis=1)
snvs['sub'] = snvs.apply(lambda r: get_sub_class(r['ref'], r['alt']), axis=1)

for sub in SUB_ORDER:
    mask = snvs['sub'] == sub
    ax.scatter(snvs.loc[mask, 'genome_pos'],
               snvs.loc[mask, 'af'],
               c=SUB_COLORS[sub], s=2, alpha=0.35, label=sub,
               rasterized=True, zorder=2)

# Overlay HIGH impact as larger markers
high = combined[(combined['impact'] == 'HIGH') & (combined['chrom'].isin(MAIN_CHROMS))].copy()
high['genome_pos'] = high.apply(lambda r: (cum_offset.get(r['chrom'], 0) + r['pos']) / 1e6, axis=1)
ax.scatter(high['genome_pos'], high['af'],
           c='none', edgecolors='#d62728', s=40, linewidths=1.2,
           alpha=0.8, label='HIGH impact', zorder=3)

# Annotate HIGH impact genes (only those with gene names, avoid overlap)
high_labeled = high[high['gene'] != ''].drop_duplicates(subset='gene').head(15)
for _, row in high_labeled.iterrows():
    ax.annotate(row['gene'],
                (row['genome_pos'], row['af']),
                fontsize=6, color='#d62728', alpha=0.8,
                xytext=(3, 5), textcoords='offset points',
                fontweight='bold')

# Chromosome labels at bottom
for c in MAIN_CHROMS:
    ax.text(chrom_centers[c], -0.07, c, ha='center', va='top',
            fontsize=8, fontweight='bold', color=CHROM_COLORS[c])

# Chromosome boundaries
for c in MAIN_CHROMS:
    bnd = cum_offset[c] / 1e6
    ax.axvline(x=bnd, color='gray', linewidth=0.3, alpha=0.3, zorder=1)

ax.set_xlim(0, GENOME_TOTAL / 1e6)
ax.set_ylim(-0.02, 1.02)
ax.set_xlabel('Chromosome', fontsize=12, labelpad=15)
ax.set_ylabel('Variant Allele Frequency (VAF)', fontsize=12)
ax.set_title('Circa Plot: Genome-Wide Somatic Variant Distribution\n'
             'Filtered variants colored by substitution type, HIGH impact circled in red',
             fontsize=14, fontweight='bold', pad=15)
ax.set_xticks([])
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)

ax.legend(title='Substitution', bbox_to_anchor=(1.01, 1), loc='upper left',
          fontsize=8, title_fontsize=9, markerscale=4, frameon=True)

plt.tight_layout()
plt.savefig(f'{OUT_DIR}/25_circa_genomewide.png')
plt.close()
print("  -> 25_circa_genomewide.png")


# =============================================================================
# 5. CIRCA-STYLE PLOT: Per-sample comparison (2x2)
# =============================================================================
print("Generating per-sample circa plots (2x2)...")

select_shorts = ['S15', 'S04', 'S16', 'S01']
select_samples = []
for short in select_shorts:
    for s in SAMPLES:
        if SAMPLE_SHORT[s] == short:
            select_samples.append(s)
            break

fig, axes = plt.subplots(2, 2, figsize=(20, 12))
fig.suptitle('Circa Plots: Per-Sample Genome-Wide Variant Distribution',
             fontsize=15, fontweight='bold', y=1.0)

for ax, sample in zip(axes.flat, select_samples):
    short = SAMPLE_SHORT[sample]
    df = all_variants[sample]
    df_main = df[df['chrom'].isin(MAIN_CHROMS)].copy()
    df_main['genome_pos'] = df_main.apply(
        lambda r: (cum_offset.get(r['chrom'], 0) + r['pos']) / 1e6, axis=1)

    # Background bands
    for i, c in enumerate(MAIN_CHROMS):
        start = cum_offset[c] / 1e6
        end = (cum_offset[c] + CHROM_SIZES[c]) / 1e6
        if i % 2 == 0:
            ax.axvspan(start, end, alpha=0.06, color='gray')
        ax.axvline(x=start, color='gray', linewidth=0.2, alpha=0.3)

    # SNVs
    snvs = df_main[df_main['var_type'] == 'SNV'].copy()
    snvs['sub'] = snvs.apply(lambda r: get_sub_class(r['ref'], r['alt']), axis=1)
    for sub in SUB_ORDER:
        mask = snvs['sub'] == sub
        ax.scatter(snvs.loc[mask, 'genome_pos'], snvs.loc[mask, 'af'],
                   c=SUB_COLORS[sub], s=1.5, alpha=0.4, rasterized=True)

    # Indels as small triangles on bottom
    indels = df_main[df_main['var_type'].isin(['Insertion', 'Deletion'])]
    ax.scatter(indels['genome_pos'], indels['af'],
               c='#3498db', s=1, alpha=0.2, marker='^', rasterized=True)

    # HIGH impact
    hi = df_main[df_main['impact'] == 'HIGH']
    ax.scatter(hi['genome_pos'], hi['af'],
               c='none', edgecolors='#d62728', s=35, linewidths=1,
               alpha=0.9, zorder=3)

    tmb = len(df) / GENOME_SIZE_MB
    indel_pct = df['var_type'].isin(['Insertion', 'Deletion']).sum() / len(df) * 100
    n_high = (df['impact'] == 'HIGH').sum()
    ax.set_title(f'{short}  (n={len(df):,}  TMB={tmb:.1f}  '
                 f'Indel={indel_pct:.0f}%  HIGH={n_high})',
                 fontsize=11, fontweight='bold', pad=8)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlim(0, GENOME_TOTAL / 1e6)
    ax.set_ylabel('VAF', fontsize=9)
    ax.set_xticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    # Chromosome labels
    for c in MAIN_CHROMS:
        mid = (cum_offset[c] + CHROM_SIZES[c] / 2) / 1e6
        label = c if c not in ('Y',) else 'Y'
        ax.text(mid, -0.06, label, ha='center', va='top',
                fontsize=5.5, color='gray')

# Legend at bottom
legend_elements = [
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=SUB_COLORS[s],
               markersize=5, label=s) for s in SUB_ORDER
] + [
    plt.Line2D([0], [0], marker='^', color='w', markerfacecolor='#3498db',
               markersize=5, label='Indel'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='none',
               markeredgecolor='#d62728', markersize=6, markeredgewidth=1.2,
               label='HIGH impact'),
]
fig.legend(handles=legend_elements, loc='lower center', ncol=8, fontsize=8,
           bbox_to_anchor=(0.5, -0.03), frameon=True, title='Variant Type')

plt.tight_layout(rect=[0, 0.03, 1, 0.97])
plt.savefig(f'{OUT_DIR}/26_circa_per_sample.png')
plt.close()
print("  -> 26_circa_per_sample.png")


# =============================================================================
# 6. CIRCA-STYLE PLOT: Variant density + impact rings (circular)
# =============================================================================
print("Generating circa-style circular plot with density rings...")

sectors = {f'chr{c}': CHROM_SIZES[c] for c in MAIN_CHROMS}
circos4 = Circos(sectors, space=4)

for sector in circos4.sectors:
    chrom = sector.name.replace('chr', '')
    chrom_size = CHROM_SIZES[chrom]
    color = CHROM_COLORS[chrom]
    n_win = int(np.ceil(chrom_size / WINDOW_SIZE))
    positions = np.linspace(0, chrom_size, n_win + 1)
    x_mid = np.array([(positions[i] + positions[i+1]) / 2 for i in range(n_win)])

    # Outer axis with chromosome color
    sector.axis(fc=color, ec='white', lw=1.2, alpha=0.7)
    sector.text(sector.name, r=125, fontsize=6, fontweight='bold')

    # Track 1: SNV density (filled area)
    snv_data = combined[(combined['chrom'] == chrom) & (combined['var_type'] == 'SNV')]
    snv_density = np.zeros(n_win)
    for pos in snv_data['pos']:
        idx = min(int(pos / WINDOW_SIZE), n_win - 1)
        snv_density[idx] += 1

    track1 = sector.add_track((85, 97))
    track1.axis(fc='none', ec='gray', lw=0.3)
    if snv_density.max() > 0:
        snv_norm = snv_density / snv_density.max() * (97 - 85)
        track1.fill_between(x_mid, snv_norm + 85, 85, fc='#E62725', alpha=0.5, ec='none')

    # Track 2: Indel density
    indel_data = combined[(combined['chrom'] == chrom) &
                          (combined['var_type'].isin(['Insertion', 'Deletion']))]
    indel_density = np.zeros(n_win)
    for pos in indel_data['pos']:
        idx = min(int(pos / WINDOW_SIZE), n_win - 1)
        indel_density[idx] += 1

    track2 = sector.add_track((72, 84))
    track2.axis(fc='none', ec='gray', lw=0.3)
    if indel_density.max() > 0:
        indel_norm = indel_density / max(indel_density.max(), 1) * (84 - 72)
        track2.fill_between(x_mid, indel_norm + 72, 72, fc='#3498db', alpha=0.5, ec='none')

    # Track 3: VAF scatter (subsample for performance)
    chrom_vars = combined[(combined['chrom'] == chrom)].copy()
    if len(chrom_vars) > 500:
        chrom_vars = chrom_vars.sample(500, random_state=42)

    track3 = sector.add_track((55, 71))
    track3.axis(fc='none', ec='gray', lw=0.3)

    # Color by impact
    for impact, color_imp in [('HIGH', '#d62728'), ('MODERATE', '#ff7f0e')]:
        imp_data = chrom_vars[chrom_vars['impact'] == impact]
        if len(imp_data) > 0:
            vaf_scaled = imp_data['af'].values * (71 - 55) + 55
            track3.scatter(imp_data['pos'].values, vaf_scaled,
                          s=8 if impact == 'HIGH' else 3,
                          c=color_imp, alpha=0.7, zorder=3)

    # Others
    other = chrom_vars[~chrom_vars['impact'].isin(['HIGH', 'MODERATE'])]
    if len(other) > 0:
        vaf_scaled = other['af'].values * (71 - 55) + 55
        track3.scatter(other['pos'].values, vaf_scaled, s=1,
                      c='#aec7e8', alpha=0.3)

fig = circos4.plotfig()
fig.text(0.5, 1.06,
         'Circa Plot: Genome-Wide Multi-Track Variant View',
         ha='center', va='top', fontsize=14, fontweight='bold')
fig.text(0.5, 1.02,
         'Outer: SNV density | Middle: Indel density | Inner: VAF scatter (colored by impact)',
         ha='center', va='top', fontsize=9, color='gray')

legend_elements = [
    Patch(fc='#E62725', alpha=0.5, label='SNV density'),
    Patch(fc='#3498db', alpha=0.5, label='Indel density'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#d62728',
               markersize=6, label='HIGH impact'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#ff7f0e',
               markersize=5, label='MODERATE impact'),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#aec7e8',
               markersize=4, label='LOW/MODIFIER'),
]
fig.legend(handles=legend_elements, loc='lower center', ncol=5,
           fontsize=8, bbox_to_anchor=(0.5, 0.01), frameon=True)

plt.savefig(f'{OUT_DIR}/27_circa_circular.png')
plt.close()
print("  -> 27_circa_circular.png")


print("\nAll chord and circa plots complete!")
