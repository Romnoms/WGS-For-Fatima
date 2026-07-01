#!/usr/bin/env python3
"""
Sucrose-Context Analysis — Additional figures for sucrose-induced tumor interpretation.
C57BL/6J mice, ad libitum sucrose water. Plain water controls were tumor-free.
GRCm39 (JAX C57BL/6J) serves as germline reference.
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
import warnings
warnings.filterwarnings('ignore')

plt.rcParams.update({
    'font.size': 11,
    'font.family': 'sans-serif',
    'axes.titlesize': 13,
    'axes.labelsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
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

SUB_COLORS = {'C>A': '#3498db', 'C>G': '#2c3e50', 'C>T': '#e74c3c',
              'T>A': '#bdc3c7', 'T>C': '#27ae60', 'T>G': '#f39c12'}


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

            indel_size = len(alt) - len(ref) if var_type != 'SNV' else 0

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
                'indel_size': indel_size,
                'gene': gene, 'consequence': consequence, 'impact': impact
            })
    return pd.DataFrame(variants)


# Parse all filtered VCFs
print("Parsing filtered somatic VCFs...")
all_variants = {}
for sample in SAMPLES:
    vcf_path = f'{DATA_DIR}/{sample}.filtered.vcf.gz'
    df = parse_filtered_vcf(vcf_path)
    df['sample'] = sample
    df['short'] = SAMPLE_SHORT[sample]
    all_variants[sample] = df

combined = pd.concat(all_variants.values(), ignore_index=True)


# =============================================================================
# Figure 1: SNV:Indel ratio comparison — highlights MSI-like phenotype
# =============================================================================
print("Generating SNV:Indel ratio figure...")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))
fig.suptitle('SNV vs Indel Composition — Evidence of Replication Slippage',
             fontsize=14, fontweight='bold')

snv_counts = []
indel_counts = []
indel_pcts = []
for s in SAMPLES:
    df = all_variants[s]
    n_snv = (df['var_type'] == 'SNV').sum()
    n_indel = df['var_type'].isin(['Insertion', 'Deletion']).sum()
    snv_counts.append(n_snv)
    indel_counts.append(n_indel)
    indel_pcts.append(n_indel / (n_snv + n_indel) * 100)

x = np.arange(len(LABELS))
width = 0.35
ax1.bar(x - width/2, snv_counts, width, label='SNVs', color='#3498db', edgecolor='white')
ax1.bar(x + width/2, indel_counts, width, label='Indels', color='#e74c3c', edgecolor='white')
ax1.set_xticks(x)
ax1.set_xticklabels(LABELS, rotation=45, ha='right')
ax1.set_ylabel('Variant Count')
ax1.set_title('SNV vs Indel Counts')
ax1.legend()
ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: f'{int(y):,}'))

colors = ['#e74c3c' if p > 60 else '#f39c12' if p > 50 else '#3498db' for p in indel_pcts]
ax2.bar(x, indel_pcts, color=colors, edgecolor='white')
ax2.axhline(y=50, color='gray', linestyle='--', alpha=0.7, label='50% (equal SNV:indel)')
ax2.set_xticks(x)
ax2.set_xticklabels(LABELS, rotation=45, ha='right')
ax2.set_ylabel('Indel Fraction (%)')
ax2.set_title('Indel Proportion — MSI-like Phenotype')
ax2.legend(fontsize=9)

# Annotate hypermutated
for i, (label, pct) in enumerate(zip(LABELS, indel_pcts)):
    if pct > 60:
        ax2.annotate(f'{pct:.0f}%', (i, pct + 1), ha='center', fontsize=9, fontweight='bold',
                     color='#e74c3c')

plt.tight_layout()
plt.savefig(f'{OUT_DIR}/10_snv_indel_ratio.png')
plt.close()
print("  -> 10_snv_indel_ratio.png")


# =============================================================================
# Figure 2: Indel size distribution — short indels = replication slippage
# =============================================================================
print("Generating indel size distribution...")
fig, axes = plt.subplots(2, 2, figsize=(13, 9))
fig.suptitle('Indel Size Distribution — Replication Slippage Signature',
             fontsize=14, fontweight='bold')

# Representative samples: typical (S01), elevated (S16), hypermutated (S04, S15)
plot_samples = [('01', 'S01 (Typical)'), ('16', 'S16 (Elevated)'),
                ('04', 'S04 (Hypermutated)'), ('15', 'S15 (Hypermutated)')]

for ax, (sid, title) in zip(axes.flat, plot_samples):
    sample = f'26034XD-04-{sid}'
    df = all_variants[sample]
    indels = df[df['var_type'].isin(['Insertion', 'Deletion'])]
    sizes = indels['indel_size'].values
    sizes = sizes[(sizes >= -30) & (sizes <= 30) & (sizes != 0)]

    ins = sizes[sizes > 0]
    dels = sizes[sizes < 0]

    bins_ins = np.arange(0.5, 31.5, 1)
    bins_del = np.arange(-30.5, 0.5, 1)

    ax.hist(dels, bins=bins_del, color='#e74c3c', alpha=0.8, label='Deletions', edgecolor='white')
    ax.hist(ins, bins=bins_ins, color='#3498db', alpha=0.8, label='Insertions', edgecolor='white')
    ax.set_title(title, fontweight='bold')
    ax.set_xlabel('Indel Size (bp)')
    ax.set_ylabel('Count')
    ax.legend(fontsize=8)
    ax.axvline(x=0, color='black', linewidth=0.5)

    # Annotate 1bp fraction
    one_bp = ((sizes == 1) | (sizes == -1)).sum()
    total = len(sizes)
    if total > 0:
        ax.text(0.98, 0.95, f'1bp indels: {one_bp/total*100:.0f}%',
                transform=ax.transAxes, ha='right', va='top', fontsize=9,
                bbox=dict(boxstyle='round,pad=0.3', facecolor='wheat', alpha=0.8))

plt.tight_layout()
plt.savefig(f'{OUT_DIR}/11_indel_size_distribution.png')
plt.close()
print("  -> 11_indel_size_distribution.png")


# =============================================================================
# Figure 3: VAF vs Depth — somatic variant confidence
# =============================================================================
print("Generating VAF vs depth scatter...")
fig, axes = plt.subplots(2, 2, figsize=(13, 10))
fig.suptitle('Variant Allele Frequency vs Read Depth — Somatic Variant Validation',
             fontsize=14, fontweight='bold')

for ax, (sid, title) in zip(axes.flat, plot_samples):
    sample = f'26034XD-04-{sid}'
    df = all_variants[sample]
    snvs = df[df['var_type'] == 'SNV']
    indels = df[df['var_type'].isin(['Insertion', 'Deletion'])]

    ax.scatter(snvs['af'], snvs['dp'], alpha=0.15, s=3, color='#3498db', label='SNVs')
    ax.scatter(indels['af'], indels['dp'], alpha=0.15, s=3, color='#e74c3c', label='Indels')
    ax.set_title(title, fontweight='bold')
    ax.set_xlabel('VAF')
    ax.set_ylabel('Read Depth')
    ax.set_xlim(0, 1.05)
    ax.set_ylim(0, min(df['dp'].quantile(0.99) * 1.2, 200))
    ax.legend(fontsize=8, markerscale=5)

    # Add median VAF annotation
    med_vaf = df['af'].median()
    ax.axvline(x=med_vaf, color='gray', linestyle='--', alpha=0.5)
    ax.text(med_vaf + 0.02, ax.get_ylim()[1] * 0.9, f'median={med_vaf:.2f}',
            fontsize=8, color='gray')

plt.tight_layout()
plt.savefig(f'{OUT_DIR}/12_vaf_vs_depth.png')
plt.close()
print("  -> 12_vaf_vs_depth.png")


# =============================================================================
# Figure 4: Mutation process summary — linking to sucrose exposure
# =============================================================================
print("Generating mutation process summary...")

# Calculate per-sample metrics
process_data = []
for s in SAMPLES:
    df = all_variants[s]
    short = SAMPLE_SHORT[s]
    total = len(df)
    snvs = df[df['var_type'] == 'SNV']
    indels = df[df['var_type'].isin(['Insertion', 'Deletion'])]

    # Substitution spectrum
    spec = Counter()
    for _, row in snvs.iterrows():
        ref, alt = row['ref'], row['alt']
        if ref in ('C', 'T'):
            spec[f'{ref}>{alt}'] += 1
        else:
            spec[f'{COMPLEMENT[ref]}>{COMPLEMENT[alt]}'] += 1

    total_snv = sum(spec.values())
    ct_pct = spec.get('C>T', 0) / total_snv * 100 if total_snv > 0 else 0
    tg_pct = spec.get('T>G', 0) / total_snv * 100 if total_snv > 0 else 0
    tc_pct = spec.get('T>C', 0) / total_snv * 100 if total_snv > 0 else 0

    # Indel metrics
    indel_pct = len(indels) / total * 100 if total > 0 else 0
    small_indels = indels[indels['indel_size'].abs() <= 5]
    small_indel_pct = len(small_indels) / len(indels) * 100 if len(indels) > 0 else 0

    tmb = total / GENOME_SIZE_MB

    process_data.append({
        'Sample': short,
        'TMB': tmb,
        'C>T %': ct_pct,
        'T>G %': tg_pct,
        'T>C %': tc_pct,
        'Indel %': indel_pct,
        'Small Indel %': small_indel_pct,
    })

pdf = pd.DataFrame(process_data)

fig = plt.figure(figsize=(16, 10))
gs = GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.3)
fig.suptitle('Mutational Process Summary — Sucrose-Induced Tumors',
             fontsize=15, fontweight='bold', y=0.98)

# TMB
ax = fig.add_subplot(gs[0, 0])
colors = ['#e74c3c' if t > 4 else '#3498db' for t in pdf['TMB']]
ax.bar(range(len(pdf)), pdf['TMB'], color=colors, edgecolor='white')
ax.set_xticks(range(len(pdf)))
ax.set_xticklabels(pdf['Sample'], rotation=45, ha='right', fontsize=8)
ax.set_ylabel('mut/Mb')
ax.set_title('Tumor Mutational Burden')

# C>T (deamination / cell division clock)
ax = fig.add_subplot(gs[0, 1])
ax.bar(range(len(pdf)), pdf['C>T %'], color='#e74c3c', edgecolor='white')
ax.set_xticks(range(len(pdf)))
ax.set_xticklabels(pdf['Sample'], rotation=45, ha='right', fontsize=8)
ax.set_ylabel('%')
ax.set_title('C>T Transitions (SBS1 — Cell Divisions)')
ax.axhline(y=pdf['C>T %'].median(), color='gray', linestyle='--', alpha=0.5)

# T>G (oxidative damage)
ax = fig.add_subplot(gs[0, 2])
colors = ['#e74c3c' if t > 15 else '#f39c12' for t in pdf['T>G %']]
ax.bar(range(len(pdf)), pdf['T>G %'], color=colors, edgecolor='white')
ax.set_xticks(range(len(pdf)))
ax.set_xticklabels(pdf['Sample'], rotation=45, ha='right', fontsize=8)
ax.set_ylabel('%')
ax.set_title('T>G Transversions (Oxidative Damage)')
ax.axhline(y=pdf['T>G %'].median(), color='gray', linestyle='--', alpha=0.5)

# Indel fraction (replication slippage)
ax = fig.add_subplot(gs[1, 0])
colors = ['#e74c3c' if p > 60 else '#f39c12' if p > 50 else '#3498db' for p in pdf['Indel %']]
ax.bar(range(len(pdf)), pdf['Indel %'], color=colors, edgecolor='white')
ax.set_xticks(range(len(pdf)))
ax.set_xticklabels(pdf['Sample'], rotation=45, ha='right', fontsize=8)
ax.set_ylabel('%')
ax.set_title('Indel Fraction (Replication Slippage)')
ax.axhline(y=50, color='gray', linestyle='--', alpha=0.5)

# Small indel fraction of all indels
ax = fig.add_subplot(gs[1, 1])
ax.bar(range(len(pdf)), pdf['Small Indel %'], color='#9b59b6', edgecolor='white')
ax.set_xticks(range(len(pdf)))
ax.set_xticklabels(pdf['Sample'], rotation=45, ha='right', fontsize=8)
ax.set_ylabel('%')
ax.set_title('Short Indels (<=5bp) of All Indels')

# Correlation: TMB vs Indel %
ax = fig.add_subplot(gs[1, 2])
scatter_colors = ['#e74c3c' if p > 60 else '#3498db' for p in pdf['Indel %']]
ax.scatter(pdf['TMB'], pdf['Indel %'], c=scatter_colors, s=80, edgecolors='white', zorder=3)
for _, row in pdf.iterrows():
    ax.annotate(row['Sample'], (row['TMB'], row['Indel %']),
                fontsize=7, ha='center', va='bottom', xytext=(0, 5),
                textcoords='offset points')
ax.set_xlabel('TMB (mut/Mb)')
ax.set_ylabel('Indel Fraction (%)')
ax.set_title('TMB vs Indel Fraction')
ax.axhline(y=50, color='gray', linestyle='--', alpha=0.3)

plt.savefig(f'{OUT_DIR}/13_mutation_process_summary.png')
plt.close()
print("  -> 13_mutation_process_summary.png")


# =============================================================================
# Figure 5: Gene impact landscape — oncoplot-style heatmap
# =============================================================================
print("Generating gene impact landscape...")

coding = combined[(combined['impact'].isin(['HIGH', 'MODERATE'])) & (combined['gene'] != '')].copy()

# Get genes mutated in >= 2 samples
gene_sample_counts = coding.groupby('gene')['short'].apply(lambda x: set(x))
genes_multi = {g: s for g, s in gene_sample_counts.items() if len(s) >= 2}

# Filter out Vmn/Zfp/Or families (likely residual germline noise)
skip_prefixes = ('Vmn', 'Zfp', 'Or', 'Gm')
genes_filtered = {g: s for g, s in genes_multi.items()
                  if not any(g.startswith(p) for p in skip_prefixes)}

# Sort by recurrence
genes_sorted = sorted(genes_filtered.items(), key=lambda x: len(x[1]), reverse=True)
if len(genes_sorted) > 20:
    genes_sorted = genes_sorted[:20]

gene_names = [g for g, _ in genes_sorted]

# Build matrix: gene x sample, value = impact severity
impact_val = {'HIGH': 3, 'MODERATE': 2, 'LOW': 1, 'MODIFIER': 0}
matrix = np.zeros((len(gene_names), len(LABELS)))

for i, gene in enumerate(gene_names):
    gene_data = coding[coding['gene'] == gene]
    for _, row in gene_data.iterrows():
        j = LABELS.index(row['short'])
        val = impact_val.get(row['impact'], 0)
        matrix[i, j] = max(matrix[i, j], val)

from matplotlib.colors import ListedColormap, BoundaryNorm
cmap = ListedColormap(['#f0f0f0', '#2ca02c', '#ff7f0e', '#d62728'])
bounds = [-0.5, 0.5, 1.5, 2.5, 3.5]
norm = BoundaryNorm(bounds, cmap.N)

fig, ax = plt.subplots(figsize=(14, max(6, len(gene_names) * 0.4 + 2)))
im = ax.imshow(matrix, aspect='auto', cmap=cmap, norm=norm)
ax.set_xticks(range(len(LABELS)))
ax.set_xticklabels(LABELS, rotation=45, ha='right')
ax.set_yticks(range(len(gene_names)))
ax.set_yticklabels(gene_names, fontsize=10)
ax.set_title('Somatic Gene Mutation Landscape — Sucrose-Induced Tumors\n'
             '(HIGH/MODERATE impact, recurrent in >=2 samples, excluding Vmn/Zfp/Or/Gm families)',
             fontsize=12, fontweight='bold')

# Add recurrence count on right
ax2 = ax.twinx()
ax2.set_ylim(ax.get_ylim())
ax2.set_yticks(range(len(gene_names)))
recurrence = [len(genes_filtered[g]) for g in gene_names]
ax2.set_yticklabels([f'{r}/16' for r in recurrence], fontsize=9)
ax2.set_ylabel('Samples Mutated', fontsize=10)

# Legend
legend_elements = [
    Patch(facecolor='#f0f0f0', edgecolor='gray', label='Wild-type'),
    Patch(facecolor='#ff7f0e', label='MODERATE'),
    Patch(facecolor='#d62728', label='HIGH'),
]
ax.legend(handles=legend_elements, loc='lower right', fontsize=9, title='Impact')

# Grid
for i in range(len(gene_names) + 1):
    ax.axhline(y=i - 0.5, color='white', linewidth=1)
for j in range(len(LABELS) + 1):
    ax.axvline(x=j - 0.5, color='white', linewidth=1)

plt.tight_layout()
plt.savefig(f'{OUT_DIR}/14_gene_landscape.png')
plt.close()
print("  -> 14_gene_landscape.png")


# =============================================================================
# Figure 6: Sucrose-relevant pathway mutations
# =============================================================================
print("Generating pathway-level figure...")

pathways = {
    'Chromatin\nRemodeling': {'genes': ['Chd5'], 'color': '#8e44ad'},
    'Protein\nFolding': {'genes': ['Pfdn2'], 'color': '#2980b9'},
    'Cell\nProliferation': {'genes': ['Cdv3', 'Pimreg'], 'color': '#27ae60'},
    'RNA\nSplicing': {'genes': ['Prpf40a', 'Polr3d'], 'color': '#e67e22'},
    'DNA Damage\nResponse': {'genes': ['Setx', 'Casp2'], 'color': '#c0392b'},
    'Metabolic\nRegulation': {'genes': ['Helz2', 'Cdv3'], 'color': '#16a085'},
    'Immune\nSignaling': {'genes': ['Il12rb2', 'Csf1r'], 'color': '#d35400'},
    'EGFR/\nInflammation': {'genes': ['Rhbdf2'], 'color': '#7f8c8d'},
}

pathway_names = list(pathways.keys())
pathway_sample_counts = []
pathway_gene_lists = []

for pw_name, pw_info in pathways.items():
    samples_affected = set()
    for gene in pw_info['genes']:
        gene_data = coding[coding['gene'] == gene]
        samples_affected.update(gene_data['short'].unique())
    pathway_sample_counts.append(len(samples_affected))
    pathway_gene_lists.append(', '.join(pw_info['genes']))

fig, ax = plt.subplots(figsize=(10, 6))
colors = [pathways[p]['color'] for p in pathway_names]
bars = ax.barh(range(len(pathway_names)), pathway_sample_counts, color=colors, edgecolor='white')
ax.set_yticks(range(len(pathway_names)))
ax.set_yticklabels(pathway_names, fontsize=10)
ax.set_xlabel('Number of Samples with Mutations')
ax.set_title('Pathway-Level Somatic Mutations in Sucrose-Induced Tumors',
             fontsize=13, fontweight='bold')
ax.invert_yaxis()
ax.set_xlim(0, 16)

# Annotate with gene names
for i, (count, genes) in enumerate(zip(pathway_sample_counts, pathway_gene_lists)):
    ax.text(count + 0.3, i, f'{genes}  ({count}/16)', va='center', fontsize=8, color='#555')

# Add vertical line at n=16
ax.axvline(x=16, color='gray', linestyle=':', alpha=0.3)

plt.tight_layout()
plt.savefig(f'{OUT_DIR}/15_pathway_mutations.png')
plt.close()
print("  -> 15_pathway_mutations.png")


# =============================================================================
# Figure 7: Germline reference validation — JAX B6J context
# =============================================================================
print("Generating reference validation figure...")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))
fig.suptitle('Somatic Variant Validation — C57BL/6J (JAX) on GRCm39 Reference',
             fontsize=14, fontweight='bold')

# VAF distribution comparison: show the shift from unfiltered to filtered
# Unfiltered data
UNFILT_DIR = os.path.expanduser('~/wgs-project/data/05.SomaticSNV')
sample_01_unfilt = f'{UNFILT_DIR}/26034XD-04-01.somatic.final.vep.vcf.gz'
unfilt_vafs = []
with gzip.open(sample_01_unfilt, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        parts = line.strip().split('\t')
        fmt_keys = parts[8].split(':')
        fmt_vals = parts[9].split(':')
        fmt_dict = dict(zip(fmt_keys, fmt_vals))
        unfilt_vafs.append(float(fmt_dict.get('AF', '0')))

filt_vafs = all_variants['26034XD-04-01']['af'].values

ax1.hist(unfilt_vafs, bins=50, range=(0, 1), alpha=0.6, color='#bdc3c7',
         label=f'Unfiltered (n={len(unfilt_vafs):,})', edgecolor='white')
ax1.hist(filt_vafs, bins=50, range=(0, 1), alpha=0.7, color='#2ecc71',
         label=f'Filtered (n={len(filt_vafs):,})', edgecolor='white')
ax1.set_xlabel('Variant Allele Frequency')
ax1.set_ylabel('Count')
ax1.set_title('S01: VAF Before vs After Filtering')
ax1.legend()
ax1.annotate('Germline peak\n(removed)', xy=(0.5, 3500), fontsize=9,
             ha='center', color='#7f8c8d',
             arrowprops=dict(arrowstyle='->', color='#7f8c8d'),
             xytext=(0.7, 4500))

# Somatic confidence: since B6J = GRCm39, show what was removed
removed_count = len(unfilt_vafs) - len(filt_vafs)
labels_pie = [f'Removed as\ngermline/artifact\n({removed_count:,})',
              f'Retained as\nsomatic\n({len(filt_vafs):,})']
sizes = [removed_count, len(filt_vafs)]
colors_pie = ['#bdc3c7', '#2ecc71']
explode = (0, 0.05)
wedges, texts, autotexts = ax2.pie(sizes, labels=labels_pie, colors=colors_pie,
                                    explode=explode, autopct='%1.1f%%',
                                    startangle=90, textprops={'fontsize': 10})
autotexts[0].set_fontweight('bold')
autotexts[1].set_fontweight('bold')
ax2.set_title('S01: Variant Classification\n(GRCm39 = JAX C57BL/6J germline)')

plt.tight_layout()
plt.savefig(f'{OUT_DIR}/16_reference_validation.png')
plt.close()
print("  -> 16_reference_validation.png")


print(f"\nAll new figures saved to {OUT_DIR}/")
print("Done!")
