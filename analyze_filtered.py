#!/usr/bin/env python3
"""
WGS Filtered Analysis — Sucrose-Induced Tumors in C57BL/6J Mice (GRCm39)
Generates figures from filtered somatic variant calls
(post-MGP and cross-sample germline subtraction).
Plain water controls remained tumor-free; GRCm39 from JAX C57BL/6J
serves as the germline reference.
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
UNFILT_DIR = os.path.expanduser('~/wgs-project/data/05.SomaticSNV')
OUT_DIR = os.path.expanduser('~/wgs-project/figures_filtered')
os.makedirs(OUT_DIR, exist_ok=True)

SAMPLES = sorted([
    os.path.basename(f).replace('.filtered.vcf.gz', '')
    for f in glob.glob(f'{DATA_DIR}/*.filtered.vcf.gz')
])
SAMPLE_SHORT = {s: s.replace('26034XD-04-', 'S') for s in SAMPLES}
GENOME_SIZE_MB = 2728.22

IMPACT_COLORS = {'HIGH': '#d62728', 'MODERATE': '#ff7f0e', 'LOW': '#2ca02c', 'MODIFIER': '#aec7e8'}
SUB_COLORS = {'C>A': '#3498db', 'C>G': '#2c3e50', 'C>T': '#e74c3c',
              'T>A': '#bdc3c7', 'T>C': '#27ae60', 'T>G': '#f39c12'}
COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def parse_vcf(vcf_path):
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
    df = parse_vcf(vcf_path)
    df['sample'] = sample
    df['short'] = SAMPLE_SHORT[sample]
    all_variants[sample] = df
    print(f"  {SAMPLE_SHORT[sample]}: {len(df)} filtered variants")

combined = pd.concat(all_variants.values(), ignore_index=True)

# Also get unfiltered counts for comparison
print("\nCounting unfiltered variants for comparison...")
unfilt_counts = {}
for sample in SAMPLES:
    vcf_path = f'{UNFILT_DIR}/{sample}.somatic.final.vep.vcf.gz'
    count = sum(1 for line in gzip.open(vcf_path, 'rt') if not line.startswith('#'))
    unfilt_counts[sample] = count


# =============================================================================
# 1. Before/After Filtering Comparison
# =============================================================================
print("Generating filtering comparison figure...")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
fig.suptitle('Impact of Germline Filtering on Somatic Variant Calls', fontsize=15, fontweight='bold')

x = range(len(SAMPLES))
labels = [SAMPLE_SHORT[s] for s in SAMPLES]
unfilt_vals = [unfilt_counts[s] for s in SAMPLES]
filt_vals = [len(all_variants[s]) for s in SAMPLES]

ax1.bar(x, unfilt_vals, color='#bdc3c7', label='Unfiltered', edgecolor='white')
ax1.bar(x, filt_vals, color='#2ecc71', label='Filtered', edgecolor='white')
ax1.set_xticks(x)
ax1.set_xticklabels(labels, rotation=45, ha='right')
ax1.set_ylabel('Variant Count')
ax1.set_title('Variant Counts: Before vs After Filtering')
ax1.legend()
ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: f'{int(y):,}'))

# Percentage retained
pct_retained = [f / u * 100 for f, u in zip(filt_vals, unfilt_vals)]
colors = ['#e74c3c' if p < 8 else '#3498db' for p in pct_retained]
ax2.bar(x, pct_retained, color=colors, edgecolor='white')
ax2.set_xticks(x)
ax2.set_xticklabels(labels, rotation=45, ha='right')
ax2.set_ylabel('Variants Retained (%)')
ax2.set_title('Percentage Retained After Filtering')
ax2.axhline(y=np.mean(pct_retained), color='gray', linestyle='--', alpha=0.7,
            label=f'Mean: {np.mean(pct_retained):.1f}%')
ax2.legend()

plt.tight_layout()
plt.savefig(f'{OUT_DIR}/01_filtering_comparison.png')
plt.close()
print("  -> 01_filtering_comparison.png")


# =============================================================================
# 2. Filtered Mutation Burden
# =============================================================================
print("Generating filtered mutation burden...")
burden = combined.groupby('short').size().reindex(labels)
burden_per_mb = burden / GENOME_SIZE_MB

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
fig.suptitle('Somatic Mutation Burden (After Germline Filtering)', fontsize=15, fontweight='bold')

colors = ['#e74c3c' if c > 10000 else '#3498db' for c in burden.values]
ax1.bar(range(len(burden)), burden.values, color=colors, edgecolor='white')
ax1.set_xticks(range(len(burden)))
ax1.set_xticklabels(burden.index, rotation=45, ha='right')
ax1.set_ylabel('Filtered Somatic Variants')
ax1.set_title('Total Variant Count')
ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: f'{int(y):,}'))

colors2 = ['#e74c3c' if c > 4 else '#3498db' for c in burden_per_mb.values]
ax2.bar(range(len(burden_per_mb)), burden_per_mb.values, color=colors2, edgecolor='white')
ax2.set_xticks(range(len(burden_per_mb)))
ax2.set_xticklabels(burden_per_mb.index, rotation=45, ha='right')
ax2.set_ylabel('Mutations / Mb')
ax2.set_title('Tumor Mutational Burden (TMB)')

plt.tight_layout()
plt.savefig(f'{OUT_DIR}/02_filtered_mutation_burden.png')
plt.close()
print("  -> 02_filtered_mutation_burden.png")


# =============================================================================
# 3. Filtered Variant Types
# =============================================================================
print("Generating filtered variant types...")
type_counts = combined.groupby(['short', 'var_type']).size().unstack(fill_value=0).reindex(labels)
type_order = ['SNV', 'Insertion', 'Deletion', 'MNV']
type_counts = type_counts.reindex(columns=[c for c in type_order if c in type_counts.columns])
type_colors = {'SNV': '#3498db', 'Insertion': '#2ecc71', 'Deletion': '#e74c3c', 'MNV': '#f39c12'}

fig, ax = plt.subplots(figsize=(12, 5))
type_counts.plot(kind='bar', stacked=True, ax=ax,
                 color=[type_colors.get(c, '#999') for c in type_counts.columns],
                 edgecolor='white', linewidth=0.5)
ax.set_title('Filtered Somatic Variant Types per Sample', fontsize=14, fontweight='bold')
ax.set_ylabel('Variant Count')
ax.set_xlabel('')
ax.set_xticklabels(type_counts.index, rotation=45, ha='right')
ax.legend(title='Variant Type', bbox_to_anchor=(1.02, 1), loc='upper left')
ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: f'{int(y):,}'))
plt.tight_layout()
plt.savefig(f'{OUT_DIR}/03_filtered_variant_types.png')
plt.close()
print("  -> 03_filtered_variant_types.png")


# =============================================================================
# 4. Filtered Mutation Spectrum
# =============================================================================
print("Generating filtered mutation spectrum...")
snvs = combined[combined['var_type'] == 'SNV'].copy()
snvs['sub_class'] = snvs.apply(lambda r: get_substitution_class(r['ref'], r['alt']), axis=1)

sub_order = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
spec_counts = snvs.groupby(['short', 'sub_class']).size().unstack(fill_value=0).reindex(labels)
spec_counts = spec_counts.reindex(columns=sub_order)
spec_pct = spec_counts.div(spec_counts.sum(axis=1), axis=0) * 100

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 9))
fig.suptitle('Filtered Somatic SNV Mutation Spectrum', fontsize=15, fontweight='bold')

spec_counts.plot(kind='bar', stacked=True, ax=ax1,
                 color=[SUB_COLORS[s] for s in sub_order], edgecolor='white', linewidth=0.5)
ax1.set_ylabel('SNV Count')
ax1.set_title('Absolute Counts')
ax1.set_xlabel('')
ax1.legend(title='Substitution', bbox_to_anchor=(1.02, 1), loc='upper left')
ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: f'{int(y):,}'))

spec_pct.plot(kind='bar', stacked=True, ax=ax2,
              color=[SUB_COLORS[s] for s in sub_order], edgecolor='white', linewidth=0.5)
ax2.set_ylabel('Proportion (%)')
ax2.set_title('Relative Proportions')
ax2.set_xlabel('')
ax2.set_xticklabels(spec_pct.index, rotation=45, ha='right')
ax2.legend(title='Substitution', bbox_to_anchor=(1.02, 1), loc='upper left')

plt.tight_layout()
plt.savefig(f'{OUT_DIR}/04_filtered_mutation_spectrum.png')
plt.close()
print("  -> 04_filtered_mutation_spectrum.png")


# =============================================================================
# 5. VEP Impact
# =============================================================================
print("Generating filtered VEP impact...")
impact_counts = combined.groupby(['short', 'impact']).size().unstack(fill_value=0).reindex(labels)
impact_order = ['HIGH', 'MODERATE', 'LOW', 'MODIFIER']
impact_counts = impact_counts.reindex(columns=[c for c in impact_order if c in impact_counts.columns])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5.5))
fig.suptitle('Filtered VEP Variant Impact Classification', fontsize=15, fontweight='bold')

impact_counts.plot(kind='bar', stacked=True, ax=ax1,
                   color=[IMPACT_COLORS[c] for c in impact_counts.columns],
                   edgecolor='white', linewidth=0.5)
ax1.set_ylabel('Variant Count')
ax1.set_title('Absolute Counts by Impact')
ax1.set_xlabel('')
ax1.legend(title='Impact', bbox_to_anchor=(1.02, 1), loc='upper left')
ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: f'{int(y):,}'))

impact_pct = impact_counts.div(impact_counts.sum(axis=1), axis=0) * 100
impact_pct.plot(kind='bar', stacked=True, ax=ax2,
                color=[IMPACT_COLORS[c] for c in impact_counts.columns],
                edgecolor='white', linewidth=0.5)
ax2.set_ylabel('Proportion (%)')
ax2.set_title('Relative Proportions by Impact')
ax2.set_xlabel('')
ax2.set_xticklabels(impact_pct.index, rotation=45, ha='right')
ax2.legend(title='Impact', bbox_to_anchor=(1.02, 1), loc='upper left')

plt.tight_layout()
plt.savefig(f'{OUT_DIR}/05_filtered_vep_impact.png')
plt.close()
print("  -> 05_filtered_vep_impact.png")

# Top consequences
all_conseq = combined['consequence'].str.split('&').explode()
top_conseq = all_conseq.value_counts().head(15)

fig, ax = plt.subplots(figsize=(10, 6))
top_conseq.plot(kind='barh', ax=ax, color='#3498db', edgecolor='white')
ax.set_xlabel('Total Variant Count')
ax.set_title('Top 15 VEP Consequence Types — Filtered Variants', fontsize=14, fontweight='bold')
ax.invert_yaxis()
ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: f'{int(y):,}'))
plt.tight_layout()
plt.savefig(f'{OUT_DIR}/06_filtered_top_consequences.png')
plt.close()
print("  -> 06_filtered_top_consequences.png")


# =============================================================================
# 6. Recurrently Mutated Genes (filtered)
# =============================================================================
print("Generating filtered recurrent genes...")
coding = combined[(combined['impact'].isin(['HIGH', 'MODERATE'])) & (combined['gene'] != '')].copy()
gene_sample_counts = coding.groupby('gene')['sample'].nunique().sort_values(ascending=False)
gene_mut_counts = coding.groupby('gene').size()

top_genes = gene_sample_counts.head(30)
top_gene_muts = gene_mut_counts.reindex(top_genes.index)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7))
fig.suptitle('Top Recurrently Mutated Genes — Filtered (HIGH/MODERATE Impact)', fontsize=15, fontweight='bold')

ax1.barh(range(len(top_genes)), top_genes.values, color='#e74c3c', edgecolor='white')
ax1.set_yticks(range(len(top_genes)))
ax1.set_yticklabels(top_genes.index, fontsize=9)
ax1.set_xlabel('Number of Samples Mutated')
ax1.set_title('Sample Recurrence')
ax1.invert_yaxis()

ax2.barh(range(len(top_gene_muts)), top_gene_muts.values, color='#3498db', edgecolor='white')
ax2.set_yticks(range(len(top_gene_muts)))
ax2.set_yticklabels(top_gene_muts.index, fontsize=9)
ax2.set_xlabel('Total Mutation Count')
ax2.set_title('Total Mutations')
ax2.invert_yaxis()

plt.tight_layout()
plt.savefig(f'{OUT_DIR}/07_filtered_recurrent_genes.png')
plt.close()
print("  -> 07_filtered_recurrent_genes.png")


# =============================================================================
# 7. VAF Distribution (filtered)
# =============================================================================
print("Generating filtered VAF distribution...")
fig, ax = plt.subplots(figsize=(12, 5))
for sample in SAMPLES:
    short = SAMPLE_SHORT[sample]
    af_vals = all_variants[sample]['af']
    ax.hist(af_vals, bins=50, range=(0, 1), alpha=0.4, label=short, histtype='step', linewidth=1.5)

ax.set_xlabel('Variant Allele Frequency (VAF)')
ax.set_ylabel('Count')
ax.set_title('Filtered Somatic VAF Distribution per Sample', fontsize=14, fontweight='bold')
ax.legend(ncol=4, fontsize=8, bbox_to_anchor=(1.02, 1), loc='upper left')
plt.tight_layout()
plt.savefig(f'{OUT_DIR}/08_filtered_vaf_distribution.png')
plt.close()
print("  -> 08_filtered_vaf_distribution.png")


# =============================================================================
# 8. Chromosomal Density (filtered)
# =============================================================================
print("Generating filtered chromosomal density...")
main_chroms = [str(i) for i in range(1, 20)] + ['X', 'Y']
chrom_var_counts = combined[combined['chrom'].isin(main_chroms)].groupby(['short', 'chrom']).size().unstack(fill_value=0)
chrom_var_counts = chrom_var_counts.reindex(columns=main_chroms, fill_value=0).reindex(labels)

fig, ax = plt.subplots(figsize=(14, 7))
im = ax.imshow(chrom_var_counts.values, aspect='auto', cmap='YlOrRd')
ax.set_xticks(range(len(main_chroms)))
ax.set_xticklabels(main_chroms)
ax.set_yticks(range(len(chrom_var_counts)))
ax.set_yticklabels(chrom_var_counts.index)
ax.set_xlabel('Chromosome')
ax.set_ylabel('Sample')
ax.set_title('Filtered Somatic Variant Density by Chromosome', fontsize=14, fontweight='bold')
cbar = plt.colorbar(im, ax=ax, shrink=0.7)
cbar.set_label('Variant Count')
plt.tight_layout()
plt.savefig(f'{OUT_DIR}/09_filtered_chrom_density.png')
plt.close()
print("  -> 09_filtered_chrom_density.png")


# =============================================================================
# 9. High impact variants table
# =============================================================================
print("Generating filtered high impact variant table...")
high_impact = combined[combined['impact'] == 'HIGH'].copy()
high_summary = high_impact.groupby(['gene', 'consequence']).agg(
    n_samples=('sample', 'nunique'),
    n_variants=('pos', 'count'),
    samples=('short', lambda x: ', '.join(sorted(set(x))))
).reset_index().sort_values(['n_samples', 'n_variants'], ascending=[False, False])
high_summary.to_csv(f'{OUT_DIR}/filtered_high_impact_variants.csv', index=False)
print(f"  -> filtered_high_impact_variants.csv ({len(high_summary)} entries)")


# =============================================================================
# 10. Summary table
# =============================================================================
print("Generating filtered summary statistics...")
summary_rows = []
for sample in SAMPLES:
    df = all_variants[sample]
    snv_count = (df['var_type'] == 'SNV').sum()
    indel_count = df['var_type'].isin(['Insertion', 'Deletion']).sum()
    high_count = (df['impact'] == 'HIGH').sum()
    mod_count = (df['impact'] == 'MODERATE').sum()
    coding_genes = df[(df['impact'].isin(['HIGH', 'MODERATE'])) & (df['gene'] != '')]['gene'].nunique()
    tmb = len(df) / GENOME_SIZE_MB

    summary_rows.append({
        'Sample': SAMPLE_SHORT[sample],
        'Unfiltered': unfilt_counts[sample],
        'Filtered': len(df),
        'SNVs': snv_count,
        'Indels': indel_count,
        'TMB (mut/Mb)': round(tmb, 2),
        'HIGH Impact': high_count,
        'MODERATE Impact': mod_count,
        'Genes Affected': coding_genes,
    })

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(f'{OUT_DIR}/filtered_variant_summary.csv', index=False)
print("  -> filtered_variant_summary.csv")
print(summary_df.to_string(index=False))
print(f"\nAll figures saved to {OUT_DIR}/")
print("Done!")
