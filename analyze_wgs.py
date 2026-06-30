#!/usr/bin/env python3
"""
WGS Analysis Pipeline — Mouse Tumor Samples (GRCm39)
Generates publication-quality figures from DRAGEN somatic variant calls.
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

DATA_DIR = os.path.expanduser('~/wgs-project/data')
OUT_DIR = os.path.expanduser('~/wgs-project/figures')
os.makedirs(OUT_DIR, exist_ok=True)

SAMPLES = sorted([
    os.path.basename(f).replace('.wgs_overall_mean_cov.csv', '')
    for f in glob.glob(f'{DATA_DIR}/01.QualityCheck/DragenStats/*.wgs_overall_mean_cov.csv')
])
SAMPLE_SHORT = {s: s.replace('26034XD-04-', 'S') for s in SAMPLES}

# --- Color palettes ---
IMPACT_COLORS = {'HIGH': '#d62728', 'MODERATE': '#ff7f0e', 'LOW': '#2ca02c', 'MODIFIER': '#aec7e8'}
SUB_COLORS = {'C>A': '#3498db', 'C>G': '#2c3e50', 'C>T': '#e74c3c',
              'T>A': '#bdc3c7', 'T>C': '#27ae60', 'T>G': '#f39c12'}
COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def parse_mapping_metrics(sample):
    """Parse Dragen mapping metrics CSV for a given sample."""
    f = f'{DATA_DIR}/01.QualityCheck/DragenStats/{sample}.mapping_metrics.csv'
    metrics = {}
    with open(f) as fh:
        for line in fh:
            if not line.startswith('MAPPING/ALIGNING SUMMARY'):
                continue
            parts = line.strip().split(',')
            key = parts[2]
            val = parts[3] if len(parts) > 3 else None
            pct = parts[4] if len(parts) > 4 else None
            metrics[key] = (val, pct)
    return metrics


def parse_coverage_metrics(sample):
    """Parse Dragen coverage metrics CSV."""
    f = f'{DATA_DIR}/01.QualityCheck/DragenStats/{sample}.wgs_coverage_metrics.csv'
    metrics = {}
    with open(f) as fh:
        for line in fh:
            parts = line.strip().split(',')
            key = parts[2]
            val = parts[3] if len(parts) > 3 else None
            metrics[key] = val
    return metrics


def parse_ploidy(sample):
    """Parse ploidy estimation."""
    f = f'{DATA_DIR}/01.QualityCheck/DragenStats/{sample}.ploidy_estimation_metrics.csv'
    with open(f) as fh:
        for line in fh:
            if 'Ploidy estimation' in line:
                return line.strip().split(',')[-1]
    return 'NA'


def parse_somatic_vcf(vcf_path):
    """Parse a somatic VCF, extracting variants with VEP annotations."""
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
            chrom, pos, _, ref, alt, _, filt, info = parts[0], int(parts[1]), parts[2], parts[3], parts[4], parts[5], parts[6], parts[7]

            # Parse FORMAT fields
            fmt_keys = parts[8].split(':')
            fmt_vals = parts[9].split(':')
            fmt_dict = dict(zip(fmt_keys, fmt_vals))

            af = float(fmt_dict.get('AF', '0'))
            dp = int(fmt_dict.get('DP', '0'))

            # Determine variant type
            if len(ref) == 1 and len(alt) == 1:
                var_type = 'SNV'
            elif len(ref) < len(alt):
                var_type = 'Insertion'
            elif len(ref) > len(alt):
                var_type = 'Deletion'
            else:
                var_type = 'MNV'

            # Parse CSQ (VEP annotation) - take the first/canonical
            csq_str = ''
            for item in info.split(';'):
                if item.startswith('CSQ='):
                    csq_str = item[4:]
                    break

            gene, consequence, impact = '', '', ''
            if csq_str and csq_fields:
                # Take the most severe consequence (first CSQ entry)
                first_csq = csq_str.split(',')[0]
                csq_vals = first_csq.split('|')
                csq_dict = dict(zip(csq_fields, csq_vals))
                gene = csq_dict.get('SYMBOL', '')
                consequence = csq_dict.get('Consequence', '')
                impact = csq_dict.get('IMPACT', '')

            variants.append({
                'chrom': chrom, 'pos': pos, 'ref': ref, 'alt': alt,
                'filter': filt, 'af': af, 'dp': dp,
                'var_type': var_type, 'gene': gene,
                'consequence': consequence, 'impact': impact
            })

    return pd.DataFrame(variants)


def get_substitution_class(ref, alt):
    """Normalize SNV to pyrimidine context (C>x or T>x)."""
    if ref in ('C', 'T'):
        return f'{ref}>{alt}'
    else:
        return f'{COMPLEMENT[ref]}>{COMPLEMENT[alt]}'


# =============================================================================
# 1. QC Summary Table + Figure
# =============================================================================
print("Generating QC summary...")
qc_data = []
for sample in SAMPLES:
    mm = parse_mapping_metrics(sample)
    cm = parse_coverage_metrics(sample)
    ploidy = parse_ploidy(sample)

    total_reads = int(mm.get('Total input reads', ('0',))[0])
    mapped_pct = float(mm.get('Mapped reads', ('0', '0'))[1])
    dup_pct = float(mm.get('Number of duplicate marked reads', ('0', '0'))[1])
    q30_pct = float(mm.get('Q30 bases', ('0', '0'))[1])
    mean_cov = float(cm.get('Average alignment coverage over genome', '0'))
    pct_20x = float(cm.get('PCT of genome with coverage [  20x: inf)', '0'))

    qc_data.append({
        'Sample': sample,
        'Short': SAMPLE_SHORT[sample],
        'Total Reads (M)': total_reads / 1e6,
        'Mapped %': mapped_pct,
        'Dup %': dup_pct,
        'Q30 %': q30_pct,
        'Mean Cov': mean_cov,
        'PCT >= 20x': pct_20x,
        'Sex': ploidy,
    })

qc_df = pd.DataFrame(qc_data)
qc_df.to_csv(f'{OUT_DIR}/qc_summary.csv', index=False)

# QC multi-panel figure
fig, axes = plt.subplots(2, 2, figsize=(14, 9))
fig.suptitle('WGS Quality Control Summary — Mouse Tumor Samples', fontsize=15, fontweight='bold')

x = range(len(qc_df))
labels = qc_df['Short'].values

# Coverage bar chart
ax = axes[0, 0]
colors = ['#e74c3c' if c < 20 else '#2ecc71' for c in qc_df['Mean Cov']]
ax.bar(x, qc_df['Mean Cov'], color=colors, edgecolor='white', linewidth=0.5)
ax.axhline(y=30, color='gray', linestyle='--', alpha=0.7, label='30x target')
ax.set_ylabel('Mean Coverage (x)')
ax.set_title('Genome Coverage')
ax.set_xticks(x)
ax.set_xticklabels(labels, rotation=45, ha='right')
ax.legend(fontsize=9)

# Mapping rate
ax = axes[0, 1]
ax.bar(x, qc_df['Mapped %'], color='#3498db', edgecolor='white', linewidth=0.5)
ax.set_ylim(97, 100)
ax.set_ylabel('Mapped Reads (%)')
ax.set_title('Mapping Rate')
ax.set_xticks(x)
ax.set_xticklabels(labels, rotation=45, ha='right')

# Duplication rate
ax = axes[1, 0]
colors = ['#e74c3c' if d > 30 else '#f39c12' if d > 20 else '#2ecc71' for d in qc_df['Dup %']]
ax.bar(x, qc_df['Dup %'], color=colors, edgecolor='white', linewidth=0.5)
ax.axhline(y=20, color='gray', linestyle='--', alpha=0.7, label='20% threshold')
ax.set_ylabel('Duplicate Rate (%)')
ax.set_title('Duplication Rate')
ax.set_xticks(x)
ax.set_xticklabels(labels, rotation=45, ha='right')
ax.legend(fontsize=9)

# Q30 bases
ax = axes[1, 1]
ax.bar(x, qc_df['Q30 %'], color='#9b59b6', edgecolor='white', linewidth=0.5)
ax.set_ylim(93, 97)
ax.set_ylabel('Q30 Bases (%)')
ax.set_title('Base Quality (Q30)')
ax.set_xticks(x)
ax.set_xticklabels(labels, rotation=45, ha='right')

plt.tight_layout()
plt.savefig(f'{OUT_DIR}/01_qc_summary.png')
plt.close()
print("  -> 01_qc_summary.png")


# =============================================================================
# 2-6. Parse all somatic VCFs
# =============================================================================
print("Parsing somatic VCFs (this may take a minute)...")
all_variants = {}
for sample in SAMPLES:
    vcf_path = f'{DATA_DIR}/05.SomaticSNV/{sample}.somatic.final.vep.vcf.gz'
    if os.path.exists(vcf_path):
        df = parse_somatic_vcf(vcf_path)
        df['sample'] = sample
        df['short'] = SAMPLE_SHORT[sample]
        all_variants[sample] = df
        print(f"  {SAMPLE_SHORT[sample]}: {len(df)} variants")

combined = pd.concat(all_variants.values(), ignore_index=True)


# =============================================================================
# 2. Mutation Burden
# =============================================================================
print("Generating mutation burden figure...")
# Mouse genome size ~2.73 Gb
GENOME_SIZE_MB = 2728.22

burden = combined.groupby('short').size().reindex([SAMPLE_SHORT[s] for s in SAMPLES])
burden_per_mb = burden / GENOME_SIZE_MB

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
fig.suptitle('Somatic Mutation Burden', fontsize=15, fontweight='bold')

# Total count
colors = ['#e74c3c' if c > 100000 else '#3498db' for c in burden.values]
ax1.bar(range(len(burden)), burden.values, color=colors, edgecolor='white')
ax1.set_xticks(range(len(burden)))
ax1.set_xticklabels(burden.index, rotation=45, ha='right')
ax1.set_ylabel('Total Somatic Variants')
ax1.set_title('Total Variant Count')
ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: f'{int(y):,}'))

# Per Mb (TMB)
colors2 = ['#e74c3c' if c > 50 else '#3498db' for c in burden_per_mb.values]
ax2.bar(range(len(burden_per_mb)), burden_per_mb.values, color=colors2, edgecolor='white')
ax2.set_xticks(range(len(burden_per_mb)))
ax2.set_xticklabels(burden_per_mb.index, rotation=45, ha='right')
ax2.set_ylabel('Mutations / Mb')
ax2.set_title('Tumor Mutational Burden (TMB)')

plt.tight_layout()
plt.savefig(f'{OUT_DIR}/02_mutation_burden.png')
plt.close()
print("  -> 02_mutation_burden.png")


# =============================================================================
# 3. Variant Type Breakdown (SNV, Indel, MNV)
# =============================================================================
print("Generating variant type breakdown...")
type_counts = combined.groupby(['short', 'var_type']).size().unstack(fill_value=0)
type_counts = type_counts.reindex([SAMPLE_SHORT[s] for s in SAMPLES])
type_order = ['SNV', 'Insertion', 'Deletion', 'MNV']
type_counts = type_counts.reindex(columns=[c for c in type_order if c in type_counts.columns])

type_colors = {'SNV': '#3498db', 'Insertion': '#2ecc71', 'Deletion': '#e74c3c', 'MNV': '#f39c12'}

fig, ax = plt.subplots(figsize=(12, 5))
type_counts.plot(kind='bar', stacked=True, ax=ax,
                 color=[type_colors.get(c, '#999') for c in type_counts.columns],
                 edgecolor='white', linewidth=0.5)
ax.set_title('Somatic Variant Types per Sample', fontsize=14, fontweight='bold')
ax.set_ylabel('Variant Count')
ax.set_xlabel('')
ax.set_xticklabels(type_counts.index, rotation=45, ha='right')
ax.legend(title='Variant Type', bbox_to_anchor=(1.02, 1), loc='upper left')
ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: f'{int(y):,}'))
plt.tight_layout()
plt.savefig(f'{OUT_DIR}/03_variant_types.png')
plt.close()
print("  -> 03_variant_types.png")


# =============================================================================
# 4. Mutation Spectrum (6-class substitution types)
# =============================================================================
print("Generating mutation spectrum...")
snvs = combined[combined['var_type'] == 'SNV'].copy()
snvs['sub_class'] = snvs.apply(lambda r: get_substitution_class(r['ref'], r['alt']), axis=1)

sub_order = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']
spec_counts = snvs.groupby(['short', 'sub_class']).size().unstack(fill_value=0)
spec_counts = spec_counts.reindex([SAMPLE_SHORT[s] for s in SAMPLES])
spec_counts = spec_counts.reindex(columns=sub_order)

# Normalize to proportions
spec_pct = spec_counts.div(spec_counts.sum(axis=1), axis=0) * 100

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 9))
fig.suptitle('Somatic SNV Mutation Spectrum', fontsize=15, fontweight='bold')

# Absolute counts stacked
spec_counts.plot(kind='bar', stacked=True, ax=ax1,
                 color=[SUB_COLORS[s] for s in sub_order],
                 edgecolor='white', linewidth=0.5)
ax1.set_ylabel('SNV Count')
ax1.set_title('Absolute Counts')
ax1.set_xlabel('')
ax1.legend(title='Substitution', bbox_to_anchor=(1.02, 1), loc='upper left')
ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: f'{int(y):,}'))

# Proportions stacked
spec_pct.plot(kind='bar', stacked=True, ax=ax2,
              color=[SUB_COLORS[s] for s in sub_order],
              edgecolor='white', linewidth=0.5)
ax2.set_ylabel('Proportion (%)')
ax2.set_title('Relative Proportions')
ax2.set_xlabel('')
ax2.set_xticklabels(spec_pct.index, rotation=45, ha='right')
ax2.legend(title='Substitution', bbox_to_anchor=(1.02, 1), loc='upper left')

plt.tight_layout()
plt.savefig(f'{OUT_DIR}/04_mutation_spectrum.png')
plt.close()
print("  -> 04_mutation_spectrum.png")


# =============================================================================
# 5. VEP Impact & Consequence Summary
# =============================================================================
print("Generating VEP impact/consequence figures...")

# Impact distribution per sample
impact_counts = combined.groupby(['short', 'impact']).size().unstack(fill_value=0)
impact_counts = impact_counts.reindex([SAMPLE_SHORT[s] for s in SAMPLES])
impact_order = ['HIGH', 'MODERATE', 'LOW', 'MODIFIER']
impact_counts = impact_counts.reindex(columns=[c for c in impact_order if c in impact_counts.columns])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5.5))
fig.suptitle('VEP Variant Impact Classification', fontsize=15, fontweight='bold')

impact_counts.plot(kind='bar', stacked=True, ax=ax1,
                   color=[IMPACT_COLORS[c] for c in impact_counts.columns],
                   edgecolor='white', linewidth=0.5)
ax1.set_ylabel('Variant Count')
ax1.set_title('Absolute Counts by Impact')
ax1.set_xlabel('')
ax1.legend(title='Impact', bbox_to_anchor=(1.02, 1), loc='upper left')
ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: f'{int(y):,}'))

# Proportion
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
plt.savefig(f'{OUT_DIR}/05_vep_impact.png')
plt.close()
print("  -> 05_vep_impact.png")

# Top consequences across all samples
all_conseq = combined['consequence'].str.split('&').explode()
top_conseq = all_conseq.value_counts().head(15)

fig, ax = plt.subplots(figsize=(10, 6))
top_conseq.plot(kind='barh', ax=ax, color='#3498db', edgecolor='white')
ax.set_xlabel('Total Variant Count')
ax.set_title('Top 15 VEP Consequence Types (All Samples)', fontsize=14, fontweight='bold')
ax.invert_yaxis()
ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: f'{int(y):,}'))
plt.tight_layout()
plt.savefig(f'{OUT_DIR}/06_top_consequences.png')
plt.close()
print("  -> 06_top_consequences.png")


# =============================================================================
# 6. Top Recurrently Mutated Genes
# =============================================================================
print("Generating recurrently mutated genes...")

# Focus on HIGH + MODERATE impact variants with gene symbols
coding = combined[(combined['impact'].isin(['HIGH', 'MODERATE'])) & (combined['gene'] != '')].copy()

# Count: number of samples each gene is mutated in
gene_sample_counts = coding.groupby('gene')['sample'].nunique().sort_values(ascending=False)
# Count: total mutations per gene
gene_mut_counts = coding.groupby('gene').size()

top_genes = gene_sample_counts.head(30)
top_gene_muts = gene_mut_counts.reindex(top_genes.index)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 7))
fig.suptitle('Top Recurrently Mutated Genes (HIGH/MODERATE Impact)', fontsize=15, fontweight='bold')

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
plt.savefig(f'{OUT_DIR}/07_recurrent_genes.png')
plt.close()
print("  -> 07_recurrent_genes.png")


# =============================================================================
# 7. Allele Frequency Distribution
# =============================================================================
print("Generating allele frequency distributions...")

fig, ax = plt.subplots(figsize=(12, 5))
for sample in SAMPLES:
    short = SAMPLE_SHORT[sample]
    af_vals = all_variants[sample]['af']
    ax.hist(af_vals, bins=50, range=(0, 1), alpha=0.4, label=short, histtype='step', linewidth=1.5)

ax.set_xlabel('Variant Allele Frequency (VAF)')
ax.set_ylabel('Count')
ax.set_title('Somatic VAF Distribution per Sample', fontsize=14, fontweight='bold')
ax.legend(ncol=4, fontsize=8, bbox_to_anchor=(1.02, 1), loc='upper left')
plt.tight_layout()
plt.savefig(f'{OUT_DIR}/08_vaf_distribution.png')
plt.close()
print("  -> 08_vaf_distribution.png")


# =============================================================================
# 8. Genome-wide Variant Density (Chromosome-level)
# =============================================================================
print("Generating genome-wide variant density...")

main_chroms = [str(i) for i in range(1, 20)] + ['X', 'Y']
chrom_var_counts = combined[combined['chrom'].isin(main_chroms)].groupby(['short', 'chrom']).size().unstack(fill_value=0)
chrom_var_counts = chrom_var_counts.reindex(columns=main_chroms)
chrom_var_counts = chrom_var_counts.reindex([SAMPLE_SHORT[s] for s in SAMPLES])

fig, ax = plt.subplots(figsize=(14, 7))
im = ax.imshow(chrom_var_counts.values, aspect='auto', cmap='YlOrRd')
ax.set_xticks(range(len(main_chroms)))
ax.set_xticklabels(main_chroms)
ax.set_yticks(range(len(chrom_var_counts)))
ax.set_yticklabels(chrom_var_counts.index)
ax.set_xlabel('Chromosome')
ax.set_ylabel('Sample')
ax.set_title('Somatic Variant Density by Chromosome', fontsize=14, fontweight='bold')
cbar = plt.colorbar(im, ax=ax, shrink=0.7)
cbar.set_label('Variant Count')
plt.tight_layout()
plt.savefig(f'{OUT_DIR}/09_chrom_density_heatmap.png')
plt.close()
print("  -> 09_chrom_density_heatmap.png")


# =============================================================================
# 9. HIGH impact variant summary table
# =============================================================================
print("Generating HIGH impact variant table...")
high_impact = combined[combined['impact'] == 'HIGH'].copy()
high_summary = high_impact.groupby(['gene', 'consequence']).agg(
    n_samples=('sample', 'nunique'),
    n_variants=('pos', 'count'),
    samples=('short', lambda x: ', '.join(sorted(set(x))))
).reset_index().sort_values(['n_samples', 'n_variants'], ascending=[False, False])

high_summary.to_csv(f'{OUT_DIR}/high_impact_variants.csv', index=False)
print(f"  -> high_impact_variants.csv ({len(high_summary)} entries)")


# =============================================================================
# 10. Summary statistics table
# =============================================================================
print("Generating summary statistics...")
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
        'Total Variants': len(df),
        'SNVs': snv_count,
        'Indels': indel_count,
        'TMB (mut/Mb)': round(tmb, 2),
        'HIGH Impact': high_count,
        'MODERATE Impact': mod_count,
        'Genes Affected': coding_genes,
    })

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(f'{OUT_DIR}/variant_summary.csv', index=False)
print("  -> variant_summary.csv")
print(summary_df.to_string(index=False))

print(f"\nAll figures saved to {OUT_DIR}/")
print("Done!")
