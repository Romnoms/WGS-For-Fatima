#!/usr/bin/env python3
"""
Circos Plots — WGS Somatic Variants in Sucrose-Induced Mouse Tumors
Uses pycirclize to generate publication-quality circos plots showing:
  - Outer ring: Mouse chromosome ideogram (GRCm39)
  - Track 1: Variant density heatmap (1 Mb windows)
  - Track 2: SNV substitution spectrum (stacked area by 6-class)
  - Track 3: Indel density
  - Track 4: HIGH impact variant markers
  - Center: Sample label + TMB

Generates:
  1. Per-sample circos (2x2 grid of representative samples)
  2. Cohort-wide circos (all samples aggregated)
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
from matplotlib.patches import Patch
import warnings
warnings.filterwarnings('ignore')

from pycirclize import Circos

plt.rcParams.update({
    'font.size': 10,
    'font.family': 'sans-serif',
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
GENOME_SIZE_MB = 2728.22
COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

SUB_COLORS = {
    'C>A': '#1EBBD7', 'C>G': '#050708', 'C>T': '#E62725',
    'T>A': '#CBCACB', 'T>C': '#A1CF64', 'T>G': '#EDC8C5',
}
SUB_ORDER = ['C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G']

# GRCm39 chromosome sizes
CHROM_SIZES = {
    '1': 195154279, '2': 181755017, '3': 159745316, '4': 156860686,
    '5': 151758149, '6': 149588044, '7': 144995196, '8': 130127694,
    '9': 124359700, '10': 130530862, '11': 121973369, '12': 120092757,
    '13': 120883175, '14': 125139656, '15': 104073951, '16': 98008968,
    '17': 95294699, '18': 90720763, '19': 61420004, 'X': 169476592,
    'Y': 91455967,
}
MAIN_CHROMS = [str(i) for i in range(1, 20)] + ['X', 'Y']
WINDOW_SIZE = 2_000_000  # 2 Mb windows


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
                'af': af, 'var_type': var_type,
                'gene': gene, 'consequence': consequence, 'impact': impact
            })
    return pd.DataFrame(variants)


def get_sub_class(ref, alt):
    if ref in ('C', 'T'):
        return f'{ref}>{alt}'
    return f'{COMPLEMENT[ref]}>{COMPLEMENT[alt]}'


def compute_windowed_density(df, chrom, chrom_size, window=WINDOW_SIZE):
    """Compute variant counts per window along a chromosome."""
    n_windows = int(np.ceil(chrom_size / window))
    density = np.zeros(n_windows)
    chrom_df = df[df['chrom'] == chrom]
    for pos in chrom_df['pos']:
        idx = min(int(pos / window), n_windows - 1)
        density[idx] += 1
    return density


def compute_windowed_spectrum(df, chrom, chrom_size, window=WINDOW_SIZE):
    """Compute per-window substitution class counts."""
    n_windows = int(np.ceil(chrom_size / window))
    spectrum = {sub: np.zeros(n_windows) for sub in SUB_ORDER}
    snvs = df[(df['chrom'] == chrom) & (df['var_type'] == 'SNV')]
    for _, row in snvs.iterrows():
        idx = min(int(row['pos'] / window), n_windows - 1)
        sub = get_sub_class(row['ref'], row['alt'])
        if sub in spectrum:
            spectrum[sub][idx] += 1
    return spectrum


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


def make_circos(df, title, filename, show_genes=True):
    """Generate a single circos plot for a given variant DataFrame."""
    # Sector data: chromosome name -> size
    sectors = {f'chr{c}': CHROM_SIZES[c] for c in MAIN_CHROMS}

    circos = Circos(sectors, space=3)

    for sector in circos.sectors:
        chrom_name = sector.name  # e.g. 'chr1'
        chrom = chrom_name.replace('chr', '')
        chrom_size = CHROM_SIZES[chrom]

        # Chromosome label
        sector.text(sector.name, r=115, fontsize=7, fontweight='bold')

        # Axis ticks on outer ring
        sector.axis(fc='#e8e8e8', ec='#333333', lw=0.8)

        # Track 1: Total variant density heatmap
        density = compute_windowed_density(df, chrom, chrom_size)
        n_win = len(density)
        positions = np.linspace(0, chrom_size, n_win + 1)

        track1 = sector.add_track((88, 98))
        track1.axis(fc='none', ec='gray', lw=0.3)

        # Normalize density for color mapping
        max_d = max(density.max(), 1)
        for i in range(n_win):
            start, end = positions[i], positions[i + 1]
            intensity = min(density[i] / max_d, 1.0)
            color = plt.cm.YlOrRd(intensity)
            track1.rect(start, end, fc=color, ec='none', lw=0)

        # Track 2: SNV substitution spectrum (stacked bars)
        spectrum = compute_windowed_spectrum(df, chrom, chrom_size)
        track2 = sector.add_track((73, 87))
        track2.axis(fc='none', ec='gray', lw=0.3)

        # Stack the substitution types
        total_per_window = np.zeros(n_win)
        for sub in SUB_ORDER:
            total_per_window += spectrum[sub]
        max_total = max(total_per_window.max(), 1)

        for i in range(n_win):
            start, end = positions[i], positions[i + 1]
            bottom = 0
            for sub in SUB_ORDER:
                val = spectrum[sub][i]
                if val > 0:
                    height_frac = val / max_total
                    bar_bottom = 73 + bottom * (87 - 73)
                    bar_top = 73 + (bottom + height_frac) * (87 - 73)
                    # Use rect with radial positioning
                    track2.rect(start, end, fc=SUB_COLORS[sub], ec='none', lw=0)
                    bottom += height_frac

        # For track 2, use line plot instead (cleaner)
        # Plot dominant substitution as line
        ct_density = spectrum.get('C>T', np.zeros(n_win))
        x_mid = np.array([(positions[i] + positions[i+1]) / 2 for i in range(n_win)])
        if ct_density.max() > 0:
            ct_norm = ct_density / max(ct_density.max(), 1) * (87 - 73)
            track2.fill_between(x_mid, ct_norm + 73, 73, fc='#E62725', alpha=0.5, ec='none')

        # Overlay T>C
        tc_density = spectrum.get('T>C', np.zeros(n_win))
        if tc_density.max() > 0:
            tc_norm = tc_density / max(ct_density.max(), 1) * (87 - 73)
            track2.fill_between(x_mid, tc_norm + 73, 73, fc='#A1CF64', alpha=0.4, ec='none')

        # Track 3: Indel density
        indels = df[(df['chrom'] == chrom) & (df['var_type'].isin(['Insertion', 'Deletion']))]
        indel_density = np.zeros(n_win)
        for pos in indels['pos']:
            idx = min(int(pos / WINDOW_SIZE), n_win - 1)
            indel_density[idx] += 1

        track3 = sector.add_track((60, 72))
        track3.axis(fc='none', ec='gray', lw=0.3)
        if indel_density.max() > 0:
            indel_norm = indel_density / max(indel_density.max(), 1) * (72 - 60)
            track3.fill_between(x_mid, indel_norm + 60, 60, fc='#3498db', alpha=0.6, ec='none')

        # Track 4: HIGH impact markers
        if show_genes:
            high = df[(df['chrom'] == chrom) & (df['impact'] == 'HIGH')]
            track4 = sector.add_track((50, 59))
            track4.axis(fc='none', ec='gray', lw=0.3)
            for _, row in high.iterrows():
                track4.scatter([row['pos']], [55], s=15, c='#d62728',
                              zorder=5, alpha=0.8)

    fig = circos.plotfig()

    # Add title
    fig.text(0.5, 0.97, title, ha='center', va='top',
             fontsize=14, fontweight='bold')

    # Add legend
    legend_elements = [
        Patch(fc='#E62725', alpha=0.5, label='C>T density'),
        Patch(fc='#A1CF64', alpha=0.4, label='T>C density'),
        Patch(fc='#3498db', alpha=0.6, label='Indel density'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#d62728',
                   markersize=6, label='HIGH impact'),
        Patch(fc=plt.cm.YlOrRd(0.2), label='Low variant density'),
        Patch(fc=plt.cm.YlOrRd(0.8), label='High variant density'),
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=3,
               fontsize=8, bbox_to_anchor=(0.5, 0.01), frameon=True,
               fancybox=True, shadow=False)

    plt.savefig(f'{OUT_DIR}/{filename}')
    plt.close()
    print(f"  -> {filename}")
    return fig


# =============================================================================
# 1. Cohort-wide circos (all samples combined)
# =============================================================================
print("\nGenerating cohort-wide circos plot...")
n_total = len(combined)
tmb_cohort = n_total / GENOME_SIZE_MB
make_circos(combined,
            f'Cohort-Wide Somatic Variant Landscape (n={n_total:,}, TMB={tmb_cohort:.1f} mut/Mb)',
            '19_circos_cohort.png')


# =============================================================================
# 2. Per-sample circos (2x2 grid)
# =============================================================================
print("\nGenerating per-sample circos plots (2x2 grid)...")

select_shorts = ['S15', 'S04', 'S16', 'S01']
select_samples = []
for short in select_shorts:
    for s in SAMPLES:
        if SAMPLE_SHORT[s] == short:
            select_samples.append(s)
            break

fig_grid, axes_grid = plt.subplots(2, 2, figsize=(20, 20),
                                    subplot_kw={'projection': 'polar'})
plt.close(fig_grid)  # We'll build individual plots instead

# Generate individual circos for each sample
for sample in select_samples:
    short = SAMPLE_SHORT[sample]
    df = all_variants[sample]
    tmb = len(df) / GENOME_SIZE_MB
    n_high = (df['impact'] == 'HIGH').sum()
    indel_pct = df['var_type'].isin(['Insertion', 'Deletion']).sum() / len(df) * 100

    make_circos(
        df,
        f'{short} — Somatic Variant Circos (n={len(df):,}, TMB={tmb:.1f}, '
        f'Indel={indel_pct:.0f}%, HIGH={n_high})',
        f'20_circos_{short}.png'
    )


# =============================================================================
# 3. Comparative multi-sample circos using concentric rings
# =============================================================================
print("\nGenerating comparative multi-sample circos...")

sectors = {f'chr{c}': CHROM_SIZES[c] for c in MAIN_CHROMS}
circos = Circos(sectors, space=3)

# Color each sample differently
sample_colors = {
    'S15': '#e74c3c',
    'S04': '#3498db',
    'S16': '#2ecc71',
    'S01': '#9b59b6',
}

for sector in circos.sectors:
    chrom = sector.name.replace('chr', '')
    chrom_size = CHROM_SIZES[chrom]
    n_win = int(np.ceil(chrom_size / WINDOW_SIZE))
    positions = np.linspace(0, chrom_size, n_win + 1)
    x_mid = np.array([(positions[i] + positions[i+1]) / 2 for i in range(n_win)])

    sector.text(sector.name, r=115, fontsize=7, fontweight='bold')
    sector.axis(fc='#f0f0f0', ec='#333333', lw=0.8)

    # One ring per sample
    ring_ranges = [(85, 97), (72, 84), (59, 71), (46, 58)]

    for idx, (sample, (r_inner, r_outer)) in enumerate(
            zip(select_samples, ring_ranges)):
        short = SAMPLE_SHORT[sample]
        df = all_variants[sample]
        color = sample_colors.get(short, '#999')

        density = compute_windowed_density(df, chrom, chrom_size)
        track = sector.add_track((r_inner, r_outer))
        track.axis(fc='none', ec='gray', lw=0.3)

        if density.max() > 0:
            density_norm = density / max(density.max(), 1) * (r_outer - r_inner)
            track.fill_between(x_mid, density_norm + r_inner, r_inner,
                              fc=color, alpha=0.5, ec='none')

fig = circos.plotfig()

fig.text(0.5, 0.97,
         'Comparative Variant Density — 4 Representative Samples',
         ha='center', va='top', fontsize=14, fontweight='bold')
fig.text(0.5, 0.94,
         'Concentric rings (outer to inner): S15, S04, S16, S01',
         ha='center', va='top', fontsize=10, color='gray')

legend_elements = [
    Patch(fc=sample_colors['S15'], alpha=0.5, label=f'S15 (hypermutated)'),
    Patch(fc=sample_colors['S04'], alpha=0.5, label=f'S04 (hypermutated)'),
    Patch(fc=sample_colors['S16'], alpha=0.5, label=f'S16 (elevated)'),
    Patch(fc=sample_colors['S01'], alpha=0.5, label=f'S01 (typical)'),
]
fig.legend(handles=legend_elements, loc='lower center', ncol=4,
           fontsize=9, bbox_to_anchor=(0.5, 0.01), frameon=True)

plt.savefig(f'{OUT_DIR}/21_circos_comparative.png')
plt.close()
print("  -> 21_circos_comparative.png")


print("\nAll circos plots complete!")
