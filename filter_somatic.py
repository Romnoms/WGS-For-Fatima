#!/usr/bin/env python3
"""
Somatic Variant Filtering Pipeline
Applies two-pass filtering to remove germline contamination from tumor-only calls:
  1. MGP (Mouse Genomes Project) allele frequency filter — removes known strain variants
  2. Cross-sample recurrence filter — removes variants present in too many samples
"""

import gzip
import os
import re
import glob
from collections import defaultdict

DATA_DIR = os.path.expanduser('~/wgs-project/data/05.SomaticSNV')
OUT_DIR = os.path.expanduser('~/wgs-project/data/filtered_somatic')
os.makedirs(OUT_DIR, exist_ok=True)

SAMPLES = sorted([
    os.path.basename(f).replace('.somatic.final.vep.vcf.gz', '')
    for f in glob.glob(f'{DATA_DIR}/*.somatic.final.vep.vcf.gz')
])

MGP_AF_THRESHOLD = 0.01  # Remove if MGP_AF >= 1%
MAX_SAMPLE_RECURRENCE = 10  # Remove if variant seen in > 10 of 16 samples (>62%)


def extract_mgp_af(info_field):
    """Extract MGP_AF from the CSQ annotation in the INFO field."""
    m = re.search(r'MGP_AF=([0-9.eE+-]+)', info_field)
    if m:
        return float(m.group(1))
    return None


def parse_variant_key(line):
    """Return (chrom, pos, ref, alt) from a VCF data line."""
    parts = line.split('\t', 5)
    return (parts[0], parts[1], parts[3], parts[4])


# ==========================================================================
# Pass 1: Index all variants by position across samples (after MGP filter)
# ==========================================================================
print("Pass 1: Building cross-sample variant index (post-MGP filter)...")
variant_samples = defaultdict(set)  # variant_key -> set of samples
sample_mgp_passed = {}  # sample -> list of (line, variant_key)

for sample in SAMPLES:
    short = sample.replace('26034XD-04-', 'S')
    vcf_path = f'{DATA_DIR}/{sample}.somatic.final.vep.vcf.gz'
    passed = []
    total = 0
    mgp_removed = 0

    with gzip.open(vcf_path, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            total += 1
            parts = line.split('\t', 8)
            info = parts[7]

            mgp_af = extract_mgp_af(info)
            if mgp_af is not None and mgp_af >= MGP_AF_THRESHOLD:
                mgp_removed += 1
                continue

            vkey = parse_variant_key(line)
            variant_samples[vkey].add(sample)
            passed.append((line, vkey))

    sample_mgp_passed[sample] = passed
    print(f"  {short}: {total} total -> {len(passed)} after MGP filter "
          f"({mgp_removed} removed, {mgp_removed/total*100:.1f}%)")


# ==========================================================================
# Pass 2: Apply cross-sample recurrence filter and write filtered VCFs
# ==========================================================================
print(f"\nPass 2: Applying cross-sample recurrence filter (max {MAX_SAMPLE_RECURRENCE} samples)...")

stats = {}
for sample in SAMPLES:
    short = sample.replace('26034XD-04-', 'S')
    vcf_in = f'{DATA_DIR}/{sample}.somatic.final.vep.vcf.gz'
    vcf_out = f'{OUT_DIR}/{sample}.filtered.vcf.gz'

    # Read header from original
    header_lines = []
    with gzip.open(vcf_in, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                header_lines.append(line)
            else:
                break

    # Add filter descriptions to header
    filter_lines = [
        '##FILTER=<ID=mgp_germline,Description="Variant has MGP allele frequency >= 0.01, likely germline">\n',
        '##FILTER=<ID=recurrent_germline,Description="Variant present in >10 of 16 samples, likely germline">\n',
        f'##filtering_command="filter_somatic.py MGP_AF>={MGP_AF_THRESHOLD}, recurrence>{MAX_SAMPLE_RECURRENCE}/16"\n',
    ]
    # Insert before #CHROM line
    header_out = header_lines[:-1] + filter_lines + [header_lines[-1]]

    recurrence_removed = 0
    kept = 0
    passed_lines = sample_mgp_passed[sample]

    with gzip.open(vcf_out, 'wt') as fh:
        for h in header_out:
            fh.write(h)

        for line, vkey in passed_lines:
            n_samples = len(variant_samples[vkey])
            if n_samples > MAX_SAMPLE_RECURRENCE:
                recurrence_removed += 1
                continue
            fh.write(line)
            kept += 1

    total_original = kept + recurrence_removed + (
        sum(1 for _ in gzip.open(vcf_in, 'rt') if not _.startswith('#')) - len(passed_lines) - recurrence_removed
    )

    stats[sample] = {
        'mgp_passed': len(passed_lines),
        'recurrence_removed': recurrence_removed,
        'final': kept,
    }
    print(f"  {short}: {len(passed_lines)} post-MGP -> {kept} final "
          f"({recurrence_removed} recurrent removed)")

# Print summary
print("\n" + "=" * 70)
print("FILTERING SUMMARY")
print("=" * 70)
print(f"{'Sample':<8} {'Original':>10} {'Post-MGP':>10} {'Final':>10} {'% Retained':>12}")
print("-" * 70)

for sample in SAMPLES:
    short = sample.replace('26034XD-04-', 'S')
    vcf_in = f'{DATA_DIR}/{sample}.somatic.final.vep.vcf.gz'
    original = sum(1 for line in gzip.open(vcf_in, 'rt') if not line.startswith('#'))
    s = stats[sample]
    pct = s['final'] / original * 100 if original > 0 else 0
    print(f"{short:<8} {original:>10,} {s['mgp_passed']:>10,} {s['final']:>10,} {pct:>11.1f}%")

print("=" * 70)
print(f"\nFiltered VCFs written to: {OUT_DIR}/")
