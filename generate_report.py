#!/usr/bin/env python3
"""Generate a self-contained HTML report for WGS analysis."""

import base64
import os
import csv

FIGURES_DIR = os.path.expanduser('~/wgs-project/figures')
OUT_PATH = os.path.expanduser('~/wgs-project/WGS_Analysis_Report.html')


def img_b64(filename):
    with open(os.path.join(FIGURES_DIR, filename), 'rb') as f:
        return base64.b64encode(f.read()).decode()


def read_csv(filename):
    with open(os.path.join(FIGURES_DIR, filename)) as f:
        return list(csv.DictReader(f))


def csv_to_table(rows, columns=None):
    if not rows:
        return ''
    if columns is None:
        columns = list(rows[0].keys())
    html = '<table>\n<thead><tr>'
    for col in columns:
        html += f'<th>{col}</th>'
    html += '</tr></thead>\n<tbody>\n'
    for row in rows:
        html += '<tr>'
        for col in columns:
            val = row.get(col, '')
            html += f'<td>{val}</td>'
        html += '</tr>\n'
    html += '</tbody></table>'
    return html


# Load data
FILT_DIR = os.path.expanduser('~/wgs-project/figures_filtered')


def read_csv_from(directory, filename):
    with open(os.path.join(directory, filename)) as f:
        return list(csv.DictReader(f))


def filt_img_b64(filename):
    with open(os.path.join(FILT_DIR, filename), 'rb') as f:
        return base64.b64encode(f.read()).decode()


qc_data = read_csv('qc_summary.csv')
variant_data = read_csv('variant_summary.csv')
high_impact = read_csv('high_impact_variants.csv')
filtered_variant_data = read_csv_from(FILT_DIR, 'filtered_variant_summary.csv')
filtered_high_impact = read_csv_from(FILT_DIR, 'filtered_high_impact_variants.csv')

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>WGS Analysis Report — Mouse Tumor Samples</title>
<style>
:root {{
    --bg: #ffffff;
    --fg: #1a1a2e;
    --accent: #2563eb;
    --accent-light: #dbeafe;
    --border: #e2e8f0;
    --section-bg: #f8fafc;
    --highlight: #fef3c7;
    --red: #dc2626;
    --green: #16a34a;
}}
* {{ margin: 0; padding: 0; box-sizing: border-box; }}
body {{
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
    color: var(--fg);
    background: var(--bg);
    line-height: 1.6;
    max-width: 1200px;
    margin: 0 auto;
    padding: 2rem;
}}
header {{
    text-align: center;
    padding: 2.5rem 1rem;
    margin-bottom: 2rem;
    background: linear-gradient(135deg, #1e3a5f 0%, #2563eb 100%);
    color: white;
    border-radius: 12px;
}}
header h1 {{ font-size: 2rem; margin-bottom: 0.5rem; font-weight: 700; }}
header p {{ font-size: 1rem; opacity: 0.9; }}
header .meta {{ font-size: 0.85rem; opacity: 0.7; margin-top: 0.5rem; }}
h2 {{
    font-size: 1.5rem;
    color: var(--accent);
    border-bottom: 2px solid var(--accent);
    padding-bottom: 0.4rem;
    margin: 2.5rem 0 1rem;
}}
h3 {{ font-size: 1.15rem; margin: 1.5rem 0 0.5rem; color: #334155; }}
p, li {{ margin-bottom: 0.5rem; }}
ul {{ padding-left: 1.5rem; }}
.figure-container {{
    background: var(--section-bg);
    border: 1px solid var(--border);
    border-radius: 8px;
    padding: 1rem;
    margin: 1.5rem 0;
    text-align: center;
}}
.figure-container img {{
    max-width: 100%;
    height: auto;
    border-radius: 4px;
}}
.figure-caption {{
    font-size: 0.9rem;
    color: #64748b;
    margin-top: 0.5rem;
    font-style: italic;
}}
table {{
    width: 100%;
    border-collapse: collapse;
    margin: 1rem 0;
    font-size: 0.85rem;
}}
th, td {{
    padding: 0.5rem 0.75rem;
    text-align: left;
    border-bottom: 1px solid var(--border);
}}
th {{
    background: var(--accent);
    color: white;
    font-weight: 600;
    position: sticky;
    top: 0;
}}
tbody tr:nth-child(even) {{ background: var(--section-bg); }}
tbody tr:hover {{ background: var(--accent-light); }}
.callout {{
    padding: 1rem 1.25rem;
    border-radius: 8px;
    margin: 1.5rem 0;
    border-left: 4px solid;
}}
.callout-warning {{
    background: #fef3c7;
    border-color: #f59e0b;
}}
.callout-info {{
    background: #dbeafe;
    border-color: #2563eb;
}}
.callout-alert {{
    background: #fee2e2;
    border-color: #dc2626;
}}
.callout strong {{ display: block; margin-bottom: 0.3rem; }}
.stats-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
    gap: 1rem;
    margin: 1.5rem 0;
}}
.stat-card {{
    background: var(--section-bg);
    border: 1px solid var(--border);
    border-radius: 8px;
    padding: 1.25rem;
    text-align: center;
}}
.stat-card .value {{ font-size: 1.8rem; font-weight: 700; color: var(--accent); }}
.stat-card .label {{ font-size: 0.85rem; color: #64748b; margin-top: 0.25rem; }}
.toc {{
    background: var(--section-bg);
    border: 1px solid var(--border);
    border-radius: 8px;
    padding: 1.5rem;
    margin: 1.5rem 0;
}}
.toc h3 {{ margin-top: 0; }}
.toc ol {{ padding-left: 1.5rem; }}
.toc a {{ color: var(--accent); text-decoration: none; }}
.toc a:hover {{ text-decoration: underline; }}
.highlight {{ background: var(--highlight); padding: 0.1rem 0.3rem; border-radius: 3px; }}
footer {{
    text-align: center;
    padding: 2rem;
    margin-top: 3rem;
    border-top: 1px solid var(--border);
    color: #94a3b8;
    font-size: 0.85rem;
}}
@media print {{
    body {{ max-width: 100%; padding: 1rem; }}
    header {{ break-after: page; }}
    .figure-container {{ break-inside: avoid; }}
    table {{ font-size: 0.75rem; }}
}}
</style>
</head>
<body>

<header>
    <h1>Whole Genome Sequencing Analysis Report</h1>
    <p>Mouse Tumor Samples &mdash; DRAGEN Somatic Variant Analysis</p>
    <div class="meta">
        Reference Genome: GRCm39 (Mus musculus) &bull;
        Sequencing Project: 26034-04 &bull;
        Report Date: June 30, 2026
    </div>
</header>

<div class="toc">
    <h3>Table of Contents</h3>
    <ol>
        <li><a href="#overview">Study Overview</a></li>
        <li><a href="#qc">Quality Control</a></li>
        <li><a href="#burden">Mutation Burden</a></li>
        <li><a href="#types">Variant Types</a></li>
        <li><a href="#spectrum">Mutation Spectrum</a></li>
        <li><a href="#impact">Functional Impact</a></li>
        <li><a href="#genes">Recurrently Mutated Genes</a></li>
        <li><a href="#vaf">Allele Frequency Distribution</a></li>
        <li><a href="#chrom">Chromosomal Distribution</a></li>
        <li><a href="#high-impact">High Impact Variants</a></li>
        <li><a href="#caveats">Caveats &amp; Next Steps</a></li>
    </ol>
    <h3 style="margin-top:1rem;">Part II: Filtered Analysis</h3>
    <ol start="12">
        <li><a href="#filtering">Germline Filtering Strategy</a></li>
        <li><a href="#filt-burden">Filtered Mutation Burden</a></li>
        <li><a href="#filt-spectrum">Filtered Mutation Spectrum</a></li>
        <li><a href="#filt-impact">Filtered Functional Impact</a></li>
        <li><a href="#filt-genes">Filtered Recurrently Mutated Genes</a></li>
        <li><a href="#filt-vaf">Filtered VAF Distribution</a></li>
        <li><a href="#filt-chrom">Filtered Chromosomal Distribution</a></li>
        <li><a href="#filt-high">Filtered High Impact Variants</a></li>
    </ol>
</div>

<!-- ================================================================== -->
<h2 id="overview">1. Study Overview</h2>

<p>
    This report summarizes whole genome sequencing (WGS) analysis of <strong>16 mouse tumor samples</strong>
    on a <strong>C57BL/6 (Black 6) background</strong> (all confirmed tumor, no matched normal
    controls available). Samples were sequenced on an Illumina platform (2x151 bp paired-end reads)
    and processed with the <strong>DRAGEN v13.021</strong> pipeline for alignment to the GRCm39
    reference genome and somatic variant calling in tumor-only mode. Variants were annotated with
    Ensembl VEP. Since GRCm39 is derived from the C57BL/6J strain, variants called against this
    reference in B6-background tumors are expected to be predominantly somatic in origin.
</p>

<div class="stats-grid">
    <div class="stat-card">
        <div class="value">16</div>
        <div class="label">Tumor Samples</div>
    </div>
    <div class="stat-card">
        <div class="value">32.3x</div>
        <div class="label">Mean Coverage</div>
    </div>
    <div class="stat-card">
        <div class="value">1.23M</div>
        <div class="label">Total Variants</div>
    </div>
    <div class="stat-card">
        <div class="value">10 XX / 6 XY</div>
        <div class="label">Sex Karyotype</div>
    </div>
</div>

<h3>Sample ID Mapping</h3>
{csv_to_table(qc_data, ['Sample', 'Short', 'Total Reads (M)', 'Mean Cov', 'Sex'])}

<!-- ================================================================== -->
<h2 id="qc">2. Quality Control</h2>

<p>
    Overall sequencing quality is high across all samples. All samples achieve &gt;98.9% mapping rate
    and &gt;94% Q30 base quality. Coverage ranges from 19.4x (S11) to 45.8x (S07), with most samples
    exceeding 30x.
</p>

<div class="callout callout-warning">
    <strong>Elevated Duplication Rates</strong>
    Samples S10 (37.0%) and S11 (38.6%) show duplication rates well above the 20% threshold,
    which may indicate low library complexity or over-sequencing. These samples also have the
    lowest effective coverage (23.2x and 19.4x respectively). S09 and S12 also exceed 20% (24.0%).
</div>

<div class="figure-container">
    <img src="data:image/png;base64,{img_b64('01_qc_summary.png')}" alt="QC Summary">
    <div class="figure-caption">
        Figure 1. Quality control metrics across all 16 samples. (A) Mean genome coverage,
        (B) mapping rate, (C) duplicate read rate, (D) Q30 base quality percentage.
        Red/orange bars indicate samples below recommended thresholds.
    </div>
</div>

<h3>Full QC Summary Table</h3>
{csv_to_table(qc_data, ['Short', 'Total Reads (M)', 'Mapped %', 'Dup %', 'Q30 %', 'Mean Cov', 'PCT >= 20x', 'Sex'])}

<!-- ================================================================== -->
<h2 id="burden">3. Mutation Burden</h2>

<p>
    Most samples carry 39,000&ndash;53,000 somatic variants (TMB 14&ndash;19 mut/Mb), consistent with
    a common strain background calling pattern in tumor-only mode. Two samples are striking outliers:
</p>

<ul>
    <li><span class="highlight">S04: 263,055 variants (96.4 mut/Mb)</span> &mdash; ~6x higher than the cohort median</li>
    <li><span class="highlight">S15: 347,161 variants (127.3 mut/Mb)</span> &mdash; ~8x higher than the cohort median</li>
</ul>

<p>
    This hypermutation phenotype may reflect mismatch repair deficiency, POLE/POLD1 mutations,
    or a fundamentally different tumor biology in these samples. However, germline contamination
    should be ruled out first (see <a href="#caveats">Caveats</a>).
</p>

<div class="figure-container">
    <img src="data:image/png;base64,{img_b64('02_mutation_burden.png')}" alt="Mutation Burden">
    <div class="figure-caption">
        Figure 2. Somatic mutation burden. Left: total variant count per sample. Right: tumor
        mutational burden (TMB) normalized to mutations per megabase. Hypermutated samples
        (S04, S15) are highlighted in red.
    </div>
</div>

<!-- ================================================================== -->
<h2 id="types">4. Variant Types</h2>

<p>
    Across all samples, SNVs comprise the majority of somatic variants (~73&ndash;78%), followed by
    indels (insertions and deletions). The ratio of SNVs to indels is consistent across
    non-hypermutated samples. S16 shows a slightly elevated proportion of indels compared to peers.
</p>

<div class="figure-container">
    <img src="data:image/png;base64,{img_b64('03_variant_types.png')}" alt="Variant Types">
    <div class="figure-caption">
        Figure 3. Variant type composition per sample, showing counts of SNVs, insertions,
        deletions, and MNVs (multi-nucleotide variants).
    </div>
</div>

<h3>Variant Summary Table</h3>
{csv_to_table(variant_data)}

<!-- ================================================================== -->
<h2 id="spectrum">5. Mutation Spectrum</h2>

<p>
    The six-class SNV substitution spectrum is shown below. Non-hypermutated samples share a
    remarkably uniform mutation profile dominated by C&gt;T and T&gt;C transitions, consistent
    with the expected strain-background germline pattern in tumor-only calling.
</p>

<div class="callout callout-info">
    <strong>Notable Spectral Differences in Hypermutated Samples</strong>
    <strong>S04</strong> shows an enrichment in C&gt;T transitions with a broader distribution of
    C&gt;A mutations, potentially indicative of COSMIC SBS1/SBS15 (deamination / MMR deficiency).
    <strong>S15</strong> has a more balanced spectrum with elevated T&gt;C, which may suggest a
    distinct mutational process. Formal mutational signature decomposition (e.g., SigProfiler)
    would be informative.
</div>

<div class="figure-container">
    <img src="data:image/png;base64,{img_b64('04_mutation_spectrum.png')}" alt="Mutation Spectrum">
    <div class="figure-caption">
        Figure 4. Somatic SNV mutation spectrum. Top: absolute counts of each substitution class.
        Bottom: relative proportions normalized per sample. Substitutions are shown in pyrimidine
        context (C&gt;x or T&gt;x).
    </div>
</div>

<!-- ================================================================== -->
<h2 id="impact">6. Functional Impact</h2>

<p>
    As expected for whole-genome somatic calls, the vast majority of variants fall in non-coding
    regions (MODIFIER impact). HIGH and MODERATE impact variants (affecting protein-coding sequences)
    represent a small but biologically important fraction.
</p>

<div class="figure-container">
    <img src="data:image/png;base64,{img_b64('05_vep_impact.png')}" alt="VEP Impact">
    <div class="figure-caption">
        Figure 5. VEP functional impact classification. Left: absolute variant counts by impact
        tier. Right: proportional distribution showing the consistent fraction of coding variants.
    </div>
</div>

<div class="figure-container">
    <img src="data:image/png;base64,{img_b64('06_top_consequences.png')}" alt="Top Consequences">
    <div class="figure-caption">
        Figure 6. Top 15 VEP consequence types across all samples. Intergenic and intronic variants
        dominate, as expected for WGS somatic calls.
    </div>
</div>

<!-- ================================================================== -->
<h2 id="genes">7. Recurrently Mutated Genes</h2>

<p>
    The most recurrently mutated genes (HIGH/MODERATE impact variants present across the most
    samples) are shown below. Several gene families dominate: zinc finger proteins (Zfp),
    vomeronasal receptors (Vmn2r), and Mroh2a.
</p>

<div class="callout callout-alert">
    <strong>Interpretation Caution</strong>
    The top recurrently mutated genes are largely from highly polymorphic gene families in mice
    (Zfp, Vmn2r, Mroh2a). Their presence across all 16 samples strongly suggests these represent
    <strong>germline strain variants</strong> rather than true somatic cancer drivers. This is
    expected in tumor-only somatic calling without a matched normal. Genuine tumor-specific
    recurrent mutations would be expected to appear in a subset of samples.
</div>

<div class="figure-container">
    <img src="data:image/png;base64,{img_b64('07_recurrent_genes.png')}" alt="Recurrent Genes">
    <div class="figure-caption">
        Figure 7. Top 30 recurrently mutated genes with HIGH or MODERATE impact variants.
        Left: number of samples in which each gene is mutated. Right: total mutation count
        across all samples.
    </div>
</div>

<!-- ================================================================== -->
<h2 id="vaf">8. Allele Frequency Distribution</h2>

<p>
    The variant allele frequency (VAF) distribution provides insight into clonal architecture
    and germline contamination. The dominant peak at VAF ~0.5 across all samples is characteristic
    of <strong>heterozygous germline variants</strong>, not subclonal somatic mutations.
</p>

<p>
    Hypermutated samples S04 and S15 show a broader VAF distribution extending to lower frequencies,
    which may contain true subclonal somatic variants mixed with the germline background.
</p>

<div class="figure-container">
    <img src="data:image/png;base64,{img_b64('08_vaf_distribution.png')}" alt="VAF Distribution">
    <div class="figure-caption">
        Figure 8. Variant allele frequency (VAF) distribution for somatic variants per sample.
        The peak at ~0.5 indicates substantial germline variant contamination in the tumor-only
        call set.
    </div>
</div>

<!-- ================================================================== -->
<h2 id="chrom">9. Chromosomal Distribution</h2>

<p>
    The genome-wide variant density heatmap reveals the distribution of somatic variants
    across chromosomes for each sample. Non-hypermutated samples show relatively uniform,
    low-level variant density.
</p>

<ul>
    <li><strong>S15</strong> shows a striking hotspot on <strong>chr12&ndash;13</strong>, suggesting
        a localized hypermutation event (kataegis) or structural variant affecting this region.</li>
    <li><strong>S04</strong> has elevated variant density across multiple chromosomes
        (chr1, 3, 6, 13, 17), consistent with genome-wide hypermutation.</li>
</ul>

<div class="figure-container">
    <img src="data:image/png;base64,{img_b64('09_chrom_density_heatmap.png')}" alt="Chromosomal Density">
    <div class="figure-caption">
        Figure 9. Chromosomal variant density heatmap. Color intensity represents the number
        of somatic variants per chromosome per sample. S04 and S15 show distinct patterns of
        elevated variant density.
    </div>
</div>

<!-- ================================================================== -->
<h2 id="high-impact">10. High Impact Variants</h2>

<p>
    A total of <strong>{len(high_impact)}</strong> unique gene&ndash;consequence combinations
    were classified as HIGH impact (stop gained, frameshift, splice site). The table below
    shows the most recurrent high-impact events across samples.
</p>

{csv_to_table(high_impact[:25], ['gene', 'consequence', 'n_samples', 'n_variants', 'samples'])}

<p style="font-size:0.85rem; color:#64748b; margin-top:0.5rem;">
    Showing top 25 of {len(high_impact)} high-impact gene&ndash;consequence pairs.
    Full table available in <code>high_impact_variants.csv</code>.
</p>

<!-- ================================================================== -->
<h2 id="caveats">11. Caveats &amp; Next Steps</h2>

<div class="callout callout-alert">
    <strong>Critical: Tumor-Only Somatic Calling</strong>
    All somatic variant calls were generated in <strong>tumor-only mode</strong> (no matched
    normal control). This means the call set contains a substantial proportion of germline
    variants that would normally be subtracted. The VAF peak at ~0.5 and the recurrence of
    highly polymorphic gene families across all samples confirm this contamination. Absolute
    variant counts and TMB values reported here are significantly inflated relative to what
    would be obtained with matched normals.
</div>

<h3>Recommended Next Steps</h3>
<ol>
    <li><strong>Obtain matched normals</strong> &mdash; All 16 samples are confirmed tumor with
        no matched normal controls. Acquiring matched normal tissue from the same mice would
        enable proper germline subtraction and dramatically improve somatic call specificity.</li>
    <li><strong>Panel of Normals (PoN)</strong> &mdash; In the absence of matched normals, construct
        a panel of normals from the same strain background to filter common germline variants.
        This is the most practical path to cleaner somatic calls with the current dataset.</li>
    <li><strong>Mutational signature analysis</strong> &mdash; Run SigProfiler or a similar
        tool to decompose the mutation spectra into COSMIC signatures, especially for
        hypermutated samples S04 and S15.</li>
    <li><strong>Investigate hypermutated samples</strong> &mdash; Check S04 and S15 for defects
        in DNA repair genes (Mlh1, Msh2, Msh6, Pms2, Pole, Pold1).</li>
    <li><strong>Copy number analysis</strong> &mdash; Use the BAM files to call somatic copy
        number alterations with CNVkit or similar tools.</li>
    <li><strong>Structural variant calling</strong> &mdash; Run Delly or Manta on the BAM files
        to identify translocations, inversions, and large-scale rearrangements.</li>
</ol>

<!-- ================================================================== -->
<!-- PART II: FILTERED ANALYSIS                                         -->
<!-- ================================================================== -->

<header style="margin-top:3rem;">
    <h1>Part II: Filtered Somatic Analysis</h1>
    <p>After MGP Germline Subtraction &amp; Cross-Sample Recurrence Filtering</p>
</header>

<h2 id="filtering">12. Germline Filtering Strategy</h2>

<p>
    To address the germline contamination inherent in tumor-only somatic calling,
    a two-pass filtering strategy was applied:
</p>

<ol>
    <li><strong>MGP Allele Frequency Filter</strong> &mdash; Variants with a Mouse Genomes Project
        (MGP) allele frequency &ge; 1% were removed as known strain germline variants. This
        eliminated ~70% of variants in most samples, and &gt;90% in hypermutated samples S04 and S15.</li>
    <li><strong>Cross-Sample Recurrence Filter</strong> &mdash; Variants present in &gt;10 of 16
        samples were removed, as true somatic mutations are unlikely to recur across independent
        tumors. This removed an additional ~8,000 shared variants per sample.</li>
</ol>

<div class="callout callout-info">
    <strong>Filtering Impact</strong>
    On average, <strong>88% of variants were removed</strong> as likely germline or shared artifacts.
    The filtered call set retains ~4,000&ndash;17,000 variants per sample, with TMB values of
    1.5&ndash;6.4 mut/Mb. Because these tumors are on a C57BL/6 background and GRCm39 is derived
    from C57BL/6J, the filtered variants are expected to be predominantly somatic in origin.
</div>

<div class="figure-container">
    <img src="data:image/png;base64,{filt_img_b64('01_filtering_comparison.png')}" alt="Filtering Comparison">
    <div class="figure-caption">
        Figure 10. Impact of germline filtering. Left: variant counts before (gray) and after
        (green) filtering. Right: percentage of variants retained per sample. S04 and S15 had
        the lowest retention rates (5&ndash;6%), indicating their excess variants were predominantly germline.
    </div>
</div>

<h3>Filtered Variant Summary</h3>
{csv_to_table(filtered_variant_data)}

<!-- ================================================================== -->
<h2 id="filt-burden">13. Filtered Mutation Burden</h2>

<p>
    After filtering, S04 and S15 remain elevated (~2&ndash;3x the cohort median) but are no longer
    the dramatic outliers seen in the unfiltered data. S04 retains 14,489 variants (5.3 mut/Mb)
    and S15 retains 17,335 (6.4 mut/Mb), compared to a cohort median of ~5,500 (~2.0 mut/Mb).
    Notably, S14 and S16 also emerge as moderately elevated.
</p>

<div class="figure-container">
    <img src="data:image/png;base64,{filt_img_b64('02_filtered_mutation_burden.png')}" alt="Filtered Mutation Burden">
    <div class="figure-caption">
        Figure 11. Filtered somatic mutation burden. TMB values now range from 1.5&ndash;6.4 mut/Mb,
        consistent with expected mouse tumor somatic mutation rates.
    </div>
</div>

<!-- ================================================================== -->
<h2 id="filt-spectrum">14. Filtered Mutation Spectrum</h2>

<p>
    With germline variants removed, sample-to-sample spectral differences become more apparent.
    The mutation spectrum is more balanced across substitution classes compared to the
    unfiltered data, where germline C&gt;T and T&gt;C transitions dominated.
</p>

<div class="figure-container">
    <img src="data:image/png;base64,{filt_img_b64('04_filtered_mutation_spectrum.png')}" alt="Filtered Mutation Spectrum">
    <div class="figure-caption">
        Figure 12. Filtered somatic SNV mutation spectrum. The proportional differences between
        samples are now more evident, enabling more informative mutational signature analysis.
    </div>
</div>

<!-- ================================================================== -->
<h2 id="filt-impact">15. Filtered Functional Impact</h2>

<div class="figure-container">
    <img src="data:image/png;base64,{filt_img_b64('05_filtered_vep_impact.png')}" alt="Filtered VEP Impact">
    <div class="figure-caption">
        Figure 13. Filtered VEP impact classification. The proportion of HIGH and MODERATE impact
        variants is reduced after filtering, as many protein-coding germline polymorphisms have
        been removed.
    </div>
</div>

<div class="figure-container">
    <img src="data:image/png;base64,{filt_img_b64('06_filtered_top_consequences.png')}" alt="Filtered Top Consequences">
    <div class="figure-caption">
        Figure 14. Top 15 VEP consequence types in filtered variants.
    </div>
</div>

<!-- ================================================================== -->
<h2 id="filt-genes">16. Filtered Recurrently Mutated Genes</h2>

<p>
    After filtering, the recurrently mutated gene list is substantially different from the
    unfiltered analysis. The highly polymorphic Zfp, Vmn2r, and Mroh2a gene families that
    dominated the unfiltered list are largely removed. The remaining recurrent genes are
    more likely to represent genuine somatic cancer-associated events.
</p>

<div class="figure-container">
    <img src="data:image/png;base64,{filt_img_b64('07_filtered_recurrent_genes.png')}" alt="Filtered Recurrent Genes">
    <div class="figure-caption">
        Figure 15. Top recurrently mutated genes after germline filtering (HIGH/MODERATE impact).
        These represent a curated list more likely to contain true somatic events.
    </div>
</div>

<!-- ================================================================== -->
<h2 id="filt-vaf">17. Filtered VAF Distribution</h2>

<p>
    The filtered VAF distribution shows a markedly different profile compared to the unfiltered
    data. While a residual peak at ~0.5 persists (some germline variants without MGP annotation
    remain), the distribution is now broader with substantially more variants at lower VAFs,
    consistent with subclonal somatic mutations.
</p>

<div class="figure-container">
    <img src="data:image/png;base64,{filt_img_b64('08_filtered_vaf_distribution.png')}" alt="Filtered VAF Distribution">
    <div class="figure-caption">
        Figure 16. Filtered VAF distribution. The shift toward lower allele frequencies compared
        to the unfiltered data is consistent with enrichment for true somatic variants.
    </div>
</div>

<!-- ================================================================== -->
<h2 id="filt-chrom">18. Filtered Chromosomal Distribution</h2>

<div class="figure-container">
    <img src="data:image/png;base64,{filt_img_b64('09_filtered_chrom_density.png')}" alt="Filtered Chromosomal Density">
    <div class="figure-caption">
        Figure 17. Chromosomal variant density after filtering. Regional hotspots are more
        clearly resolved without the uniform germline background.
    </div>
</div>

<!-- ================================================================== -->
<h2 id="filt-high">19. Filtered High Impact Variants</h2>

<p>
    After filtering, <strong>{len(filtered_high_impact)}</strong> unique gene&ndash;consequence
    combinations remain classified as HIGH impact. This curated list is substantially smaller
    and more likely to contain genuine somatic loss-of-function events.
</p>

{csv_to_table(filtered_high_impact[:25], ['gene', 'consequence', 'n_samples', 'n_variants', 'samples'])}

<p style="font-size:0.85rem; color:#64748b; margin-top:0.5rem;">
    Full table available in <code>filtered_high_impact_variants.csv</code>.
</p>

<footer>
    WGS Analysis Report &bull; Project 26034-04 &bull; Generated June 30, 2026<br>
    Reference: GRCm39 (Mus musculus) &bull; Pipeline: DRAGEN v13.021 + Ensembl VEP<br>
    Germline filtering: MGP AF &ge; 1% + cross-sample recurrence &gt; 10/16
</footer>

</body>
</html>"""

with open(OUT_PATH, 'w') as f:
    f.write(html)

print(f"Report written to {OUT_PATH}")
print(f"File size: {os.path.getsize(OUT_PATH) / 1024 / 1024:.1f} MB")
