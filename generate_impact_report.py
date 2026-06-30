#!/usr/bin/env python3
"""Generate a biological impact report for WGS somatic mutations in sucrose-fed mice."""

import base64
import os
import csv

FIGURES_DIR = os.path.expanduser('~/wgs-project/figures_filtered')
OUT_PATH = os.path.expanduser('~/wgs-project/WGS_Biological_Impact_Report.html')


def img_b64(filename):
    with open(os.path.join(FIGURES_DIR, filename), 'rb') as f:
        return base64.b64encode(f.read()).decode()


def read_csv_from(directory, filename):
    with open(os.path.join(directory, filename)) as f:
        return list(csv.DictReader(f))


filtered_variant_data = read_csv_from(FIGURES_DIR, 'filtered_variant_summary.csv')
filtered_high_impact = read_csv_from(FIGURES_DIR, 'filtered_high_impact_variants.csv')


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


html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Biological Impact Report — Somatic Mutations in Sucrose-Fed C57BL/6J Mice</title>
<style>
:root {{
    --bg: #ffffff;
    --fg: #1a1a2e;
    --accent: #0e7c61;
    --accent-light: #d1fae5;
    --border: #e2e8f0;
    --section-bg: #f8fafc;
    --highlight: #fef3c7;
    --red: #dc2626;
    --orange: #ea580c;
}}
* {{ margin: 0; padding: 0; box-sizing: border-box; }}
body {{
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
    color: var(--fg);
    background: var(--bg);
    line-height: 1.7;
    max-width: 1100px;
    margin: 0 auto;
    padding: 2rem;
}}
header {{
    text-align: center;
    padding: 2.5rem 1.5rem;
    margin-bottom: 2rem;
    background: linear-gradient(135deg, #064e3b 0%, #0e7c61 100%);
    color: white;
    border-radius: 12px;
}}
header h1 {{ font-size: 1.8rem; margin-bottom: 0.5rem; font-weight: 700; }}
header p {{ font-size: 1rem; opacity: 0.9; }}
header .meta {{ font-size: 0.85rem; opacity: 0.7; margin-top: 0.5rem; }}
h2 {{
    font-size: 1.4rem;
    color: var(--accent);
    border-bottom: 2px solid var(--accent);
    padding-bottom: 0.4rem;
    margin: 2.5rem 0 1rem;
}}
h3 {{ font-size: 1.15rem; margin: 1.5rem 0 0.5rem; color: #334155; }}
p, li {{ margin-bottom: 0.6rem; }}
ul, ol {{ padding-left: 1.5rem; }}
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
}}
tbody tr:nth-child(even) {{ background: var(--section-bg); }}
tbody tr:hover {{ background: var(--accent-light); }}
.callout {{
    padding: 1rem 1.25rem;
    border-radius: 8px;
    margin: 1.5rem 0;
    border-left: 4px solid;
}}
.callout-key {{ background: #d1fae5; border-color: #0e7c61; }}
.callout-warning {{ background: #fef3c7; border-color: #f59e0b; }}
.callout-insight {{ background: #dbeafe; border-color: #2563eb; }}
.callout-alert {{ background: #fee2e2; border-color: #dc2626; }}
.callout strong {{ display: block; margin-bottom: 0.3rem; }}
.gene-card {{
    background: var(--section-bg);
    border: 1px solid var(--border);
    border-radius: 8px;
    padding: 1rem 1.25rem;
    margin: 0.75rem 0;
}}
.gene-card .gene-name {{
    font-weight: 700;
    color: var(--accent);
    font-size: 1.05rem;
}}
.gene-card .gene-full {{ font-size: 0.9rem; color: #64748b; font-style: italic; }}
.gene-card .gene-detail {{ margin-top: 0.4rem; }}
.highlight {{ background: var(--highlight); padding: 0.1rem 0.3rem; border-radius: 3px; }}
.stats-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
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
.two-col {{
    display: grid;
    grid-template-columns: 1fr 1fr;
    gap: 1.5rem;
    margin: 1rem 0;
}}
@media (max-width: 768px) {{
    .two-col {{ grid-template-columns: 1fr; }}
}}
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
    .figure-container {{ break-inside: avoid; }}
    .gene-card {{ break-inside: avoid; }}
}}
</style>
</head>
<body>

<header>
    <h1>Biological Impact of Somatic Mutations</h1>
    <p>Tumors from C57BL/6J Mice with Ad Libitum Sucrose Water</p>
    <div class="meta">
        16 Tumor Samples &bull; WGS on GRCm39 &bull; Filtered Somatic Variants &bull; June 30, 2026
    </div>
</header>

<div class="toc">
    <h3>Table of Contents</h3>
    <ol>
        <li><a href="#summary">Executive Summary</a></li>
        <li><a href="#context">Experimental Context</a></li>
        <li><a href="#landscape">Mutational Landscape</a></li>
        <li><a href="#spectrum-interp">Mutation Spectrum &amp; Mutational Processes</a></li>
        <li><a href="#hypermut">Hypermutated Samples: S04 &amp; S15</a></li>
        <li><a href="#drivers">Recurrent Somatic Events</a></li>
        <li><a href="#private">Sample-Private High-Impact Mutations</a></li>
        <li><a href="#pathways">Pathway-Level Interpretation</a></li>
        <li><a href="#sucrose">Relevance to Sucrose Exposure</a></li>
        <li><a href="#limitations">Limitations</a></li>
        <li><a href="#next">Recommended Next Steps</a></li>
    </ol>
</div>


<!-- ================================================================== -->
<h2 id="summary">1. Executive Summary</h2>

<div class="stats-grid">
    <div class="stat-card">
        <div class="value">1.5&ndash;6.4</div>
        <div class="label">TMB (mut/Mb)</div>
    </div>
    <div class="stat-card">
        <div class="value">0</div>
        <div class="label">Classic Driver Mutations</div>
    </div>
    <div class="stat-card">
        <div class="value">2 / 16</div>
        <div class="label">Hypermutated Samples</div>
    </div>
    <div class="stat-card">
        <div class="value">37</div>
        <div class="label">High Impact Gene Events</div>
    </div>
</div>

<p>
    Whole genome sequencing of 16 tumors from sucrose-fed C57BL/6J mice reveals a moderate
    somatic mutation burden (median TMB ~2.0 mut/Mb after germline filtering) with
    <strong>no mutations in canonical cancer driver genes</strong> (Trp53, Kras, Pten, Rb1,
    Apc, Braf, etc.). Two samples (S04 and S15) display a hypermutation phenotype
    characterized by a striking <strong>excess of small indels</strong> (67&ndash;75% of
    filtered variants), suggestive of replication slippage or mismatch repair stress.
</p>

<p>
    The dominant mutation signature across all samples is C&gt;T transitions (~28&ndash;42%),
    consistent with spontaneous cytosine deamination (COSMIC SBS1), a clock-like process that
    accumulates with cell divisions. The absence of classic oncogenic drivers suggests these
    tumors may be driven by non-genetic mechanisms (epigenetic, metabolic, or inflammatory)
    with the observed somatic mutations representing passenger events accumulated during
    clonal expansion.
</p>


<!-- ================================================================== -->
<h2 id="context">2. Experimental Context</h2>

<p>
    C57BL/6J mice were provided <strong>sucrose water ad libitum</strong> and subsequently
    developed tumors. This experimental design positions the study at the intersection of
    <strong>dietary sugar metabolism and tumorigenesis</strong>. High sucrose intake in rodent
    models is associated with:
</p>

<ul>
    <li><strong>Chronic hyperinsulinemia</strong> &mdash; elevated insulin and IGF-1 signaling
        promote cell proliferation and suppress apoptosis via PI3K/Akt/mTOR</li>
    <li><strong>Enhanced Warburg effect</strong> &mdash; high glucose availability can select
        for glycolysis-dependent tumor phenotypes</li>
    <li><strong>Oxidative stress</strong> &mdash; excess fructose metabolism generates reactive
        oxygen species (ROS), contributing to DNA damage</li>
    <li><strong>Chronic low-grade inflammation</strong> &mdash; elevated pro-inflammatory
        cytokines (TNF-alpha, IL-6) in adipose tissue create a tumor-promoting microenvironment</li>
    <li><strong>Dysbiosis</strong> &mdash; altered gut microbiome composition affecting systemic
        immune surveillance and metabolite production</li>
    <li><strong>De novo lipogenesis</strong> &mdash; increased hepatic fat synthesis via
        SREBP-1c/FASN axis, potentially supporting membrane biogenesis in proliferating cells</li>
</ul>


<!-- ================================================================== -->
<h2 id="landscape">3. Mutational Landscape</h2>

<p>
    After removing germline variants (MGP allele frequency filtering + cross-sample recurrence),
    the filtered somatic mutation burden ranges from <strong>1.49 mut/Mb (S11)</strong> to
    <strong>6.35 mut/Mb (S15)</strong>. The majority of samples cluster around 1.5&ndash;2.5
    mut/Mb, which is within the expected range for murine tumors without known exogenous
    mutagen exposure.
</p>

<div class="figure-container">
    <img src="data:image/png;base64,{img_b64('02_filtered_mutation_burden.png')}" alt="Filtered Mutation Burden">
    <div class="figure-caption">
        Figure 1. Filtered somatic mutation burden across all 16 tumor samples.
        S04 and S15 show elevated TMB at 5.3 and 6.4 mut/Mb respectively.
    </div>
</div>

<h3>Variant Composition</h3>
<p>
    A notable feature of this cohort is the <strong>high proportion of indels</strong> relative
    to SNVs. In typical WGS somatic calling, SNVs outnumber indels by 5&ndash;10:1. In this
    dataset, the SNV:indel ratio is approximately 1:1 for most samples, and drops to
    <strong>1:2 for S04</strong> and <strong>1:2 for S15</strong>. This indel excess is
    enriched for short (&le;5 bp) insertions and deletions, consistent with
    <strong>replication slippage at microsatellite loci</strong>.
</p>

<div class="figure-container">
    <img src="data:image/png;base64,{img_b64('03_filtered_variant_types.png')}" alt="Filtered Variant Types">
    <div class="figure-caption">
        Figure 2. Variant type composition after filtering. Note the elevated indel fraction
        across all samples, with S04 and S15 showing a dominant indel phenotype.
    </div>
</div>

{csv_to_table(filtered_variant_data)}


<!-- ================================================================== -->
<h2 id="spectrum-interp">4. Mutation Spectrum &amp; Mutational Processes</h2>

<p>
    The six-class SNV substitution spectrum provides insight into the mutational processes
    active in these tumors.
</p>

<div class="figure-container">
    <img src="data:image/png;base64,{img_b64('04_filtered_mutation_spectrum.png')}" alt="Filtered Mutation Spectrum">
    <div class="figure-caption">
        Figure 3. Filtered somatic SNV mutation spectrum showing absolute counts and
        relative proportions per sample.
    </div>
</div>

<h3>Dominant Signatures</h3>

<div class="callout callout-insight">
    <strong>C&gt;T Transitions (28&ndash;42% of SNVs)</strong>
    The dominant substitution class across all samples, consistent with <strong>COSMIC SBS1</strong>
    (spontaneous deamination of 5-methylcytosine). SBS1 is a clock-like signature that
    correlates with the number of cell divisions, not a specific mutagen. Its dominance
    suggests these tumors have undergone substantial clonal expansion, potentially accelerated
    by the proliferative effects of chronic hyperinsulinemia.
</div>

<div class="callout callout-insight">
    <strong>T&gt;C Transitions (19&ndash;30% of SNVs)</strong>
    The second most common substitution. Elevated T&gt;C transitions are associated with
    <strong>COSMIC SBS5</strong>, another clock-like signature of unknown etiology that
    accumulates with age. Together, the SBS1 + SBS5 pattern indicates these tumors arose
    through endogenous mutational processes rather than exogenous carcinogen exposure.
</div>

<h3>Sample-Specific Patterns</h3>

<div class="two-col">
    <div class="gene-card">
        <div class="gene-name">S16: Elevated C&gt;T (42.1%)</div>
        <div class="gene-detail">
            The highest C&gt;T proportion in the cohort, potentially indicating enhanced
            deamination or a more pronounced SBS1/aging signature. S16 also carries the
            most HIGH impact mutations (9) and affected genes (43) among non-hypermutated
            samples.
        </div>
    </div>
    <div class="gene-card">
        <div class="gene-name">S14: Elevated T&gt;G (18.2%)</div>
        <div class="gene-detail">
            S14 shows an unusual enrichment of T&gt;G transversions (nearly double the cohort
            average of ~9.5%). T&gt;G mutations can be associated with <strong>oxidative DNA
            damage</strong> (8-oxoguanine lesions), which is notable in the context of
            fructose-induced ROS production from chronic sucrose intake.
        </div>
    </div>
</div>


<!-- ================================================================== -->
<h2 id="hypermut">5. Hypermutated Samples: S04 &amp; S15</h2>

<div class="callout callout-alert">
    <strong>Indel-Dominant Hypermutation Phenotype</strong>
    S04 and S15 retain 2&ndash;3x more variants than the cohort median even after aggressive
    germline filtering. Critically, their excess is overwhelmingly composed of
    <strong>indels</strong> (74.5% and 67.5% respectively), not SNVs. This is the hallmark
    of <strong>microsatellite instability (MSI)</strong>.
</div>

<table>
<thead><tr><th>Sample</th><th>Total</th><th>SNVs</th><th>Indels</th><th>Indel %</th><th>Short (&le;5bp)</th><th>Long (&gt;5bp)</th></tr></thead>
<tbody>
<tr><td>S01 (typical)</td><td>4,779</td><td>2,496</td><td>2,283</td><td>47.8%</td><td>1,722</td><td>561</td></tr>
<tr style="background:#fee2e2;"><td>S04</td><td>14,489</td><td>3,697</td><td>10,792</td><td>74.5%</td><td>7,860</td><td>2,932</td></tr>
<tr style="background:#fee2e2;"><td>S15</td><td>17,335</td><td>5,632</td><td>11,703</td><td>67.5%</td><td>8,206</td><td>3,497</td></tr>
<tr><td>S16 (elevated)</td><td>9,031</td><td>4,893</td><td>4,138</td><td>45.8%</td><td>3,354</td><td>784</td></tr>
</tbody>
</table>

<h3>Biological Interpretation</h3>

<p>
    The short indel excess in S04 and S15 is consistent with <strong>replication slippage
    errors at repetitive sequences</strong> that were not corrected by mismatch repair (MMR).
    While no coding mutations in canonical MMR genes (Mlh1, Msh2, Msh6, Pms2) were detected,
    MSI can also arise from:
</p>

<ul>
    <li><strong>Epigenetic silencing</strong> of MMR genes (e.g., Mlh1 promoter methylation) &mdash;
        not detectable by WGS alone</li>
    <li><strong>Replication stress</strong> from rapid proliferation outpacing repair capacity &mdash;
        potentially exacerbated by insulin-driven growth signaling</li>
    <li><strong>Metabolic disruption of nucleotide pools</strong> &mdash; high glucose/fructose
        flux through the pentose phosphate pathway can alter dNTP balance, increasing
        replication errors</li>
</ul>

<div class="callout callout-key">
    <strong>Connection to Sucrose Exposure</strong>
    The MSI-like phenotype in S04 and S15 is potentially the most significant finding linking
    the sucrose diet to tumor mutagenesis. Chronic sucrose intake drives hyperinsulinemia and
    rapid cell proliferation, which increases replication stress. Combined with fructose-mediated
    oxidative damage and potential epigenetic changes from metabolic reprogramming, this could
    create conditions for microsatellite instability even without genetic MMR deficiency.
</div>


<!-- ================================================================== -->
<h2 id="drivers">6. Recurrent Somatic Events</h2>

<p>
    No mutations were detected in canonical oncogenes or tumor suppressors (Trp53, Kras, Pten,
    Rb1, Apc, Braf, Egfr, Myc, Pik3ca, Nf1, etc.). The most recurrently affected genes with
    HIGH or MODERATE impact are:
</p>

<div class="figure-container">
    <img src="data:image/png;base64,{img_b64('07_filtered_recurrent_genes.png')}" alt="Filtered Recurrent Genes">
    <div class="figure-caption">
        Figure 4. Top recurrently mutated genes after germline filtering.
    </div>
</div>

<h3>Genes of Potential Biological Interest</h3>

<div class="gene-card">
    <div class="gene-name">Chd5 <span class="gene-full">(Chromodomain Helicase DNA Binding Protein 5)</span></div>
    <div class="gene-detail">
        Mutated in <strong>4 samples</strong> (missense variants). Chd5 is a well-established
        <strong>tumor suppressor</strong> involved in chromatin remodeling, p53-mediated apoptosis,
        and cell cycle regulation. Loss of Chd5 function has been implicated in neuroblastoma,
        breast cancer, and colorectal cancer. Its recurrence across 25% of samples is notable
        and may represent a genuine somatic selection event in these tumors.
    </div>
</div>

<div class="gene-card">
    <div class="gene-name">Pfdn2 <span class="gene-full">(Prefoldin Subunit 2)</span></div>
    <div class="gene-detail">
        Splice donor variant in <strong>9 samples</strong> (56%). Prefoldin is a molecular
        chaperone involved in cytoskeletal protein folding (actin and tubulin). Disruption
        of prefoldin function can impair cytoskeletal integrity and mitotic fidelity, potentially
        contributing to chromosomal instability. Its high recurrence warrants investigation
        as a potential early event.
    </div>
</div>

<div class="gene-card">
    <div class="gene-name">Cdv3 <span class="gene-full">(Carnitine Deficiency-Associated Gene 3)</span></div>
    <div class="gene-detail">
        Frameshift in <strong>8 samples</strong> (50%). Cdv3 is involved in
        <strong>cell proliferation</strong> and interacts with the c-Myc transcriptional network.
        Its association with carnitine metabolism is relevant to the metabolic context of this
        study, as carnitine is critical for fatty acid oxidation &mdash; a pathway that competes
        with glycolysis for cellular energy production.
    </div>
</div>

<div class="gene-card">
    <div class="gene-name">Prpf40a <span class="gene-full">(Pre-mRNA Processing Factor 40 Homolog A)</span></div>
    <div class="gene-detail">
        Splice region missense in <strong>8 samples</strong> (50%). PRPF40A is a component of
        the spliceosome involved in pre-mRNA splicing. Recurrent splicing factor mutations are
        increasingly recognized in cancer, where they can cause widespread alternative splicing
        and generate tumor-promoting transcript isoforms.
    </div>
</div>

<div class="gene-card">
    <div class="gene-name">Rhbdf2 <span class="gene-full">(Rhomboid Family Member 2 / iRhom2)</span></div>
    <div class="gene-detail">
        Frameshift in <strong>S16</strong>. RHBDF2 regulates the maturation and trafficking
        of ADAM17 (TACE), which processes EGFR ligands and TNF-alpha. Loss of RHBDF2 function
        affects <strong>EGFR signaling and inflammatory cytokine release</strong> &mdash;
        directly relevant to both tumorigenesis and the inflammatory consequences of high
        sugar intake.
    </div>
</div>

<div class="gene-card">
    <div class="gene-name">Helz2 <span class="gene-full">(Helicase with Zinc Finger 2)</span></div>
    <div class="gene-detail">
        Frameshift in <strong>S16</strong>. HELZ2 functions as a transcriptional coactivator
        of <strong>PPAR-gamma</strong>, a master regulator of adipogenesis and lipid metabolism.
        PPAR-gamma also plays roles in glucose homeostasis, inflammation, and cell differentiation.
        Disruption of this coactivator in the context of chronic sucrose feeding could impact
        metabolic reprogramming of tumor cells.
    </div>
</div>

<div class="gene-card">
    <div class="gene-name">Casp2 <span class="gene-full">(Caspase-2)</span></div>
    <div class="gene-detail">
        Frameshift in <strong>S07</strong>. Caspase-2 is a unique initiator caspase with
        tumor-suppressive functions. It responds to metabolic stress including
        <strong>oxidative damage</strong> and <strong>DNA damage from replication stress</strong>.
        Recent studies show Caspase-2 acts as a metabolic sensor that triggers apoptosis in
        response to excess lipid accumulation, making its loss particularly relevant in a
        high-sugar feeding model.
    </div>
</div>

<div class="gene-card">
    <div class="gene-name">Setx <span class="gene-full">(Senataxin)</span></div>
    <div class="gene-detail">
        Frameshift in <strong>S03</strong>. Senataxin is a DNA/RNA helicase essential for
        <strong>R-loop resolution</strong> during transcription. R-loop accumulation causes
        replication-transcription conflicts and genomic instability. SETX deficiency increases
        sensitivity to oxidative stress, connecting to the ROS-generating effects of fructose
        metabolism.
    </div>
</div>

<div class="gene-card">
    <div class="gene-name">Il12rb2 <span class="gene-full">(Interleukin-12 Receptor Beta 2)</span></div>
    <div class="gene-detail">
        Stop-gained in <strong>S15</strong>. IL-12Rbeta2 is critical for IL-12 signaling, which
        drives Th1 immune responses and NK cell activation &mdash; key anti-tumor immune
        mechanisms. Loss of IL-12 receptor function could impair <strong>immune surveillance</strong>
        of the tumor.
    </div>
</div>


<!-- ================================================================== -->
<h2 id="private">7. Sample-Private High-Impact Mutations</h2>

<p>
    Several high-impact mutations are private to individual samples, potentially reflecting
    unique tumor evolution trajectories:
</p>

{csv_to_table(filtered_high_impact, ['gene', 'consequence', 'n_samples', 'n_variants', 'samples'])}


<!-- ================================================================== -->
<h2 id="pathways">8. Pathway-Level Interpretation</h2>

<p>
    While no single canonical pathway is recurrently disrupted by driver mutations, the
    aggregate mutational profile suggests several themes:
</p>

<h3>8.1 Chromatin Remodeling &amp; Epigenetic Regulation</h3>
<p>
    Chd5 (4 samples) is a chromatin remodeler with tumor-suppressive function. Combined with
    the MSI-like phenotype in S04/S15 (which can arise from epigenetic silencing of MMR genes),
    this suggests <strong>epigenetic dysregulation</strong> may be an important feature of
    tumor development in this model.
</p>

<h3>8.2 RNA Processing &amp; Splicing</h3>
<p>
    Recurrent mutations in Prpf40a (8 samples, spliceosome) and Polr3d (RNA Polymerase III)
    suggest disruption of RNA processing. Splicing factor mutations are an emerging category
    of cancer drivers that can generate tumor-promoting transcript diversity.
</p>

<h3>8.3 Protein Quality Control</h3>
<p>
    Pfdn2 (9 samples) functions in protein folding. The prefoldin complex is critical for
    maintaining proteostasis under proliferative stress. Its disruption may reflect the
    <strong>increased protein synthesis demands</strong> of rapidly dividing tumor cells
    fueled by high glucose availability.
</p>

<h3>8.4 Metabolic Reprogramming</h3>
<p>
    Cdv3 (carnitine metabolism, 8 samples), Helz2 (PPAR-gamma coactivator, S16), and the
    broader metabolic context suggest tumor cells may have undergone metabolic adaptations.
    Carnitine is essential for fatty acid beta-oxidation; its disruption in the context of
    excess dietary sugar may force cells toward glycolytic dependence (Warburg effect),
    providing a survival advantage for tumor cells.
</p>

<h3>8.5 Immune Evasion</h3>
<p>
    Il12rb2 loss (S15) and Csf1r mutation (S14) affect immune signaling pathways. The
    tumor microenvironment in chronically inflamed, metabolically dysregulated tissue may
    select for immune evasion mutations.
</p>

<h3>8.6 DNA Damage Response</h3>
<p>
    Setx (R-loop resolution, S03) and Casp2 (metabolic stress-induced apoptosis, S07)
    are involved in detecting and responding to DNA damage. Their loss could increase
    tolerance for genomic instability, permitting the accumulation of further mutations.
</p>


<!-- ================================================================== -->
<h2 id="sucrose">9. Relevance to Sucrose Exposure</h2>

<div class="callout callout-key">
    <strong>Key Findings in the Context of Ad Libitum Sucrose</strong>
    <ol style="margin-top:0.5rem;">
        <li>The <strong>absence of classic driver mutations</strong> suggests sucrose-associated
            tumorigenesis may operate through non-mutational mechanisms (metabolic reprogramming,
            chronic inflammation, epigenetic changes) rather than direct genotoxicity.</li>
        <li>The <strong>clock-like mutation signatures</strong> (SBS1/SBS5 dominance) indicate
            accelerated cell proliferation rather than chemical mutagenesis &mdash; consistent
            with insulin/IGF-1-driven growth.</li>
        <li>The <strong>MSI-like indel phenotype</strong> in S04/S15 may represent replication
            stress from rapid proliferation outpacing DNA repair, exacerbated by metabolic
            disruption of nucleotide pools.</li>
        <li>The <strong>T&gt;G signature elevation in S14</strong> (18.2%) is consistent with
            oxidative DNA damage, potentially from fructose-induced ROS production.</li>
        <li><strong>Metabolically relevant gene mutations</strong> (Cdv3/carnitine metabolism,
            Helz2/PPAR-gamma, Casp2/metabolic stress sensing) suggest selective pressure
            related to the metabolic environment created by chronic sucrose feeding.</li>
    </ol>
</div>

<h3>Proposed Model</h3>
<p>
    Based on the mutational evidence, the following model of sucrose-associated tumorigenesis
    is supported:
</p>
<ol>
    <li><strong>Initiation:</strong> Chronic sucrose intake leads to hyperinsulinemia, elevated
        IGF-1, and increased cell proliferation rates. This accelerates the accumulation of
        endogenous mutations (SBS1 deamination) through increased cell divisions.</li>
    <li><strong>Promotion:</strong> Fructose metabolism generates ROS (seen as T&gt;G mutations),
        while chronic inflammation creates a tumor-promoting microenvironment. Epigenetic
        changes (potentially including MMR gene silencing) contribute to genomic instability.</li>
    <li><strong>Progression:</strong> Tumors adapt through metabolic reprogramming (Cdv3, Helz2
        mutations favoring glycolytic dependence), evasion of metabolic stress checkpoints
        (Casp2 loss), and immune escape (Il12rb2 loss). Chromatin remodeling defects (Chd5)
        may enable transcriptional plasticity.</li>
</ol>


<!-- ================================================================== -->
<h2 id="limitations">10. Limitations</h2>

<ul>
    <li><strong>Tumor-only calling:</strong> Without matched normal tissue, some residual
        germline variants may persist despite filtering. The B6J background on GRCm39
        mitigates but does not eliminate this concern.</li>
    <li><strong>No control group data:</strong> Without WGS from tumors in non-sucrose-fed mice,
        it is not possible to determine which mutational features are specific to sucrose
        exposure versus baseline tumor mutagenesis in C57BL/6J.</li>
    <li><strong>Tumor type unknown:</strong> The tissue of origin and histological classification
        of these tumors are not specified, limiting pathway interpretation.</li>
    <li><strong>Epigenetic data absent:</strong> The hypothesis of epigenetic MMR silencing
        cannot be tested with WGS data alone; bisulfite sequencing or methylation arrays
        would be needed.</li>
    <li><strong>No transcriptomic data:</strong> The functional impact of splice-site and
        regulatory mutations cannot be confirmed without RNA-seq.</li>
    <li><strong>Formal signature analysis pending:</strong> SigProfiler decomposition has not
        been performed; spectrum interpretations are based on the six-class SNV profile only.</li>
</ul>


<!-- ================================================================== -->
<h2 id="next">11. Recommended Next Steps</h2>

<ol>
    <li><strong>Confirm tumor histology</strong> &mdash; Pathological classification of the
        tumor type would greatly aid biological interpretation.</li>
    <li><strong>Formal mutational signature analysis</strong> &mdash; Run SigProfiler or
        MutationalPatterns on the filtered VCFs to decompose into COSMIC SBS, DBS, and ID
        signatures. The ID signature catalog is particularly relevant for characterizing the
        indel-dominant phenotype in S04/S15.</li>
    <li><strong>MSI testing</strong> &mdash; Formally assess microsatellite instability status
        using MSIsensor or similar tools on the BAM files.</li>
    <li><strong>Methylation analysis</strong> &mdash; Investigate MLH1 promoter methylation
        in S04/S15 to test the epigenetic MMR silencing hypothesis.</li>
    <li><strong>Include a control cohort</strong> &mdash; WGS of tumors (if any develop) from
        age-matched C57BL/6J mice on normal water would establish baseline mutation rates
        and distinguish sucrose-specific effects.</li>
    <li><strong>RNA-seq integration</strong> &mdash; Transcriptomic profiling would reveal
        whether splicing mutations (Pfdn2, Prpf40a) produce functional transcript changes
        and whether metabolic gene expression is altered.</li>
    <li><strong>Copy number analysis</strong> &mdash; Run CNVkit on the BAM files to identify
        chromosomal gains/losses that may harbor driver events not captured by SNV/indel
        calling.</li>
    <li><strong>Pathway enrichment</strong> &mdash; Perform Gene Ontology and KEGG pathway
        enrichment on all coding mutations to identify over-represented biological processes.</li>
</ol>


<footer>
    Biological Impact Report &bull; Project 26034-04 &bull; June 30, 2026<br>
    C57BL/6J Tumor Samples &bull; Ad Libitum Sucrose Water &bull; GRCm39 + DRAGEN v13.021 + VEP
</footer>

</body>
</html>"""

with open(OUT_PATH, 'w') as f:
    f.write(html)

print(f"Report written to {OUT_PATH}")
print(f"File size: {os.path.getsize(OUT_PATH) / 1024 / 1024:.1f} MB")
