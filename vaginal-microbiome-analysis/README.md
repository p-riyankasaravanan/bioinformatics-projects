# Vaginal Microbiome 16S rRNA Bioinformatics Analysis

> 🚧 **WORK IN PROGRESS**
>
> This is an ongoing independent bioinformatics project. The current repository documents the development and testing of a 16S rRNA amplicon sequencing workflow and a pilot analysis of publicly available vaginal microbiome data. The analysis is being progressively expanded to additional samples and datasets.

## Project overview

This project explores publicly available vaginal microbiome 16S rRNA sequencing data and develops a reproducible workflow for amplicon sequence processing and exploratory microbial community analysis.

The current work includes:

- Public dataset identification and metadata curation
- FASTQ quality assessment using FastQC
- Quality filtering using DADA2
- Sequencing error-rate modelling
- Amplicon sequence variant (ASV) inference
- Paired-read merging
- Chimera removal
- SILVA-based taxonomic assignment
- Exploratory taxonomic and diversity analysis
- Data visualisation using R and ggplot2

## Current status

**Status: Work in progress**

### Completed

- Identified and curated publicly available vaginal microbiome 16S rRNA datasets.
- Harmonised metadata for the Pramanik et al. and Tandon et al. datasets.
- Performed FastQC quality assessment on retrieved paired-end sequencing data.
- Performed DADA2 quality filtering on the Pramanik dataset.
- Developed and tested the complete DADA2 workflow on a representative sample.
- Generated 199 non-chimeric ASVs from the pilot sample.
- Assigned taxonomy using the SILVA reference database.
- Performed exploratory taxonomic composition, richness, dominance and diversity analyses.
- Generated figures documenting the pilot analysis.

### Currently in progress

- Extending the DADA2 workflow to additional samples.
- Expanding downstream analysis beyond the pilot sample.
- Integrating additional sample-level metadata.
- Developing cohort-level comparative analyses.

## Pilot analysis

The current downstream analysis is based on:

**SRA accession:** `SRR12976656`

The pilot sample produced:

| Metric | Result |
|---|---:|
| Non-chimeric ASVs | 199 |
| Phyla | 11 |
| Classes | 17 |
| Orders | 36 |
| Families | 61 |
| Genera | 72 |
| Shannon diversity | 3.935 |
| Simpson diversity | 0.971 |

These measurements are descriptive results from a single pilot sample and should not be interpreted as population-level findings.

## Data sources

Publicly available sequencing data were obtained from the NCBI Sequence Read Archive (SRA).

- Pramanik et al. — **PRJNA674451** — V3–V4 16S rRNA
- Tandon et al. — **PRJNA832047** — V3–V4 16S rRNA

Harmonised metadata are provided in the `metadata/` directory.

## Workflow

```text
Public SRA data
      ↓
FASTQ quality assessment
      ↓
DADA2 quality filtering
      ↓
Error-rate modelling
      ↓
ASV inference
      ↓
Paired-read merging
      ↓
Chimera removal
      ↓
SILVA taxonomic assignment
      ↓
Exploratory analysis
      ↓
Visualisation

vaginal-microbiome-16S-analysis/
│
├── README.md
├── metadata/
├── scripts/
├── results/
├── figures/
└── data/

Tools

R · DADA2 · FastQC · ggplot2 · SILVA · NCBI SRA/ENA
Scope and limitations

This repository represents an ongoing independent project rather than a completed cohort-level study.

The current downstream results are based on a representative pilot sample. Therefore, no cohort-level statistical comparisons, differential abundance analysis, or clinical-group comparisons are claimed at this stage.

Raw FASTQ files and reference databases are not included because of their size. Source BioProject and accession information are provided to identify the publicly available datasets used.

Planned extensions

The project will be extended toward:

Cohort-level ASV inference
Alpha- and beta-diversity analysis across samples
Microbial community composition comparisons
Integration of demographic and clinical metadata
Differential abundance analysis
Chimera removal
      ↓
SILVA taxonomic assignment
      ↓
Exploratory analysis
      ↓
Visualisation
