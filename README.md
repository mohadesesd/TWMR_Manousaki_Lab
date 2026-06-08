# Multi-Tissue Mendelian Randomization (MR) Nextflow Pipeline

A comprehensive Nextflow pipeline for performing Mendelian Randomization analysis across different tissues, including eQTL/sQTL preprocessing, GWAS harmonization, MR testing, and colocalization analysis.

## Table of Contents

- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
- [Parameters](#parameters)
- [Output](#output)
- [Execution Profiles](#execution-profiles)
- [Advanced Usage](#advanced-usage)
- [Troubleshooting](#troubleshooting)

## Features

✅ **Multi-tissue support** - Run analysis for any tissue type (e.g., Whole Blood, Heart, Liver, etc.)

✅ **Automated preprocessing** - Quality control and formatting of exposure and outcome data

✅ **Parallel processing** - Process chromosomes in parallel for improved efficiency

✅ **Multiple MR methods** - Implements IVW, MR Egger, Weighted Median, Mode-based estimation

✅ **Colocalization analysis** - Bayesian colocalization testing with H0-H4 hypotheses

✅ **Comprehensive reporting** - Automatic generation of plots and summary statistics

✅ **Flexible execution** - Support for local, HPC cluster, Docker, and Singularity environments

## Requirements

- **Nextflow** >= 22.10.0
- **R** >= 4.0 with the following packages:
  - TwoSampleMR
  - data.table
  - dplyr
  - tidyr
  - stringr
  - reshape2
  - LDlinkR
  - coloc
  - susieR
  - ggplot2
  - gridExtra

- **Container runtime** (one of):
  - Docker
  - Singularity/Apptainer

### Installation

1. **Install Nextflow**:
```bash
# Install Nextflow (requires Java 11+)
curl -s https://get.nextflow.io | bash
chmod +x nextflow
mv nextflow /usr/local/bin/  # or add to your PATH
```

2. **Clone or download the pipeline**:
```bash[
git clone https://github.com/mohadesesd/TWMR_Manousaki_Lab.git
cd mr-pipeline
```

3. **Install required R packages** (if not using containers):
```R
# In R console
install.packages(c("data.table", "dplyr", "tidyr", "stringr", "reshape2", 
                   "ggplot2", "gridExtra"))

# Install Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("TwoSampleMR", "coloc"))

# Install from GitHub
devtools::install_github("explodecomputer/LDlinkR")
devtools::install_github("stephenslab/susieR")
```

## Quick Start

### Basic execution with default parameters:

```bash
nextflow run nextflow_mr_pipeline.nf \
  --tissue "Whole_Blood" \
  --exposure_file "path/to/exposure.txt.gz" \
  --outcome_file "path/to/outcome.tsv" \
  --outdir "results" \
  -profile standard
```

### Running multiple tissues in sequence:

```bash
# Create a sample sheet
cat > tissues.txt << EOF
Whole_Blood
Heart_Left_Ventricle
Liver
Brain_Cerebellum
EOF

# Run for each tissue
while read tissue; do
  nextflow run nextflow_mr_pipeline.nf \
    --tissue "$tissue" \
    --exposure_file "data/eqtl_${tissue}.txt.gz" \
    --outcome_file "data/outcome_gwas.tsv" \
    --outdir "results/${tissue}" \
    -profile standard \
    -resume
done < tissues.txt
```

## Usage

### Minimal example:

```bash
nextflow run nextflow_mr_pipeline.nf \
  --tissue "MyTissue" \
  --exposure_file "exposure.txt.gz" \
  --outcome_file "outcome.tsv"
```

### Full example with custom parameters:

```bash
nextflow run nextflow_mr_pipeline.nf \
  --tissue "Heart_Left_Ventricle" \
  --exposure_file "/data/eqtl/Heart_LV.txt.gz" \
  --outcome_file "/data/gwas/BMI_GCST90502757.tsv" \
  --sample_size 370 \
  --fstat_threshold 10 \
  --p_threshold 0.05 \
  --ld_threshold 0.8 \
  --coloc_h4_threshold 0.8 \
  --window_size 500000 \
  --outdir "results/Heart_Analysis" \
  -profile cluster \
  -with-report \
  -with-timeline
```

## Parameters

### Required Parameters

| Parameter | Description |
|-----------|-------------|
| `--exposure_file` | Path to exposure data (eQTL/sQTL file) |
| `--outcome_file` | Path to outcome GWAS data |
| `--tissue` | Tissue name for analysis labeling |

### Optional Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--outdir` | results | Output directory for results |
| `--sample_size` | 670 | Sample size for F-statistic calculation |
| `--fstat_threshold` | 10 | Minimum F-statistic for instrument strength |
| `--p_threshold` | 0.05 | P-value threshold for significance (multiple testing corrected) |
| `--lddlink_token` | 23a23... | LDlink API token for proxy searching |
| `--ld_threshold` | 0.8 | Linkage disequilibrium R² threshold for proxies |
| `--coloc_h4_threshold` | 0.8 | H4 posterior probability threshold for colocalization |
| `--window_size` | 500000 | Window size (±bp) around lead SNP for colocalization |
| `--chromosomes` | 1..22 | Chromosomes to process (default: all) |

### Input File Formats

**Exposure file** (eQTL/sQTL data):
```
gene_id     variant_id      maf     pval_nominal    slope   slope_se    ref_factor
ENSG00001   chr1_100_A_T    0.05    0.001           0.5     0.1         1
...
```

**Outcome file** (GWAS data):
```
variant_id  chromosome  base_pair_location  effect_allele   other_allele    eaf     pval    beta    se  sample_size
chr1:100    chr1        100                 T               A               0.3     0.05    0.1     0.05    500000
...
```

## Output

The pipeline generates organized results:

```
results/
├── <TISSUE>/
│   ├── exposure/
│   │   ├── exposure_raw.rds
│   │   └── exposure_by_chr_*.rds
│   ├── outcome/
│   │   ├── outcome_raw.rds
│   │   └── outcome_chr_*.rds
│   ├── mr_results/
│   │   ├── harmonized_chr_*.rds
│   │   ├── mr_raw_chr_*.rds
│   │   ├── mr_filtered_chr_*.rds
│   │   ├── all_mr_results_combined.rds
│   │   └── significant_results_combined.rds
│   ├── colocalization/
│   │   ├── coloc_results_summary.rds
│   │   ├── passing_coloc_genes.rds
│   │   └── coloc_genes_for_proxy.rds
│   ├── plots/
│   │   ├── 01_mr_manhattan_plot.pdf
│   │   ├── 02_coloc_h4_distribution.pdf
│   │   ├── 03_mr_method_comparison.pdf
│   │   ├── 04_coloc_hypothesis_summary.pdf
│   │   └── *_summary_stats.txt
│   └── analysis_report_<TISSUE>.md
├── execution_report.html
├── execution_timeline.html
└── pipeline_dag.svg
```

## Execution Profiles

### Local Execution (standard)

```bash
nextflow run nextflow_mr_pipeline.nf \
  --tissue "Whole_Blood" \
  --exposure_file "exposure.txt.gz" \
  --outcome_file "outcome.tsv" \
  -profile standard
```

### HPC Cluster (SLURM)

```bash
nextflow run nextflow_mr_pipeline.nf \
  --tissue "Whole_Blood" \
  --exposure_file "exposure.txt.gz" \
  --outcome_file "outcome.tsv" \
  -profile cluster
```

### Docker Container

```bash
nextflow run nextflow_mr_pipeline.nf \
  --tissue "Whole_Blood" \
  --exposure_file "exposure.txt.gz" \
  --outcome_file "outcome.tsv" \
  -profile docker
```

### Singularity Container

```bash
nextflow run nextflow_mr_pipeline.nf \
  --tissue "Whole_Blood" \
  --exposure_file "exposure.txt.gz" \
  --outcome_file "outcome.tsv" \
  -profile singularity
```

## Advanced Usage

### Resume failed runs

```bash
# Continue from last successful step
nextflow run nextflow_mr_pipeline.nf \
  --tissue "Whole_Blood" \
  --exposure_file "exposure.txt.gz" \
  --outcome_file "outcome.tsv" \
  -resume
```

### Dry run (preview without execution)

```bash
nextflow run nextflow_mr_pipeline.nf \
  --tissue "Whole_Blood" \
  --exposure_file "exposure.txt.gz" \
  --outcome_file "outcome.tsv" \
  -preview
```

### Generate execution DAG visualization

```bash
nextflow run nextflow_mr_pipeline.nf \
  --tissue "Whole_Blood" \
  --exposure_file "exposure.txt.gz" \
  --outcome_file "outcome.tsv" \
  -with-dag workflow_dag.html
```

### Parallel processing for multiple tissues

```bash
# Use GNU parallel
cat tissues.txt | parallel -j 3 "nextflow run nextflow_mr_pipeline.nf \
  --tissue {} \
  --exposure_file data/eqtl_{}.txt.gz \
  --outcome_file data/outcome.tsv \
  --outdir results/{}"
```

### Custom configuration

Create a custom `my_config.config`:

```groovy
// Custom configuration
process {
    memory = '32 GB'
    cpus = 8
    time = '24h'
    
    withName: 'COLOCALIZATION' {
        memory = '64 GB'
        cpus = 16
        time = '48h'
    }
}

// Custom container images
docker {
    runOptions = '-u $(id -u):$(id -g) --entrypoint /bin/bash'
}

// Custom output organization
params {
    outdir = "/scratch/results"
    save_intermediate_files = true
}
```

Run with custom config:

```bash
nextflow run nextflow_mr_pipeline.nf \
  --tissue "Whole_Blood" \
  --exposure_file "exposure.txt.gz" \
  --outcome_file "outcome.tsv" \
  -c my_config.config
```

## Troubleshooting

### Issue: "Container image not found"

**Solution**: Pull the container image first:
```bash
docker pull rocker/tidyverse:latest
# OR
singularity pull docker://rocker/tidyverse:latest
```

### Issue: "LDlinkR token invalid"

**Solution**: Update your LDlinkR token:
```bash
nextflow run nextflow_mr_pipeline.nf \
  --tissue "Whole_Blood" \
  --exposure_file "exposure.txt.gz" \
  --outcome_file "outcome.tsv" \
  --lddlink_token "your_valid_token"
```

Get a token from: https://ldlink.nci.nih.gov/

### Issue: Out of memory errors

**Solution**: Increase memory allocation:
```bash
nextflow run nextflow_mr_pipeline.nf \
  --tissue "Whole_Blood" \
  --exposure_file "exposure.txt.gz" \
  --outcome_file "outcome.tsv" \
  -profile cluster  # Uses higher memory settings
```

Or modify `nextflow.config` with custom resource limits.

### Issue: No harmonized SNPs found

**Possible causes**:
- Exposure and outcome datasets are not overlapping
- SNP naming conventions don't match (CHR:POS vs rsid)
- Alleles are misaligned (check strand flipping)

**Solution**:
- Verify input file formats
- Check chromosome naming (chr1 vs 1)
- Examine sample QC steps

### Accessing help

```bash
nextflow run nextflow_mr_pipeline.nf --help
```


## License

MIT License - See LICENSE file for details

## Support

For issues, questions, or contributions:
- Open an issue on GitHub
- Contact: mohadese.sayahiandehkordi@mail.mcgill.ca
