# Multi-Tissue MR Pipeline - Workflow Guide

## Complete Pipeline Workflow

```
┌─────────────────────────────────────────────────────────────────┐
│                    INPUT DATA                                   │
│         Exposure (eQTL/sQTL)    Outcome (GWAS)                 │
└────────────┬──────────────────────────────────┬─────────────────┘
             │                                  │
             ▼                                  │
    ┌────────────────────────┐                 │
    │ EXPOSURE PREPROCESSING │                 │
    │ - F-statistic filter   │                 │
    │ - Remove indels        │                 │
    │ - Format columns       │                 │
    │ - By chromosome        │                 │
    └──────┬─────────────────┘                 │
           │                                  │
           ▼                                  │
    ┌────────────────────────┐                │
    │ FIND MISSING SNPs      │ ◄──────────────┘
    │ Compare exposure vs    │
    │ outcome SNPs           │
    └──────┬─────────────────┘
           │
           ▼
    ┌────────────────────────┐
    │ PROXY BATCH LDLINK     │
    │ - Query LDlink API     │
    │ - Find LD proxies      │
    │ - R2 > 0.8 threshold   │
    └──────┬─────────────────┘
           │
           ▼
    ┌────────────────────────┐
    │ PROCESS PROXY RESULTS  │
    │ - Filter by R2         │
    │ - Match to outcome     │
    │ - Handle allele flip   │
    └──────┬─────────────────┘
           │
           ▼
    ┌────────────────────────┐
    │ MERGE PROXIED DATA     │
    │ - Combine with outcome │
    │ - Prepare for merge    │
    └──────┬─────────────────┘
           │
           ▼
    ┌────────────────────────┐
    │ OUTCOME PREPROCESSING  │
    │ - Merge proxy data     │
    │ - Remove indels        │
    │ - Format columns       │
    │ - By chromosome        │
    └──────┬─────────────────┘
           │
           ▼
    ┌────────────────────────────────────────┐
    │ MR ANALYSIS (parallel by chromosome)   │
    │ - Harmonize data                       │
    │ - IVW, Egger, Weighted Median          │
    │ - Multiple testing correction          │
    └──────┬──────────────────────────────────┘
           │
           ▼
    ┌────────────────────────┐
    │ COLOCALIZATION ANALYSIS│
    │ - Bayesian coloc test  │
    │ - H0-H4 hypotheses     │
    │ - Filter H4 > 0.8      │
    └──────┬─────────────────┘
           │
           ▼
    ┌────────────────────────┐
    │ RESULTS VISUALIZATION  │
    │ - Manhattan plots      │
    │ - H4 distribution      │
    │ - Summary statistics   │
    │ - Report generation    │
    └──────┬─────────────────┘
           │
           ▼
    ┌─────────────────────────────────────────┐
    │    OUTPUT RESULTS & VISUALIZATIONS      │
    │ exposure/ outcome/ mr_results/          │
    │ colocalization/ plots/                  │
    └─────────────────────────────────────────┘
```

## Key Workflow Steps

### 1. EXPOSURE PREPROCESSING
**Input:** Raw eQTL/sQTL file (gzipped)
**Output:** Chromosome-separated, QC'd exposure files

```
Processing by chromosome:
- F-statistic filtering (default: Fstat > 10)
- Remove indels (keep only ATCG)
- Remove palindromic SNPs (AT/TA, CG/GC)
- Format columns for harmonization
```

### 2. FIND MISSING SNPs
**Input:** Preprocessed exposure, raw outcome
**Output:** List of SNPs in exposure but missing from outcome

```
Logic:
For each SNP in exposure:
  IF SNP NOT IN outcome:
    Add to missing_snps list
    
This identifies which SNPs need proxies
```

### 3. PROXY BATCH LDlink Query
**Input:** Missing SNPs list
**Output:** LDlink proxy search results (combined_query_snp_list_*.txt)

```
For each missing SNP:
  Query LDlink for correlated SNPs
  Filter by R2 > threshold (default: 0.8)
  Save results in batch format
  
Uses LDproxy_batch() from LDlinkR package
Queries 1000 Genomes phase 3 LD data
```

### 4. PROCESS PROXY RESULTS
**Input:** LDlink proxy results, outcome GWAS
**Output:** Processed proxies matching outcome data

```
Filter proxies:
- R2 > 0.8
- Present in outcome GWAS
- Select top proxy per SNP

Handle allele flipping:
- If proxy ≠ original SNP (Distance ≠ 0)
- Flip alleles according to Correlated_Alleles
```

### 5. MERGE PROXIED DATA
**Input:** Processed proxies
**Output:** Outcome data with proxy variants merged

```
Add to outcome GWAS:
- All proxied SNPs
- Original proxy information retained
- Alleles already flipped if needed
```

### 6. OUTCOME PREPROCESSING
**Input:** Outcome GWAS + merged proxy data
**Output:** Chromosome-separated, formatted outcome files

```
Processing by chromosome:
- Include original + proxy variants
- Remove indels
- Remove palindromic SNPs
- Format columns for harmonization
```

### 7. MR ANALYSIS (Parallel by Chromosome)
**Input:** Exposure and outcome by chromosome
**Output:** Harmonized data + MR results

```
Per chromosome:
1. Harmonize data (TwoSampleMR)
   - Match SNPs between exposure/outcome
   - Align alleles
   - Remove discordant alleles

2. Run MR tests:
   - Inverse Variance Weighted (IVW)
   - MR Egger
   - Weighted Median
   - Mode-based estimation

3. Apply multiple testing correction:
   - Bonferroni correction per gene
   - p_threshold = 0.05 / n_unique_genes
```

### 8. COLOCALIZATION ANALYSIS
**Input:** Harmonized data, significant MR results
**Output:** Colocalization test results with H0-H4 posterior probabilities

```
Per significant gene:
1. Extract summary statistics
2. Format for coloc package:
   - SNP positions
   - Effect sizes (beta)
   - Variance of effect (varbeta)
   
3. Run Bayesian colocalization test
   
4. Get posterior probabilities:
   - H0: No association with either trait
   - H1: Association with exposure only
   - H2: Association with outcome only
   - H3: Independent associations
   - H4: Shared causal variant (colocalization)
   
5. Filter results (H4 > 0.8 typically indicates colocalization)
```

### 9. RESULTS VISUALIZATION
**Input:** Combined MR results, colocalization results
**Output:** Publication-ready plots and report

```
Plots generated:
- Manhattan plot: Effect sizes vs p-values
- H4 distribution: Colocalization posterior
- Method comparison: MR method performance
- Hypothesis summary: H0-H4 distribution
- Summary statistics table
- Analysis report (markdown)
```

## Proxy SNP Workflow Details

The proxy batch process is critical for studies where not all exposure SNPs are present in the outcome GWAS. This commonly occurs when:

- Using tissue-specific eQTL data with population GWAS
- Using rare variants not captured in GWAS arrays
- Using different SNP panels across studies

### Why Proxies Matter

When an exposure SNP is missing from outcome:
1. **Information loss**: Cannot directly test causality for that SNP
2. **Power reduction**: Fewer instruments available
3. **Solution**: Use correlated SNPs (proxies) that ARE in outcome

### Proxy Selection Criteria

```
SNP in outcome is a good proxy if:
- In strong LD with original (R2 > 0.8)
- Actually present in outcome GWAS
- Not palindromic (AT/TA, CG/GC)
- Sufficient sample overlap
```

### Allele Handling for Proxies

When proxy ≠ original SNP position:

```
Original:  chr1:1000 A/G
Proxy:     chr1:1050 T/C (R2=0.95 with original)

If effect is aligned to G in original,
must flip alleles in proxy to maintain direction
```

## Parameter Guide

```
--tissue              # Tissue name for labeling
--exposure_file       # eQTL/sQTL file path
--outcome_file        # GWAS summary stats file
--sample_size         # N for F-stat calculation
--fstat_threshold     # Min F-stat for instruments (default: 10)
--p_threshold         # P-value significance threshold (default: 0.05)
--lddlink_token       # LDlink API token (get from ldlink.nci.nih.gov)
--ld_threshold        # Min R2 for proxies (default: 0.8)
--coloc_h4_threshold  # H4 threshold for colocalization (default: 0.8)
--window_size         # ±bp for coloc windows (default: 500000)
--chromosomes         # Chromosomes to process (default: 1..22)
```

## Output Files Structure

```
results/
└── <TISSUE>/
    ├── exposure/
    │   ├── exposure_raw.rds              # Original exposure data
    │   ├── exposure_by_chr_1.rds         # Preprocessed, chr 1
    │   ├── exposure_by_chr_2.rds         # Preprocessed, chr 2
    │   └── ... (chr 3-22)
    │
    ├── proxies/
    │   ├── missing_snps_list.rds         # SNPs needing proxies
    │   ├── missing_snps_summary.txt      # Proxy stats
    │   ├── processed_proxies.rds         # Processed proxy data
    │   ├── proxy_stats.txt               # Proxy matching stats
    │   └── outcome_with_proxies.rds      # Merged outcome data
    │
    ├── outcome/
    │   ├── outcome_raw.rds               # Original outcome data
    │   ├── outcome_chr_1.rds             # Preprocessed, chr 1
    │   ├── outcome_chr_2.rds             # Preprocessed, chr 2
    │   └── ... (chr 3-22)
    │
    ├── mr_results/
    │   ├── harmonized_chr_1.rds          # Harmonized data
    │   ├── mr_raw_chr_1.rds              # All MR results
    │   ├── mr_filtered_chr_1.rds         # Significant results
    │   ├── all_mr_results_combined.rds   # Combined all chr
    │   └── significant_results_combined.rds
    │
    ├── colocalization/
    │   ├── coloc_results_summary.rds     # All coloc results
    │   ├── passing_coloc_genes.rds       # H4 > 0.8
    │   └── coloc_genes_for_proxy.rds     # H4 ≤ 0.8
    │
    ├── plots/
    │   ├── 01_mr_manhattan_plot.pdf
    │   ├── 01_mr_manhattan_plot.png
    │   ├── 02_coloc_h4_distribution.pdf
    │   ├── 02_coloc_h4_distribution.png
    │   ├── 03_mr_method_comparison.pdf
    │   ├── 03_mr_method_comparison.png
    │   ├── 04_coloc_hypothesis_summary.pdf
    │   ├── 04_coloc_hypothesis_summary.png
    │   └── <TISSUE>_summary_stats.txt
    │
    ├── analysis_report_<TISSUE>.md       # Summary report
    │
    ├── execution_report.html             # Nextflow report
    ├── execution_timeline.html           # Execution timeline
    ├── execution_trace.txt               # Task details
    └── pipeline_dag.svg                  # Workflow diagram
```

## Running Multiple Tissues

```bash
# Create parallel runs
for tissue in "Whole_Blood" "Heart" "Liver" "Brain"; do
  nextflow run nextflow_mr_pipeline.nf \
    --tissue "$tissue" \
    --exposure_file "data/eqtl_${tissue}.txt.gz" \
    --outcome_file "data/outcome.tsv" \
    --outdir "results/${tissue}" \
    -profile standard &
done
wait
```

## Performance Notes

- **Proxy batch query**: Most time-consuming step (depends on LDlink API)
- **Colocalization**: Memory-intensive with many genes
- **Parallelization**: Chromosome-level parallelization reduces runtime by ~20x
- **Typical runtime**: 2-4 hours for full analysis (local machine)

## Troubleshooting

| Issue | Cause | Solution |
|-------|-------|----------|
| No proxies found | No missing SNPs or all filtered out | Check LD threshold, exposure/outcome overlap |
| Low proxy success | Few proxies in outcome | Reduce R2 threshold, check SNP notation |
| Memory errors | Large outcome file | Increase memory in config, process by chromosome |
| LDlink timeout | API overloaded | Retry, reduce batch size |
| No coloc results | No significant MR results | Check p_threshold, sample size parameters |

