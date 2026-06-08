# Proxy Batch Module Integration Guide

## What's New

You were absolutely right! I've now added the **complete proxy batch processing workflow** that was missing from the initial pipeline.

### New Components Added

#### 1. **proxy_batch.nf** - Complete proxy processing module
This module includes 5 processes:

```
FIND_MISSING_SNPS
    ↓
PROXY_BATCH_LDLINK
    ↓
PROCESS_PROXY_RESULTS
    ↓
MERGE_PROXIED_DATA
    ↓
OUTCOME_PREPROCESSING
```

### Process Details

#### FIND_MISSING_SNPS
- Identifies SNPs present in exposure but absent from outcome
- Generates missing SNPs list and summary statistics
- Output: `missing_snps_list.rds`

#### PROXY_BATCH_LDLINK
- Queries LDlink's batch proxy search API
- Finds all correlated SNPs (proxies) for missing SNPs
- Uses 1000 Genomes phase 3 LD data
- Output: `ldlink_proxy_results_*.txt`

#### PROCESS_PROXY_RESULTS
- Filters proxies by R² > threshold (default: 0.8)
- Matches proxies to outcome GWAS SNPs
- Handles allele flipping for non-zero distance proxies
- Output: `processed_proxies.rds`

#### MERGE_PROXIED_DATA
- Combines proxies with outcome data
- Preserves original SNP information
- Prepares for outcome preprocessing
- Output: `outcome_with_proxies.rds`

## Directory Structure

Set up your pipeline directory like this:

```
mr-pipeline/
├── nextflow_mr_pipeline.nf          # Main workflow
├── nextflow.config                   # Configuration
├── README.md                         # Installation & usage
├── WORKFLOW_GUIDE.md                 # Detailed workflow explanation
├── PROXY_BATCH_INTEGRATION.md        # This file
│
└── modules/
    ├── exposure_preprocessing.nf
    ├── proxy_batch.nf                # ← NEW!
    ├── outcome_preprocessing.nf
    ├── proxy_search.nf
    ├── mr_analysis.nf
    ├── colocalization.nf
    └── results_plotting.nf
```

## Updated Pipeline Workflow

The main pipeline (`nextflow_mr_pipeline.nf`) now includes all proxy batch steps:

```nextflow
workflow {
    main:
        // Step 1: Exposure preprocessing
        exposure_preprocessed = EXPOSURE_PREPROCESSING(...)
        
        // Step 2: Find SNPs in exposure missing from outcome
        missing_snps = FIND_MISSING_SNPS(...)
        
        // Step 3: Query LDlink for proxies (BATCH MODE)
        proxy_results = PROXY_BATCH_LDLINK(...)
        
        // Step 4: Process proxy results
        proxies_processed = PROCESS_PROXY_RESULTS(...)
        
        // Step 5: Merge proxies with outcome
        proxied_outcome_data = MERGE_PROXIED_DATA(...)
        
        // Step 6: Outcome preprocessing (with proxies)
        outcome_preprocessed = OUTCOME_PREPROCESSING(...)
        
        // Steps 7-9: MR Analysis → Colocalization → Plotting
        ...
}
```

## Key Features of Proxy Batch Processing

### 1. Batch Efficiency
- Uses `LDproxy_batch()` instead of individual queries
- Massively reduces API calls to LDlink
- Faster processing for large SNP lists

### 2. Smart Filtering
```
Missing SNPs → LDlink Query → R² filtering → Outcome matching
     ↓              ↓              ↓               ↓
    100        1000+ proxies    ~200 (R²>0.8)  ~50-100 matched
```

### 3. Allele Handling
- Automatically detects when proxy ≠ original SNP position
- Flips alleles according to `Correlated_Alleles` field
- Maintains direction of effects

### 4. Progress Tracking
- `missing_snps_summary.txt` - Shows % of SNPs needing proxies
- `proxy_stats.txt` - Shows proxy matching success rate
- Detailed logging for debugging

## Usage Example

```bash
# Basic run with proxy batch processing
nextflow run nextflow_mr_pipeline.nf \
  --tissue "Whole_Blood" \
  --exposure_file "eqtl_data.txt.gz" \
  --outcome_file "gwas_data.tsv" \
  --lddlink_token "YOUR_TOKEN_HERE" \
  -profile standard

# With custom LD threshold
nextflow run nextflow_mr_pipeline.nf \
  --tissue "Heart_Left_Ventricle" \
  --exposure_file "eqtl_heart.txt.gz" \
  --outcome_file "gwas_bmi.tsv" \
  --ld_threshold 0.8 \
  --lddlink_token "YOUR_TOKEN_HERE" \
  -profile cluster
```

## Output Files from Proxy Batch

Located in `results/<TISSUE>/proxies/`:

```
proxies/
├── missing_snps_list.rds             # SNPs needing proxies
├── missing_snps_summary.txt          # Proxy statistics
│                                     
├── ldlink_proxy_results_*.txt        # Raw LDlink output (intermediate)
│   
├── processed_proxies.rds             # Filtered & processed proxies
├── proxy_stats.txt                   # Matching statistics
│
└── outcome_with_proxies.rds          # Final outcome data with proxies
```

### Example Output Statistics

```
Missing SNPs Summary
====================
Total SNPs in exposure: 50000
Total SNPs in outcome: 45000
Missing SNPs (in exposure but not outcome): 5000
Missing SNPs percent: 10.00%

Proxy Search Statistics
======================
Total LDlink proxies found: 12500
After R2 filtering (>0.8): 4200
Matching outcome GWAS: 3500
Success rate: 83.33%
```

## Important Parameters

```
--lddlink_token
  Description: LDlink API token (REQUIRED for proxy search)
  Get from: https://ldlink.nci.nih.gov/
  Example: "23a23730fd0b"

--ld_threshold
  Description: Minimum R² for proxy selection
  Default: 0.8
  Range: 0.0-1.0
  Note: Higher = more stringent (fewer proxies)

--proxy_genome_build
  Description: Genome build for proxy queries
  Default: "grch38"
  Options: "grch37", "grch38"

--proxy_population
  Description: Population for LD calculation
  Default: "ALL"
  Options: "EUR", "AFR", "AMR", "EAS", "SAS", "ALL"
```

## Handling Proxy Issues

### No Proxies Found
**Possible reasons:**
- LD threshold too high (R² > 0.8 is very stringent)
- Exposure and outcome SNP sets don't overlap well
- LDlink API temporarily down

**Solutions:**
```bash
# Reduce LD threshold
--ld_threshold 0.7

# Check token validity
--lddlink_token "your_valid_token"

# Use specific population
--proxy_population "EUR"
```

### Low Success Rate
**Symptoms:** Many proxies found but few match outcome GWAS

**Possible causes:**
- Different SNP panels between exposure and outcome
- SNP naming differences (CHR:POS vs rsid)
- Limited SNP overlap

**Solutions:**
```bash
# Pre-process to align SNP naming
# Check if SNPs in outcome data

# Lower LD threshold
--ld_threshold 0.7

# Verify input file formats match
```

### LDlink API Errors
**Symptoms:** "LDlink API timeout" or "Connection failed"

**Solutions:**
```bash
# Retry with -resume flag
nextflow run ... -resume

# Reduce batch size (modify in proxy_batch.nf if needed)

# Try different population
--proxy_population "EUR"
```

## Comparing to Your Original Scripts

### Your Original Workflow:
```R
# 1-exposure_preprocessing.R
# Process exposure, save by chromosome

# 3-proxy_batch.R
# Find missing SNPs
missing_snps <- exposure[!exposure$SNP %in% outcome$SNP]

# Use LDproxy_batch()
LDproxy_batch(missing_snps$SNP, ...)

# 4-proxy_search.R
# Process results and merge with outcome

# 5-outcome_preprocessing.R
# Preprocess outcome with merged proxies
```

### New Nextflow Pipeline:
```nextflow
// All steps integrated into a single, reproducible workflow
// Automatic parallelization
// Better error handling
// Containerized execution
// Parameter tracking & reporting
```

## Nextflow-Specific Advantages

1. **Reproducibility**: Exact same results every run
2. **Scalability**: Same code runs on laptop or HPC cluster
3. **Fault Tolerance**: Automatic retries on failures
4. **Monitoring**: Real-time execution tracking
5. **Parallelization**: Chromosome-level parallelism
6. **Reporting**: Automatic HTML reports and DAGs

## Integration Verification

To verify the proxy batch module is properly integrated:

```bash
# 1. Check pipeline syntax
nextflow run nextflow_mr_pipeline.nf -preview

# 2. View workflow DAG
nextflow run nextflow_mr_pipeline.nf -with-dag workflow.html

# 3. Dry run (no execution)
nextflow run nextflow_mr_pipeline.nf \
  --tissue "Test" \
  --exposure_file "test.txt.gz" \
  --outcome_file "test.tsv" \
  -profile standard \
  -dsl2
```

## Next Steps

1. **Place all files** in correct directory structure
2. **Get LDlink token** from https://ldlink.nci.nih.gov/
3. **Test with sample data** to verify proxy batch works
4. **Run full analysis** for your tissues of interest

## Files Checklist

✅ `nextflow_mr_pipeline.nf` - Updated with proxy batch imports
✅ `nextflow.config` - Resource configuration
✅ `README.md` - Installation & usage guide
✅ `WORKFLOW_GUIDE.md` - Detailed workflow explanation
✅ `modules/exposure_preprocessing.nf`
✅ `modules/proxy_batch.nf` - NEW!
✅ `modules/outcome_preprocessing.nf`
✅ `modules/proxy_search.nf`
✅ `modules/mr_analysis.nf`
✅ `modules/colocalization.nf`
✅ `modules/results_plotting.nf`

---

**Thank you for catching the missing proxy batch! It's now fully integrated and ready to use.**
