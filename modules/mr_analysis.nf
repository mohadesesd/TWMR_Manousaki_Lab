process MR_ANALYSIS {
    container 'rocker/tidyverse:latest'
    publishDir "${params.outdir}/${params.tissue}/mr_results", mode: 'copy'
    
    input:
    tuple val(chr), path(exposure_by_chr), path(outcome_by_chr)
    val tissue
    val p_threshold
    val outdir
    
    output:
    path "harmonized_chr_${chr}.rds", emit: harmonized_data
    path "mr_raw_chr_${chr}.rds", emit: all_results
    path "mr_filtered_chr_${chr}.rds", emit: filtered_results
    
    script:
    """
    #!/usr/bin/env Rscript
    
    suppressMessages({
        library(TwoSampleMR)
        library(data.table)
        library(dplyr)
    })
    
    exposure <- readRDS("${exposure_by_chr}")
    outcome <- readRDS("${outcome_by_chr}")
    
    cat(sprintf("Chr ${chr}: Exposure = %d, Outcome = %d\\n", 
                nrow(exposure), nrow(outcome)))
    
    # Set gene_id as the exposure identifier (each gene = separate exposure)
    if ("gene_id" %in% colnames(exposure)) {
        exposure\$id.exposure <- exposure\$gene_id
        exposure\$exposure <- exposure\$gene_id
    }
    
    # Remove duplicate SNPs
    exposure <- exposure %>% distinct(id.exposure, SNP, .keep_all = TRUE)
    outcome <- outcome %>% distinct(SNP, .keep_all = TRUE)
    
    # Remove palindromic SNPs (A/T, C/G) - cannot strand-resolve without allele frequencies
    is_palindromic <- function(a1, a2) {
        (a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") |
        (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C")
    }
    n_before <- nrow(exposure)
    exposure <- exposure %>% 
        filter(!is_palindromic(effect_allele.exposure, other_allele.exposure))
    cat(sprintf("Chr ${chr}: Removed %d palindromic SNPs (%d remain)\\n", 
                n_before - nrow(exposure), nrow(exposure)))
    
    # FIX: Add dummy eaf values to prevent harmonise_data crash with NA eaf
    # (Safe because palindromic SNPs already removed - eaf not needed for strand resolution)
    exposure\$eaf.exposure <- 0.3
    outcome\$eaf.outcome <- 0.3
    
    cat(sprintf("Chr ${chr}: Genes = %d\\n", length(unique(exposure\$id.exposure))))
    
    # Harmonize
    harmonised <- tryCatch({
        harmonise_data(exposure_dat = exposure, outcome_dat = outcome, action = 2)
    }, error = function(e) {
        cat(sprintf("Harmonize error chr ${chr}: %s\\n", e\$message))
        return(data.frame())
    })
    
    if (is.null(harmonised) || nrow(harmonised) == 0) {
        cat(sprintf("Chr ${chr}: No harmonized SNPs\\n"))
        saveRDS(data.frame(), "harmonized_chr_${chr}.rds")
        saveRDS(data.frame(), "mr_raw_chr_${chr}.rds")
        saveRDS(data.frame(), "mr_filtered_chr_${chr}.rds")
    } else {
        cat(sprintf("Chr ${chr}: Harmonized = %d SNPs, %d genes\\n", 
                    nrow(harmonised), length(unique(harmonised\$id.exposure))))
        saveRDS(harmonised, "harmonized_chr_${chr}.rds")
        
        # Run MR
        mr_results <- tryCatch({
            mr(harmonised)
        }, error = function(e) {
            cat(sprintf("MR error chr ${chr}: %s\\n", e\$message))
            return(data.frame())
        })
        
        if (!is.null(mr_results) && nrow(mr_results) > 0) {
            cat(sprintf("Chr ${chr}: MR results = %d\\n", nrow(mr_results)))
            saveRDS(mr_results, "mr_raw_chr_${chr}.rds")
            
            unique_genes <- length(unique(harmonised\$id.exposure))
            p_cutoff <- ${p_threshold} / unique_genes
            mr_significant <- mr_results %>% filter(pval < p_cutoff)
            
            cat(sprintf("Chr ${chr}: Significant (p < %.2e) = %d\\n", 
                        p_cutoff, nrow(mr_significant)))
            saveRDS(mr_significant, "mr_filtered_chr_${chr}.rds")
        } else {
            cat(sprintf("Chr ${chr}: No MR results\\n"))
            saveRDS(data.frame(), "mr_raw_chr_${chr}.rds")
            saveRDS(data.frame(), "mr_filtered_chr_${chr}.rds")
        }
    }
    """
}

process COMBINE_MR_RESULTS {
    container 'rocker/tidyverse:latest'
    publishDir "${params.outdir}/${params.tissue}/mr_results", mode: 'copy'
    
    input:
    path 'mr_results_*.rds'
    val tissue
    
    output:
    path "all_mr_results_combined.rds", emit: combined
    path "significant_results_combined.rds", emit: significant
    
    script:
    """
    #!/usr/bin/env Rscript
    
    suppressMessages({ library(dplyr); library(data.table) })
    
    mr_files <- list.files(".", pattern = "mr_raw_chr_.*\\\\.rds\$")
    mr_filtered_files <- list.files(".", pattern = "mr_filtered_chr_.*\\\\.rds\$")
    
    all_results <- data.frame()
    significant_results <- data.frame()
    
    for (file in mr_files) {
        df <- readRDS(file)
        if (nrow(df) > 0) all_results <- rbind(all_results, df)
    }
    for (file in mr_filtered_files) {
        df <- readRDS(file)
        if (nrow(df) > 0) significant_results <- rbind(significant_results, df)
    }
    
    cat(sprintf("Combined: Total = %d, Significant = %d\\n", 
                nrow(all_results), nrow(significant_results)))
    saveRDS(all_results, "all_mr_results_combined.rds")
    saveRDS(significant_results, "significant_results_combined.rds")
    """
}
