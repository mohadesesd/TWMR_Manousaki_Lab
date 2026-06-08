process OUTCOME_PREPROCESSING {
    container 'rocker/tidyverse:latest'
    publishDir "${params.outdir}/${params.tissue}/outcome", mode: 'copy'
    
    input:
    path outcome_file
    path proxied_data
    val tissue
    each chr
    val outdir
    
    output:
    path "outcome_chr_${chr}.rds", emit: outcome_by_chr
    path "outcome_raw.rds", optional: true, emit: outcome_raw
    
    script:
    """
    #!/usr/bin/env Rscript
    
    library(data.table)
    library(dplyr)
    library(tidyr)
    
    cat("Loading outcome GWAS...\\n")
    outcome_df <- as.data.frame(fread("${outcome_file}"))
    cn <- colnames(outcome_df)
    cat("Outcome columns:", paste(cn, collapse=", "), "\\n")
    
    # ---- Robust column standardization (works across GWAS formats) ----
    # standard error -> se
    if (!"se" %in% cn) {
        if ("standard_error" %in% cn) outcome_df\$se <- outcome_df\$standard_error
        else if ("SE" %in% cn)        outcome_df\$se <- outcome_df\$SE
    }
    # p-value -> pval
    if (!"pval" %in% cn) {
        if ("p_value" %in% cn)   outcome_df\$pval <- outcome_df\$p_value
        else if ("P" %in% cn)    outcome_df\$pval <- outcome_df\$P
        else if ("pvalue" %in% cn) outcome_df\$pval <- outcome_df\$pvalue
    }
    # sample size -> sample_size  (auto-detect among common names)
    if (!"sample_size" %in% colnames(outcome_df)) {
        for (alt in c("TotalSampleSize","n","N","SampleSize","n_complete_samples","Neff")) {
            if (alt %in% cn) { outcome_df\$sample_size <- outcome_df[[alt]]; break }
        }
    }
    if (!"sample_size" %in% colnames(outcome_df)) {
        cat("WARNING: no sample-size column found; setting NA\\n")
        outcome_df\$sample_size <- NA_real_
    }
    
    if ("${chr}" == "1") {
        saveRDS(outcome_df, "outcome_raw.rds")
    }
    
    chr_pattern <- sprintf("chr%s", ${chr})
    outcome_chr <- outcome_df %>% 
        filter(chromosome == chr_pattern) %>%
        mutate(
            effect_allele = toupper(effect_allele),
            other_allele = toupper(other_allele),
            SNP = paste(chromosome, ":", base_pair_location, sep="")
        ) %>%
        select(variant_id, SNP, chromosome, base_pair_location,
               effect_allele, other_allele, beta, se, pval, sample_size)
    
    cat(sprintf("Chromosome ${chr}: %d variants\\n", nrow(outcome_chr)))
    
    proxied <- tryCatch({
        prox <- readRDS("${proxied_data}")
        prox %>% 
            filter(chromosome == chr_pattern) %>%
            select(variant_id, SNP, chromosome, base_pair_location,
                   effect_allele, other_allele, beta, se, pval, sample_size)
    }, error = function(e) {
        data.frame()
    })
    
    if (nrow(proxied) > 0) {
        outcome_chr <- rbind(outcome_chr, proxied)
    }
    
    outcome_chr <- outcome_chr[
        outcome_chr\$other_allele %in% c('A','T','C','G') & 
        outcome_chr\$effect_allele %in% c('A','T','C','G'), ]
    
    outcome_formatted <- outcome_chr %>%
        rename(
            other_allele.outcome = other_allele,
            effect_allele.outcome = effect_allele,
            beta.outcome = beta,
            se.outcome = se,
            pval.outcome = pval,
            samplesize.outcome = sample_size
        ) %>%
        select(variant_id, SNP, chromosome, base_pair_location,
               other_allele.outcome, effect_allele.outcome, 
               beta.outcome, se.outcome, pval.outcome, samplesize.outcome)
    
    outcome_formatted <- outcome_formatted[
        !((outcome_formatted\$effect_allele.outcome == "A" & outcome_formatted\$other_allele.outcome == "T") |
          (outcome_formatted\$effect_allele.outcome == "T" & outcome_formatted\$other_allele.outcome == "A") |
          (outcome_formatted\$effect_allele.outcome == "C" & outcome_formatted\$other_allele.outcome == "G") |
          (outcome_formatted\$effect_allele.outcome == "G" & outcome_formatted\$other_allele.outcome == "C")), ]
    
    outcome_formatted\$outcome <- "GWAS"
    outcome_formatted\$id.outcome <- "GWAS"
    
    saveRDS(outcome_formatted, "outcome_chr_${chr}.rds")
    cat(sprintf("✓ Saved chromosome ${chr}\\n"))
    """
}
