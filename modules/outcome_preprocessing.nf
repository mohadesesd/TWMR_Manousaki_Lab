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
    
    if ("${chr}" == "1") {
        saveRDS(outcome_df, "outcome_raw.rds")
    }
    
    outcome_df <- outcome_df %>%
        mutate(
            se = if ("standard_error" %in% colnames(.)) standard_error else se,
            pval = if ("p_value" %in% colnames(.)) p_value else pval,
            sample_size = if ("TotalSampleSize" %in% colnames(.)) TotalSampleSize else sample_size
        )
    
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
