process COLOCALIZATION {
    container 'rocker/tidyverse:latest'
    publishDir "${params.outdir}/${params.tissue}/colocalization", mode: 'copy'
    
    maxRetries 1
    errorStrategy 'ignore'
    
    input:
    path harmonized_files
    path mr_filtered_files
    path exposure_raw
    path outcome_raw
    val tissue
    val window_size
    val h4_threshold
    val outdir
    
    output:
    path "coloc_results_summary.rds", emit: coloc_table
    path "passing_coloc_genes.rds", emit: passing_genes
    path "coloc_genes_for_proxy.rds", emit: proxy_genes
    
    script:
    """
    #!/usr/bin/env Rscript
    
    library(coloc)
    library(susieR)
    library(gridExtra)
    library(tidyverse)
    library(LDlinkR)
    library(reshape)
    library(TwoSampleMR)
    library(data.table)
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(reshape2)
    library(ieugwasr)
    library(httr)
    library(purrr)
    library(jsonlite)
    
    # Load raw data
    exposure <- readRDS("${exposure_raw}")
    outcome <- readRDS("${outcome_raw}")
    # Robust outcome column standardization (handles n / TotalSampleSize / etc.)
    if (!"se" %in% colnames(outcome) && "standard_error" %in% colnames(outcome)) outcome\$se <- outcome\$standard_error
    if (!"pval" %in% colnames(outcome) && "p_value" %in% colnames(outcome)) outcome\$pval <- outcome\$p_value
    if (!"sample_size" %in% colnames(outcome)) {
        for (alt in c("TotalSampleSize","n","N","SampleSize","n_complete_samples")) {
            if (alt %in% colnames(outcome)) { outcome\$sample_size <- outcome[[alt]]; break }
        }
    }
    if (!"sample_size" %in% colnames(outcome)) outcome\$sample_size <- NA_real_
    
    # Load and combine harmonized data from all chromosomes
    har_files <- list.files(".", pattern = "harmonized_.*\\\\.rds\$")
    har_total <- data.frame()
    
    for (i in seq_along(har_files)) {
        df <- readRDS(har_files[i])
        if (nrow(df) > 0) {
            chr_num <- as.numeric(gsub("harmonized_chr_(.*)\\\\.rds", "\\\\1", har_files[i]))
            df\$chromosome <- sprintf("%s", chr_num)
            har_total <- rbind(har_total, df)
        }
    }
    
    # Load significant results
    mr_files <- list.files(".", pattern = "mr_filtered_chr_.*\\\\.rds\$")
    signif_total <- data.frame()
    
    for (file in mr_files) {
        df <- readRDS(file)
        if (nrow(df) > 0) {
            signif_total <- rbind(signif_total, df)
        }
    }
    
    cat(sprintf("Harmonized SNPs: %d\\n", nrow(har_total)))
    cat(sprintf("Significant genes: %d\\n", length(unique(signif_total\$id.exposure))))
    
    if (nrow(signif_total) == 0) {
        cat("No significant results found for colocalization\\n")
        
        tab_res <- data.frame(
            Gene_Name = character(),
            nSNP = numeric(),
            H0 = numeric(),
            H1 = numeric(),
            H2 = numeric(),
            H3 = numeric(),
            H4 = numeric()
        )
    } else {
        merged_df <- inner_join(har_total, signif_total, 
                               by = c("id.exposure" = "id.exposure"))
        
        tab_res <- data.frame()
        
        for (i in 1:min(length(unique(merged_df\$id.exposure)), 100)) {
            tryCatch({
                gene_id <- unique(merged_df\$id.exposure)[i]
                cat(sprintf("Processing gene %d: %s\\n", i, gene_id))
                
                entry <- merged_df[merged_df\$id.exposure == gene_id, ][1, ]
                
                if (nrow(entry) == 0) next
                
                chr <- as.character(entry\$chromosome)
                
                if ("position" %in% names(entry)) {
                    position <- entry\$position
                } else if ("base_pair_location" %in% names(entry)) {
                    position <- entry\$base_pair_location
                } else {
                    next
                }
                
                position <- as.numeric(position)
                
                pos_min <- position - ${window_size}
                pos_max <- position + ${window_size}
                
                chr_label <- paste0("chr", chr)
                
                # ---- Exposure (GTEx eGenes): filter to cis-window + format for TwoSampleMR ----
                expo_subset <- exposure %>%
                    filter(sub("_.*", "", variant_id) == chr_label,
                           variant_pos >= pos_min, variant_pos <= pos_max) %>%
                    transmute(
                        SNP = paste0(chr, ":", variant_pos),
                        effect_allele.exposure = toupper(alt),
                        other_allele.exposure  = toupper(ref),
                        beta.exposure = slope,
                        se.exposure   = slope_se,
                        pval.exposure = pval_nominal,
                        id.exposure = gene_id,
                        exposure = gene_id
                    )
                
                # ---- Outcome (GWAS): filter to cis-window + format for TwoSampleMR ----
                outcome_subset <- outcome %>%
                    filter(chromosome == chr_label,
                           base_pair_location >= pos_min, base_pair_location <= pos_max) %>%
                    transmute(
                        SNP = paste0(chromosome, ":", base_pair_location),
                        base_pair_location = base_pair_location,
                        effect_allele.outcome = toupper(effect_allele),
                        other_allele.outcome  = toupper(other_allele),
                        beta.outcome = beta,
                        se.outcome   = se,
                        pval.outcome = pval,
                        samplesize.outcome = sample_size,
                        id.outcome = "GWAS",
                        outcome = "GWAS"
                    )
                
                if (nrow(expo_subset) == 0 || nrow(outcome_subset) == 0) next
                
                # Drop palindromic SNPs + set dummy eaf (same trick that fixed the MR step;
                # prevents the harmonise NA-eaf crash, and coloc here uses sdY not MAF)
                is_pal <- function(a1, a2) (a1=="A"&a2=="T")|(a1=="T"&a2=="A")|(a1=="C"&a2=="G")|(a1=="G"&a2=="C")
                expo_subset <- expo_subset[!is_pal(expo_subset\$effect_allele.exposure, expo_subset\$other_allele.exposure), ]
                expo_subset\$eaf.exposure <- 0.3
                outcome_subset\$eaf.outcome <- 0.3
                
                harm_data <- tryCatch({
                    harmonise_data(expo_subset, outcome_subset, action = 2)
                }, error = function(e) { cat(sprintf("  harmonise err: %s\\n", e\$message)); data.frame() })
                
                if (nrow(harm_data) == 0) next
                
                harm_data <- harm_data[!duplicated(harm_data\$SNP), ]
                
                if (nrow(harm_data) < 3) next
                
                D1_expo <- list(
                    type = "quant",
                    snp = harm_data\$SNP,
                    position = as.numeric(harm_data\$base_pair_location),
                    beta = as.numeric(harm_data\$beta.exposure),
                    varbeta = as.numeric(harm_data\$se.exposure)^2,
                    sdY = rep(1, nrow(harm_data))
                )
                
                D2_outc <- list(
                    type = "quant",
                    snp = harm_data\$SNP,
                    position = as.numeric(harm_data\$base_pair_location),
                    beta = as.numeric(harm_data\$beta.outcome),
                    varbeta = as.numeric(harm_data\$se.outcome)^2,
                    sdY = rep(1, nrow(harm_data))
                )
                
                res <- coloc.abf(dataset1 = D1_expo, dataset2 = D2_outc)
                
                gene_name <- ifelse("geneName" %in% names(entry), entry\$geneName, gene_id)
                tab_res <- rbind(tab_res, c(as.character(gene_name), 
                                           as.vector(res\$summary)))
                
            }, error = function(e) {
                cat(sprintf("Error processing gene %d: %s\\n", i, e\$message))
            })
        }
    }
    
    if (nrow(tab_res) > 0) {
        colnames(tab_res) <- c("Gene_Name", "nSNP", "H0", "H1", "H2", "H3", "H4")
        tab_res\$H0 <- as.numeric(tab_res\$H0)
        tab_res\$H1 <- as.numeric(tab_res\$H1)
        tab_res\$H2 <- as.numeric(tab_res\$H2)
        tab_res\$H3 <- as.numeric(tab_res\$H3)
        tab_res\$H4 <- as.numeric(tab_res\$H4)
    }
    
    saveRDS(tab_res, "coloc_results_summary.rds")
    
    passing <- tab_res[tab_res\$H4 > ${h4_threshold}, ]
    failing <- tab_res[tab_res\$H4 <= ${h4_threshold}, ]
    
    saveRDS(passing, "passing_coloc_genes.rds")
    saveRDS(failing, "coloc_genes_for_proxy.rds")
    
    cat(sprintf("Colocalization complete: %d genes tested, %d passing (H4 > %s)\\n", 
                nrow(tab_res), nrow(passing), ${h4_threshold}))
    """
}
