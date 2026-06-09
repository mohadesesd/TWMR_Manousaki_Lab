process FIND_MISSING_SNPS {
    container 'rocker/tidyverse:latest'
    publishDir "${params.outdir}/${params.tissue}/proxies", mode: 'copy'
    
    input:
    path 'exposure_by_chr_*.rds'
    path outcome_file
    val tissue
    
    output:
    path "missing_snps_list.rds", emit: missing_snps
    path "missing_snps_summary.txt", emit: summary
    
    script:
    """
    #!/usr/bin/env Rscript
    
    library(data.table)
    library(dplyr)
    
    outcome_df <- read.delim(gzfile("${outcome_file}"), stringsAsFactors=FALSE)
    outcome_df\$SNP <- paste(outcome_df\$chromosome, ":", outcome_df\$base_pair_location, sep="")
    outcome_snps <- unique(outcome_df\$SNP)
    
    cat(sprintf("Outcome GWAS contains %d unique SNPs\\n", length(outcome_snps)))
    
    exposure_files <- list.files(".", pattern = "exposure_by_chr_.*\\\\.rds\$")
    exposure_combined <- data.frame()
    for (file in exposure_files) {
        expo <- readRDS(file)
        if (nrow(expo) > 0) exposure_combined <- rbind(exposure_combined, expo)
    }
    
    cat(sprintf("Exposure data contains %d SNPs\\n", nrow(exposure_combined)))
    
    missing_snps <- exposure_combined[!exposure_combined\$SNP %in% outcome_snps, ]
    missing_snps_unique <- unique(missing_snps\$SNP)
    
    cat(sprintf("Found %d unique SNPs in exposure missing from outcome\\n", length(missing_snps_unique)))
    
    missing_df <- data.frame(SNP = missing_snps_unique, in_exposure = TRUE, in_outcome = FALSE)
    saveRDS(missing_df, "missing_snps_list.rds")
    
    summary_text <- sprintf(
        "Missing SNPs Summary\\n====================\\nTotal SNPs in exposure: %d\\nTotal SNPs in outcome: %d\\nMissing SNPs (in exposure but not outcome): %d\\nMissing SNPs percent: %.2f%%\\n",
        nrow(exposure_combined), length(outcome_snps), nrow(missing_df),
        (nrow(missing_df) / nrow(exposure_combined) * 100))
    writeLines(summary_text, "missing_snps_summary.txt")
    cat(summary_text)
    """
}

process PROXY_BATCH_LDLINK {
    container 'rocker/tidyverse:latest'
    publishDir "${params.outdir}/${params.tissue}/proxies", mode: 'copy'
    
    maxRetries 0
    errorStrategy 'finish'
    
    input:
    path missing_snps
    val lddlink_token
    val genome_build
    val pop
    
    output:
    path "proxy_batch_log.txt", emit: log
    path "combined_query_snp_list.txt", optional: true, emit: proxy_results
    
    script:
    """
    #!/usr/bin/env Rscript
    
    library(LDlinkR)
    library(data.table)
    library(dplyr)
    
    missing_df <- readRDS("${missing_snps}")
    snps_to_query <- missing_df\$SNP
    
    cat(strrep("=", 60), "\\n")
    cat("LDLINK PROXY SEQUENTIAL QUERY (with retry & backoff)\\n")
    cat(strrep("=", 60), "\\n\\n")
    
    cat("Configuration:\\n")
    cat("  Token:", nchar("${lddlink_token}"), "characters\\n")
    cat("  Genome Build: ${genome_build}\\n")
    cat("  Population: ${pop}\\n")
    cat("  SNPs to query:", length(snps_to_query), "\\n")
    cat("  Max retries: 3\\n")
    cat("  Exponential backoff: 2s, 4s, 8s\\n\\n")
    
    log_file <- file("proxy_batch_log.txt", open = "w")
    
    if (length(snps_to_query) == 0) {
        cat("No SNPs to query!\\n")
        writeLines("No SNPs to query", log_file)
        close(log_file)
    } else {
        cat("Starting sequential LDproxy queries with retry logic...\\n\\n")
        writeLines("Starting sequential LDproxy queries with retry logic", log_file)
        
        all_proxies <- list()
        success_count <- 0
        error_count <- 0
        retry_count_total <- 0
        
        for (i in seq_along(snps_to_query)) {
            snp <- snps_to_query[i]
            if (i %% 10 == 0) cat(sprintf("Progress: %d/%d SNPs queried\\n", i, length(snps_to_query)))
            
            result <- NULL
            retry_attempt <- 0
            max_retries <- 3
            
            while (is.null(result) && retry_attempt < max_retries) {
                tryCatch({
                    result <- LDproxy(
                        snp = snp,
                        pop = "${pop}",
                        token = "${lddlink_token}",
                        genome_build = "${genome_build}",
                        r2d = "r2",
                        timeout = 30
                    )
                    
                    # Validate result
                    if (!is.data.frame(result) || nrow(result) == 0 || grepl("error", result[1,1], ignore.case=TRUE)) {
                        result <- NULL
                        retry_attempt <- retry_attempt + 1
                        if (retry_attempt < max_retries) {
                            wait_time <- 2^retry_attempt
                            cat("  Retry", retry_attempt, "for", snp, "- waiting", wait_time, "seconds...\\n")
                            writeLines(paste("  Retry", retry_attempt, "for", snp), log_file)
                            Sys.sleep(wait_time)
                        }
                    } else {
                        result\$query_snp <- snp
                        all_proxies[[i]] <- result
                        success_count <- success_count + 1
                    }
                    
                }, error = function(e) {
                    retry_attempt <<- retry_attempt + 1
                    result <<- NULL
                    
                    if (grepl("timeout|connection|killed", e\$message, ignore.case=TRUE)) {
                        if (retry_attempt < max_retries) {
                            wait_time <- 2^retry_attempt
                            cat("  [Connection error] Retry", retry_attempt, "for", snp, "- waiting", wait_time, "seconds...\\n")
                            writeLines(paste("  [Connection error] Retry", retry_attempt, "for", snp, ":", e\$message), log_file)
                            Sys.sleep(wait_time)
                        } else {
                            writeLines(paste("  [Connection failed after retries] Skipping", snp, ":", e\$message), log_file)
                        }
                    } else {
                        writeLines(paste("  Error for", snp, ":", e\$message), log_file)
                        if (retry_attempt < max_retries) {
                            wait_time <- 2^retry_attempt
                            Sys.sleep(wait_time)
                        }
                    }
                })
            }
            
            # If all retries exhausted, count as error
            if (is.null(result)) {
                error_count <- error_count + 1
                retry_count_total <- retry_count_total + retry_attempt
            }
            
            # Standard sleep between queries (increased from 0.5 to 1.5 seconds)
            Sys.sleep(1.5)
        }
        
        cat(sprintf("\\nResults: %d successful, %d errors (with %d total retry attempts)\\n", 
                   success_count, error_count, retry_count_total))
        writeLines(sprintf("Results: %d successful, %d errors (with %d total retry attempts)", 
                          success_count, error_count, retry_count_total), log_file)
        
        if (success_count > 0) {
            combined_results <- do.call(rbind, all_proxies)
            write.table(combined_results, "combined_query_snp_list.txt", sep="\\t", quote=FALSE, row.names=FALSE)
            cat(sprintf("✓ Saved %d proxy records\\n", nrow(combined_results)))
            writeLines(sprintf("✓ Saved %d proxy records", nrow(combined_results)), log_file)
        } else {
            cat("✗ No successful proxy queries\\n")
            writeLines("✗ No successful proxy queries", log_file)
        }
    }
    close(log_file)
    cat(strrep("=", 60), "\\n")
    cat("Proxy query complete\\n")
    cat(strrep("=", 60), "\\n")
    """
}

process PROCESS_PROXY_RESULTS {
    container 'rocker/tidyverse:latest'
    publishDir "${params.outdir}/${params.tissue}/proxies", mode: 'copy'
    
    input:
    path proxy_file
    path outcome_file
    val ld_threshold
    
    output:
    path "processed_proxies.rds", emit: proxies_processed
    path "proxy_stats.txt", emit: stats
    
    script:
    """
    #!/usr/bin/env Rscript
    
    library(data.table)
    library(dplyr)
    
    cat("Processing proxy results from LDlink...\\n\\n")
    proxy_file <- "${proxy_file}"
    
    if (!file.exists(proxy_file)) {
        cat("✗ Proxy file not found:", proxy_file, "\\n")
        saveRDS(data.frame(), "processed_proxies.rds")
        writeLines("Proxy file not found", "proxy_stats.txt")
    } else {
        cat("✓ Found proxy file:", proxy_file, "\\n")
        tryCatch({
            proxy_data <- fread(proxy_file, stringsAsFactors=FALSE)
            cat(sprintf("✓ Loaded %d proxy records\\n", nrow(proxy_data)))
            if (nrow(proxy_data) == 0) {
                saveRDS(data.frame(), "processed_proxies.rds")
                writeLines("No proxy data in file", "proxy_stats.txt")
            } else {
                if ("R2" %in% colnames(proxy_data)) {
                    proxy_data\$R2 <- as.numeric(proxy_data\$R2)
                    proxy_filtered <- proxy_data %>% filter(R2 > ${ld_threshold})
                    cat(sprintf("After R2 filtering (>${ld_threshold}): %d proxies\\n", nrow(proxy_filtered)))
                } else {
                    proxy_filtered <- proxy_data
                }
                saveRDS(proxy_filtered, "processed_proxies.rds")
                stats_text <- sprintf("Proxy Processing Complete\\n========================\\nTotal proxies loaded: %d\\nProxies after filtering: %d\\n",
                    nrow(proxy_data), nrow(proxy_filtered))
                writeLines(stats_text, "proxy_stats.txt")
                cat(stats_text)
            }
        }, error = function(e) {
            cat("✗ Error reading proxy file:\\n", e\$message, "\\n")
            saveRDS(data.frame(), "processed_proxies.rds")
            writeLines(paste("Error:", e\$message), "proxy_stats.txt")
        })
    }
    """
}

process MERGE_PROXIED_DATA {
    container 'rocker/tidyverse:latest'
    publishDir "${params.outdir}/${params.tissue}/proxies", mode: 'copy'
    
    input:
    path processed_proxies
    path outcome_file
    val ld_threshold
    
    output:
    path "outcome_with_proxies.rds", emit: outcome_with_proxies
    path "proxies_merged_summary.txt", emit: merge_summary
    
    script:
    """
    #!/usr/bin/env Rscript
    
    library(TwoSampleMR)
    library(data.table)
    library(dplyr)
    library(tidyr)
    library(reshape2)
    
    cat("====== PROXY MERGE WITH ALLELE FLIPPING ======\\n\\n")
    
    outcome_df <- read.delim(gzfile("${outcome_file}"), stringsAsFactors=FALSE)
    proxies <- readRDS("${processed_proxies}")
    
    # Robust sample-size detection
    if (!"TotalSampleSize" %in% colnames(outcome_df)) {
        for (alt in c("n","N","sample_size","SampleSize","n_complete_samples")) {
            if (alt %in% colnames(outcome_df)) { outcome_df\$TotalSampleSize <- outcome_df[[alt]]; break }
        }
    }
    if (!"TotalSampleSize" %in% colnames(outcome_df)) outcome_df\$TotalSampleSize <- NA_real_
    
    cat("Loaded outcome GWAS:", nrow(outcome_df), "variants\\n")
    cat("Loaded proxies:", nrow(proxies), "SNPs\\n\\n")
    
    outcome_df\$SNP <- paste(outcome_df\$chromosome, ":", outcome_df\$base_pair_location, sep="")
    outcome_df\$effect_allele <- toupper(outcome_df\$effect_allele)
    outcome_df\$other_allele <- toupper(outcome_df\$other_allele)
    outcome_df <- outcome_df[outcome_df\$other_allele %in% c('A','T','C','G') & outcome_df\$effect_allele %in% c('A','T','C','G'), ]
    
    outcome_df\$SNP_Original <- NA_character_
    outcome_df\$RS_Number_Proxy <- NA_character_
    outcome_df\$LD_R2 <- NA_real_
    outcome_df\$LD_Distance <- NA_real_
    
    cat("After filtering indels:", nrow(outcome_df), "variants\\n")
    outcome_snps <- unique(outcome_df\$SNP)
    
    outcomeList <- list()
    proxy_matches <- 0; proxy_failures <- 0; allele_flips <- 0
    
    if (nrow(proxies) > 0 && "query_snp" %in% colnames(proxies)) {
        for (idx in seq_len(nrow(proxies))) {
            proxy_row <- proxies[idx, ]
            if (!proxy_row\$Coord %in% outcome_snps) { proxy_failures <- proxy_failures + 1; next }
            outcome_row <- outcome_df[outcome_df\$SNP == proxy_row\$Coord, ]
            if (nrow(outcome_row) == 0) { proxy_failures <- proxy_failures + 1; next }
            original_snp <- proxy_row\$query_snp
            
            if (proxy_row\$Distance > 0) {
                tryCatch({
                    allele_pair <- strsplit(proxy_row\$Correlated_Alleles, ",")[[1]]
                    if (length(allele_pair) == 2) {
                        ref_alleles <- strsplit(allele_pair[1], "=")[[1]]
                        alt_alleles <- strsplit(allele_pair[2], "=")[[1]]
                        if (length(ref_alleles) == 2 && length(alt_alleles) == 2) {
                            exp_ref <- ref_alleles[1]; out_ref <- ref_alleles[2]
                            exp_alt <- alt_alleles[1]; out_alt <- alt_alleles[2]
                            if (outcome_row\$effect_allele[1] == out_alt) {
                                outcome_row\$effect_allele <- exp_alt; outcome_row\$other_allele <- exp_ref
                                allele_flips <- allele_flips + 1
                            } else if (outcome_row\$effect_allele[1] == out_ref) {
                                outcome_row\$effect_allele <- exp_ref; outcome_row\$other_allele <- exp_alt
                                allele_flips <- allele_flips + 1
                            }
                        }
                    }
                }, error = function(e) cat("Warning: Allele parsing failed for", original_snp, "\\n"))
            }
            
            # NOTE: label the proxied record with the ORIGINAL missing SNP so it
            # pairs with the exposure instrument during harmonization.
            outcome_row\$SNP <- original_snp
            outcome_row\$SNP_Original <- proxy_row\$Coord
            outcome_row\$RS_Number_Proxy <- proxy_row\$RS_Number
            outcome_row\$LD_R2 <- as.numeric(proxy_row\$R2)
            outcome_row\$LD_Distance <- proxy_row\$Distance
            
            outcomeList[[length(outcomeList) + 1]] <- outcome_row
            proxy_matches <- proxy_matches + 1
            if (proxy_matches %% 50 == 0) cat("Progress: Matched", proxy_matches, "proxies\\n")
        }
    }
    
    if (length(outcomeList) > 0) {
        outcome_proxied <- as.data.frame(data.table::rbindlist(outcomeList, fill=TRUE))
    } else {
        outcome_proxied <- outcome_df[0, ]
    }
    
    # Standardize columns expected by OUTCOME_PREPROCESSING and output PROXY ROWS ONLY
    if (nrow(outcome_proxied) > 0) {
        if (!"se" %in% colnames(outcome_proxied) && "standard_error" %in% colnames(outcome_proxied)) outcome_proxied\$se <- outcome_proxied\$standard_error
        if (!"pval" %in% colnames(outcome_proxied) && "p_value" %in% colnames(outcome_proxied)) outcome_proxied\$pval <- outcome_proxied\$p_value
        if (!"sample_size" %in% colnames(outcome_proxied)) outcome_proxied\$sample_size <- outcome_proxied\$TotalSampleSize
        saveRDS(outcome_proxied, "outcome_proxied_annotated.rds")
        keep <- c("variant_id","SNP","chromosome","base_pair_location","effect_allele","other_allele","beta","se","pval","sample_size")
        outcome_with_proxies <- outcome_proxied[, keep]
    } else {
        outcome_with_proxies <- data.frame()
    }
    saveRDS(outcome_with_proxies, "outcome_with_proxies.rds")
    
    summary_text <- sprintf(
        "PROXY MERGE SUMMARY\\n===================\\nOriginal outcome variants: %d\\nProxies matched to outcome: %d\\nProxies failed to match: %d\\nAllele flips applied: %d\\nProxy records emitted: %d\\nMean LD R2: %.4f\\n",
        nrow(outcome_df), proxy_matches, proxy_failures, allele_flips, nrow(outcome_with_proxies),
        ifelse(nrow(outcome_proxied) > 0, mean(as.numeric(outcome_proxied\$LD_R2), na.rm=TRUE), 0))
    cat(summary_text)
    writeLines(summary_text, "proxies_merged_summary.txt")
    cat("\\n✓ Merged data saved\\n")
    """
}
