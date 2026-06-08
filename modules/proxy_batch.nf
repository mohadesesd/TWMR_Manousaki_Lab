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
        if (nrow(expo) > 0) {
            exposure_combined <- rbind(exposure_combined, expo)
        }
    }
    
    cat(sprintf("Exposure data contains %d SNPs\\n", nrow(exposure_combined)))
    
    missing_snps <- exposure_combined[!exposure_combined\$SNP %in% outcome_snps, ]
    missing_snps_unique <- unique(missing_snps\$SNP)
    
    cat(sprintf("Found %d unique SNPs in exposure missing from outcome\\n", length(missing_snps_unique)))
    
    missing_df <- data.frame(
        SNP = missing_snps_unique,
        in_exposure = TRUE,
        in_outcome = FALSE
    )
    
    saveRDS(missing_df, "missing_snps_list.rds")
    
    summary_text <- sprintf(
        "Missing SNPs Summary\\n====================\\nTotal SNPs in exposure: %d\\nTotal SNPs in outcome: %d\\nMissing SNPs (in exposure but not outcome): %d\\nMissing SNPs percent: %.2f%%\\n",
        nrow(exposure_combined),
        length(outcome_snps),
        nrow(missing_df),
        (nrow(missing_df) / nrow(exposure_combined) * 100)
    )
    
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
    cat("LDLINK PROXY SEQUENTIAL QUERY\\n")
    cat(strrep("=", 60), "\\n\\n")
    
    cat("Configuration:\\n")
    cat("  Token:", nchar("${lddlink_token}"), "characters\\n")
    cat("  Genome Build: ${genome_build}\\n")
    cat("  Population: ${pop}\\n")
    cat("  SNPs to query:", length(snps_to_query), "\\n\\n")
    
    log_file <- file("proxy_batch_log.txt", open = "w")
    
    if (length(snps_to_query) == 0) {
        cat("No SNPs to query!\\n")
        writeLines("No SNPs to query", log_file)
        close(log_file)
    } else {
        cat("Starting sequential LDproxy queries...\\n\\n")
        
        writeLines("Starting sequential LDproxy queries", log_file)
        
        all_proxies <- list()
        success_count <- 0
        error_count <- 0
        
        for (i in seq_along(snps_to_query)) {
            snp <- snps_to_query[i]
            
            if (i %% 10 == 0) {
                cat(sprintf("Progress: %d/%d SNPs queried\\n", i, length(snps_to_query)))
            }
            
            tryCatch({
                result <- LDproxy(
                    snp = snp,
                    pop = "${pop}",
                    token = "${lddlink_token}",
                    r2d = "r2"
                )
                
                if (is.data.frame(result) && nrow(result) > 0 && !grepl("error", result[1,1], ignore.case=TRUE)) {
                    all_proxies[[i]] <- result
                    success_count <- success_count + 1
                } else {
                    error_count <- error_count + 1
                }
                
                Sys.sleep(0.5)
                
            }, error = function(e) {
                error_count <<- error_count + 1
                writeLines(paste("Error for", snp, ":", e\$message), log_file)
            })
        }
        
        cat(sprintf("\\nResults: %d successful, %d errors\\n", success_count, error_count))
        writeLines(sprintf("Results: %d successful, %d errors", success_count, error_count), log_file)
        
        if (success_count > 0) {
            combined_results <- do.call(rbind, all_proxies)
            
            write.table(combined_results, "combined_query_snp_list.txt", 
                       sep="\\t", quote=FALSE, row.names=FALSE)
            
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
                cat("✗ No proxy data in file\\n")
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
                
                stats_text <- sprintf(
                    "Proxy Processing Complete\\n========================\\nTotal proxies loaded: %d\\nProxies after filtering: %d\\n",
                    nrow(proxy_data),
                    nrow(proxy_filtered)
                )
                writeLines(stats_text, "proxy_stats.txt")
                cat(stats_text)
                cat("✓ Processed proxies saved\\n")
            }
        }, error = function(e) {
            cat("✗ Error reading proxy file:\\n")
            cat(e\$message, "\\n\\n")
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
    
    # Load data
    outcome_df <- read.delim(gzfile("${outcome_file}"), stringsAsFactors=FALSE)
    proxies <- readRDS("${processed_proxies}")
    
    cat("Loaded outcome GWAS:", nrow(outcome_df), "variants\\n")
    cat("Loaded proxies:", nrow(proxies), "SNPs\\n\\n")
    
    # Prepare outcome data
    outcome_df\$SNP <- paste(outcome_df\$chromosome, ":", outcome_df\$base_pair_location, sep="")
    outcome_df\$effect_allele <- toupper(outcome_df\$effect_allele)
    outcome_df\$other_allele <- toupper(outcome_df\$other_allele)
    
    # Filter for single nucleotide variants (remove indels)
    outcome_df <- outcome_df[
        outcome_df\$other_allele %in% c('A','T','C','G') & 
        outcome_df\$effect_allele %in% c('A','T','C','G'), 
    ]
    
    cat("After filtering indels:", nrow(outcome_df), "variants\\n")
    outcome_snps <- unique(outcome_df\$SNP)
    cat("Unique outcome SNPs:", length(outcome_snps), "\\n\\n")
    
    # Initialize output
    outcomeList <- list()
    proxy_matches <- 0
    proxy_failures <- 0
    allele_flips <- 0
    
    cat("Starting proxy search with allele handling...\\n")
    cat(strrep("-", 60), "\\n")
    
    # Process each proxy
    if (nrow(proxies) > 0) {
        for (idx in seq_len(nrow(proxies))) {
            proxy_row <- proxies[idx, ]
            
            # Check if proxy SNP exists in outcome
            if (!proxy_row\$Coord %in% outcome_snps) {
                proxy_failures <- proxy_failures + 1
                next
            }
            
            # Get outcome row for this proxy coordinate
            outcome_row <- outcome_df[outcome_df\$SNP == proxy_row\$Coord, ]
            
            if (nrow(outcome_row) == 0) {
                proxy_failures <- proxy_failures + 1
                next
            }
            
            # Extract the original SNP being proxied
            original_snp <- proxy_row\$query_snp
            
            # Handle allele flipping if distance > 0
            if (proxy_row\$Distance > 0) {
                tryCatch({
                    # Parse correlated alleles (format: "A=A,T=C")
                    allele_pair <- strsplit(proxy_row\$Correlated_Alleles, ",")[[1]]
                    
                    if (length(allele_pair) == 2) {
                        ref_alleles <- strsplit(allele_pair[1], "=")[[1]]
                        alt_alleles <- strsplit(allele_pair[2], "=")[[1]]
                        
                        if (length(ref_alleles) == 2 && length(alt_alleles) == 2) {
                            exp_ref <- ref_alleles[1]
                            out_ref <- ref_alleles[2]
                            exp_alt <- alt_alleles[1]
                            out_alt <- alt_alleles[2]
                            
                            # Check allele match and flip if needed
                            if (outcome_row\$effect_allele[1] == out_alt) {
                                outcome_row\$effect_allele <- exp_alt
                                outcome_row\$other_allele <- exp_ref
                                allele_flips <- allele_flips + 1
                            } else if (outcome_row\$effect_allele[1] == out_ref) {
                                outcome_row\$effect_allele <- exp_ref
                                outcome_row\$other_allele <- exp_alt
                                allele_flips <- allele_flips + 1
                            }
                        }
                    }
                }, error = function(e) {
                    cat("Warning: Allele parsing failed for", original_snp, "\\n")
                })
            }
            
            # Update SNP identifiers
            outcome_row\$SNP <- proxy_row\$Coord
            outcome_row\$SNP_Original <- original_snp
            outcome_row\$RS_Number_Proxy <- proxy_row\$RS_Number
            outcome_row\$LD_R2 <- as.numeric(proxy_row\$R2)
            outcome_row\$LD_Distance <- proxy_row\$Distance
            
            outcomeList[[length(outcomeList) + 1]] <- outcome_row
            proxy_matches <- proxy_matches + 1
            
            if (proxy_matches %% 50 == 0) {
                cat("Progress: Matched", proxy_matches, "proxies\\n")
            }
        }
    }
    
    cat(strrep("-", 60), "\\n\\n")
    
    # Bind all matched proxies
    if (length(outcomeList) > 0) {
        outcome_proxied <- data.table::rbindlist(outcomeList, fill=TRUE)
    } else {
        outcome_proxied <- outcome_df[0, ]
    }
    
    # Combine original outcome + proxied outcomes
    outcome_with_proxies <- rbind(outcome_df, outcome_proxied, fill=TRUE)
    outcome_with_proxies <- as.data.frame(outcome_with_proxies)
    
    # Save results
    saveRDS(outcome_with_proxies, "outcome_with_proxies.rds")
    
    # Generate summary
    summary_text <- sprintf(
        "PROXY MERGE SUMMARY\\n===================\\nOriginal outcome variants: %d\\nProxies matched to outcome: %d\\nProxies failed to match: %d\\nAllele flips applied: %d\\nFinal merged variants: %d\\n\\nStatistics:\\n- Mean LD R²: %.4f\\n- Proxies with allele flipping: %.1f%%\\n",
        nrow(outcome_df),
        proxy_matches,
        proxy_failures,
        allele_flips,
        nrow(outcome_with_proxies),
        ifelse(nrow(outcome_proxied) > 0, mean(outcome_proxied\$LD_R2, na.rm=TRUE), 0),
        ifelse(proxy_matches > 0, (allele_flips / proxy_matches * 100), 0)
    )
    
    cat(summary_text)
    writeLines(summary_text, "proxies_merged_summary.txt")
    cat("\\n✓ Merged data saved to outcome_with_proxies.rds\\n")
    """
}
