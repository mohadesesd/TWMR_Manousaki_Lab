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
        "Missing SNPs Summary\\n
        ====================\\n
        Total SNPs in exposure: %d\\n
        Total SNPs in outcome: %d\\n
        Missing SNPs (in exposure but not outcome): %d\\n
        Missing SNPs percent: %.2f%%\\n",
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
                    "Proxy Processing Complete\\n
                    ========================\\n
                    Total proxies loaded: %d\\n
                    Proxies after filtering: %d\\n",
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
    
    output:
    path "outcome_with_proxies.rds", emit: outcome_with_proxies
    
    script:
    """
    #!/usr/bin/env Rscript
    
    library(data.table)
    library(dplyr)
    
    cat("Merging proxies with outcome GWAS...\\n\\n")
    
    proxies <- readRDS("${processed_proxies}")
    outcome_df <- read.delim(gzfile("${outcome_file}"), stringsAsFactors=FALSE)
    
    cat("Original outcome rows:", nrow(outcome_df), "\\n")
    
    if (nrow(proxies) == 0) {
        cat("No proxies found - using outcome GWAS only\\n")
        outcome_with_proxies <- outcome_df
    } else {
        cat("Proxies loaded:", nrow(proxies), "\\n\\n")
        
        outcome_formatted <- outcome_df %>%
            mutate(
                SNP = paste0(chromosome, ":", base_pair_location),
                effect_allele = tolower(effect_allele),
                other_allele = tolower(other_allele),
                pval = p_value,
                se = standard_error,
                sample_size = TotalSampleSize
            ) %>%
            select(variant_id, SNP, chromosome, base_pair_location,
                   effect_allele, other_allele, beta, se, pval, sample_size)
        
        proxies_formatted <- proxies %>%
            mutate(
                allele1 = substr(Alleles, 2, 2),
                allele2 = substr(Alleles, 4, 4),
                effect_allele = tolower(allele1),
                other_allele = tolower(allele2),
                SNP = Coord,
                variant_id = RS_Number,
                chromosome = sub(":.*", "", Coord),
                base_pair_location = as.numeric(sub(".*:", "", Coord)),
                beta = NA_real_,
                se = NA_real_,
                pval = NA_real_,
                sample_size = NA_real_
            ) %>%
            select(variant_id, SNP, chromosome, base_pair_location,
                   effect_allele, other_allele, beta, se, pval, sample_size)
        
        outcome_with_proxies <- rbind(outcome_formatted, proxies_formatted)
        
        cat("===== MERGE SUMMARY =====\\n")
        cat("Outcome GWAS: ", nrow(outcome_formatted), " variants\\n")
        cat("Added proxies: ", nrow(proxies_formatted), " variants\\n")
        cat("Total: ", nrow(outcome_with_proxies), " variants\\n")
        cat("========================\\n\\n")
    }
    
    saveRDS(outcome_with_proxies, "outcome_with_proxies.rds")
    cat("✓ Merged data saved\\n")
    """
}
