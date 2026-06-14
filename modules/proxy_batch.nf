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
    
    missing_rows <- exposure_combined[!exposure_combined\$SNP %in% outcome_snps, ]
    missing_rows <- missing_rows[!duplicated(missing_rows\$SNP), ]
    
    cat(sprintf("Found %d unique SNPs in exposure missing from outcome\\n", nrow(missing_rows)))
    
    # Detect exposure allele columns (effect allele = ALT in eGenes)
    ea_col <- intersect(c("effect_allele.exposure","effect_allele","alt","ALT","Alt"), colnames(missing_rows))[1]
    oa_col <- intersect(c("other_allele.exposure","other_allele","ref","REF","Ref"), colnames(missing_rows))[1]
    cat("Exposure effect-allele column:", ifelse(is.na(ea_col),"<none>",ea_col),
        "| other-allele column:", ifelse(is.na(oa_col),"<none>",oa_col), "\\n")
    
    missing_df <- data.frame(
        SNP               = missing_rows\$SNP,
        effect_allele_exp = if (!is.na(ea_col)) toupper(as.character(missing_rows[[ea_col]])) else NA_character_,
        other_allele_exp  = if (!is.na(oa_col)) toupper(as.character(missing_rows[[oa_col]])) else NA_character_,
        in_exposure = TRUE, in_outcome = FALSE,
        stringsAsFactors = FALSE
    )
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
    path "combined_query_snp_list_${genome_build}.txt", optional: true, emit: proxy_results

    script:
    """
    #!/usr/bin/env Rscript

    library(LDlinkR)
    library(data.table)

    missing_df <- readRDS("${missing_snps}")
    snps_to_query <- unique(missing_df\$SNP)

    cat(strrep("=", 60), "\\n")
    cat("LDLINK PROXY BATCH QUERY (all SNPs in one call)\\n")
    cat(strrep("=", 60), "\\n\\n")
    cat("  Genome Build: ${genome_build}\\n")
    cat("  Population: ${pop}\\n")
    cat("  SNPs to query:", length(snps_to_query), "\\n\\n")

    log_con <- file("proxy_batch_log.txt", open = "w")
    writeLines(sprintf("LDproxy_batch on %d SNPs (build=${genome_build}, pop=${pop})",
                       length(snps_to_query)), log_con)

    if (length(snps_to_query) == 0) {
        cat("No SNPs to query!\\n")
        writeLines("No SNPs to query", log_con)
        close(log_con)
    } else {
        # LDproxy_batch with append=TRUE writes a single combined file:
        #   combined_query_snp_list.txt   (includes a query_snp column)
        ok <- tryCatch({
            LDproxy_batch(
                snp          = snps_to_query,
                pop          = "${pop}",
                r2d          = "r2",
                token        = "${lddlink_token}",
                append       = TRUE,
                genome_build = "${genome_build}"
            )
            TRUE
        }, error = function(e) {
            cat("LDproxy_batch error:", e\$message, "\\n")
            writeLines(paste("LDproxy_batch error:", e\$message), log_con)
            FALSE
        })

        out_file <- paste0("combined_query_snp_list_", "${genome_build}", ".txt")
        if (file.exists(out_file)) {
            n <- tryCatch(nrow(fread(out_file)), error = function(e) NA)
            cat(sprintf("✓ %s written (%s rows)\\n", out_file, n))
            writeLines(sprintf("Saved %s (%s rows)", out_file, n), log_con)
        } else {
            cat("✗ No", out_file, "produced\\n")
            writeLines(paste("No", out_file, "produced"), log_con)
        }
        close(log_con)
    }

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
    path missing_snps
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
    missing_info <- readRDS("${missing_snps}")

    # Exposure effect/other allele lookup keyed by the missing SNP (query_snp)
    exp_ea <- setNames(toupper(as.character(missing_info\$effect_allele_exp)), missing_info\$SNP)
    exp_oa <- setNames(toupper(as.character(missing_info\$other_allele_exp)),  missing_info\$SNP)
    
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
            
            # ---- Exposure-relative allele alignment + beta sign flip ----
            # Label the proxied record with the EXPOSURE effect allele (eGenes ALT);
            # flip beta when that allele pairs (in LD) with the proxy's OTHER allele.
            ea_exp   <- exp_ea[[original_snp]]
            ea_out   <- toupper(outcome_row\$effect_allele[1])   # proxy GWAS effect allele
            oa_out   <- toupper(outcome_row\$other_allele[1])    # proxy GWAS other allele
            beta_use <- as.numeric(outcome_row\$beta[1])
            eff_final <- NA_character_; oth_final <- NA_character_
            aligned <- FALSE

            if (proxy_row\$Distance == 0) {
                # same locus as the (missing) query SNP: keep GWAS alleles/beta
                eff_final <- ea_out; oth_final <- oa_out
                aligned   <- !is.na(beta_use)
            } else if (!is.na(ea_exp)) {
                ok <- tryCatch({
                    cors <- strsplit(proxy_row\$Correlated_Alleles, ",")[[1]]
                    res <- FALSE
                    if (length(cors) == 2) {
                        p1 <- strsplit(cors[1], "=")[[1]]   # c(exp_ref, out_ref)
                        p2 <- strsplit(cors[2], "=")[[1]]   # c(exp_alt, out_alt)
                        if (length(p1) == 2 && length(p2) == 2) {
                            exp_ref <- toupper(p1[1]); out_ref <- toupper(p1[2])
                            exp_alt <- toupper(p2[1]); out_alt <- toupper(p2[2])
                            if (ea_exp %in% c(exp_ref, exp_alt)) {
                                # exposure other allele + the proxy allele paired with ea_exp
                                if (ea_exp == exp_ref) { oth_final <- exp_alt; corr_proxy <- out_ref }
                                else                   { oth_final <- exp_ref; corr_proxy <- out_alt }
                                eff_final <- ea_exp
                                if (corr_proxy == ea_out) {
                                    res <- TRUE                      # pairs with proxy EFFECT allele -> beta as-is
                                } else if (corr_proxy == oa_out) {
                                    beta_use <- -beta_use            # pairs with proxy OTHER allele -> FLIP beta
                                    allele_flips <- allele_flips + 1
                                    res <- TRUE
                                }
                                # else: GWAS alleles disagree w/ Correlated_Alleles (strand) -> res stays FALSE
                            }
                        }
                    }
                    res
                }, error = function(e) { cat("Allele parse failed for", original_snp, "\\n"); FALSE })
                aligned <- isTRUE(ok)
            }

            if (!aligned || is.na(eff_final) || is.na(oth_final) || is.na(beta_use)) {
                proxy_failures <- proxy_failures + 1
                next
            }

            outcome_row\$effect_allele   <- eff_final
            outcome_row\$other_allele    <- oth_final
            outcome_row\$beta            <- beta_use
            outcome_row\$SNP             <- original_snp
            outcome_row\$SNP_Original    <- proxy_row\$Coord
            outcome_row\$RS_Number_Proxy <- proxy_row\$RS_Number
            outcome_row\$LD_R2           <- as.numeric(proxy_row\$R2)
            outcome_row\$LD_Distance     <- proxy_row\$Distance
            
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
