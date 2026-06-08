process EXPOSURE_PREPROCESSING {
    container 'rocker/tidyverse:latest'
    publishDir "${params.outdir}/${params.tissue}/exposure", mode: 'copy'
    
    input:
    path exposure_file
    val tissue
    val sample_size
    val fstat_threshold
    val outdir
    
    output:
    path "exposure_raw.rds", emit: exposure_raw
    path "exposure_by_chr_*.rds", emit: exposure_by_chr
    path "column_info.txt", emit: col_info
    
    script:
    """
    #!/usr/bin/env Rscript
    
    library(TwoSampleMR)
    library(data.table)
    library(dplyr)
    library(tidyr)
    library(stringr)
    library(reshape2)
    
    # Load exposure data - handle both gzipped and plain text
    exposure_file <- "${exposure_file}"
    cat("Reading exposure file:", exposure_file, "\\n")
    
    # Use read.delim with gzfile for gzipped files
    if (grepl("\\\\.gz\$", exposure_file)) {
        exposure <- read.delim(gzfile(exposure_file), stringsAsFactors=FALSE)
    } else {
        exposure <- read.delim(exposure_file, stringsAsFactors=FALSE)
    }
    
    cat("Loaded exposure data with", nrow(exposure), "rows\\n")
    cat("Columns:", paste(colnames(exposure), collapse=", "), "\\n\\n")
    
    # Save raw data
    saveRDS(exposure, "exposure_raw.rds")
    
    # Write column info for user reference
    writeLines(paste("Exposure file columns:", paste(colnames(exposure), collapse=", ")), "column_info.txt")
    
    # Identify key columns - more specific matching
    cols_lower <- tolower(colnames(exposure))
    
    # Gene column
    gene_col <- colnames(exposure)[cols_lower == "gene_id"][1]
    if (is.na(gene_col)) gene_col <- colnames(exposure)[grepl("gene", cols_lower)][1]
    
    # Variant column
    variant_col <- colnames(exposure)[cols_lower == "variant_id"][1]
    if (is.na(variant_col)) variant_col <- colnames(exposure)[grepl("variant|snp", cols_lower)][1]
    
    # Slope/Beta - prefer exact "slope" match
    slope_col <- colnames(exposure)[cols_lower == "slope"][1]
    if (is.na(slope_col)) slope_col <- colnames(exposure)[grepl("^beta\$|^effect\$", cols_lower)][1]
    
    # SE - prefer slope_se specifically
    se_col <- colnames(exposure)[cols_lower == "slope_se"][1]
    if (is.na(se_col)) se_col <- colnames(exposure)[grepl("se|std", cols_lower)][1]
    
    # P-value - prefer pval_nominal
    pval_col <- colnames(exposure)[cols_lower == "pval_nominal"][1]
    if (is.na(pval_col)) pval_col <- colnames(exposure)[grepl("pval|p.value|p_value", cols_lower)][1]
    
    cat(sprintf("Detected columns:\\n"))
    cat(sprintf("  Gene: %s\\n", ifelse(is.na(gene_col), "NOT FOUND", gene_col)))
    cat(sprintf("  Variant: %s\\n", ifelse(is.na(variant_col), "NOT FOUND", variant_col)))
    cat(sprintf("  Slope/Beta: %s\\n", ifelse(is.na(slope_col), "NOT FOUND", slope_col)))
    cat(sprintf("  SE: %s\\n", ifelse(is.na(se_col), "NOT FOUND", se_col)))
    cat(sprintf("  P-value: %s\\n\\n", ifelse(is.na(pval_col), "NOT FOUND", pval_col)))
    
    # Check if we have the essential columns
    if (is.na(gene_col) || is.na(variant_col) || is.na(slope_col) || is.na(se_col) || is.na(pval_col)) {
        cat("ERROR: Could not find required columns!\\n")
        cat("Looking for columns containing: gene, variant, slope, se, pval\\n")
        cat("Your file has:", paste(colnames(exposure), collapse=", "), "\\n")
        stop("Missing required columns")
    }
    
    # Rename columns for processing
    exposure <- exposure %>% 
        rename(
            gene_id = !!gene_col,
            variant_id = !!variant_col,
            slope = !!slope_col,
            slope_se = !!se_col,
            pval_nominal = !!pval_col
        )
    
    cat("Successfully renamed columns\\n")
    cat("Sample data:\\n")
    print(head(exposure[, c("gene_id", "variant_id", "slope", "slope_se", "pval_nominal")], 3))
    cat("\\n")
    
    # Calculate F-statistics for instrument strength filtering
    exposure\$R2 <- pmax(0, (get_r_from_bsen(exposure\$slope, exposure\$slope_se, ${sample_size}))^2)
    exposure\$Fstat <- (exposure\$R2/1)/((1-exposure\$R2)/(${sample_size}-2))
    
    cat(sprintf("F-statistic range: %.2f - %.2f\\n", min(exposure\$Fstat, na.rm=T), max(exposure\$Fstat, na.rm=T)))
    
    # Filter by F-statistic threshold
    exposure_filt <- exposure[exposure\$Fstat > ${fstat_threshold}, ]
    cat(sprintf("After F-stat filtering: %d SNPs (from %d)\\n", nrow(exposure_filt), nrow(exposure)))
    exposure <- exposure_filt
    
    if (nrow(exposure) == 0) {
        stop("No SNPs passed F-statistic threshold of ${fstat_threshold}")
    }
    
    # Extract gene IDs
    g_id <- colsplit(exposure\$gene_id, "\\\\.", c("genes_id", "transcript_id"))
    exposure\$genes_id <- g_id\$genes_id
    
    # Process by chromosome
    for (i in 1:22) {
        pattern <- sprintf("^chr%s_", i)
        
        # Filter for chromosome
        expo_chr <- exposure %>%
            mutate(variant_id = as.character(variant_id)) %>%
            filter(str_detect(variant_id, pattern))
        
        if (nrow(expo_chr) == 0) {
            cat(sprintf("Chromosome %s: 0 SNPs\\n", i))
            # Create empty file
            empty_df <- data.frame()
            saveRDS(empty_df, sprintf("exposure_by_chr_%s.rds", i))
            next
        }
        
        # Parse variant IDs: chr1_100_A_T format
        v_id <- colsplit(expo_chr\$variant_id, "_", 
                        c("chromosome", "base_pair_location", "other_allele", "effect_allele", "name"))
        v_id\$SNP <- paste(v_id\$chromosome, v_id\$base_pair_location, sep=":")
        
        # Construct formatted dataframe
        expo_formatted <- data.frame(
            gene_id = expo_chr\$gene_id,
            variant_id = expo_chr\$variant_id,
            SNP = v_id\$SNP,
            other_allele.exposure = v_id\$other_allele,
            effect_allele.exposure = v_id\$effect_allele,
            eaf.exposure = NA_real_,  # Set to NA if not available
            pval.exposure = expo_chr\$pval_nominal,
            beta.exposure = expo_chr\$slope,
            se.exposure = expo_chr\$slope_se
        )
        
        # Remove indels
        expo_formatted <- expo_formatted[
            expo_formatted\$other_allele.exposure %in% c('A','T','C','G') & 
            expo_formatted\$effect_allele.exposure %in% c('A','T','C','G'), ]
        
        # Remove palindromic SNPs
        expo_formatted <- expo_formatted[
            !((expo_formatted\$effect_allele.exposure == "A" & expo_formatted\$other_allele.exposure == "T") |
              (expo_formatted\$effect_allele.exposure == "T" & expo_formatted\$other_allele.exposure == "A") |
              (expo_formatted\$effect_allele.exposure == "C" & expo_formatted\$other_allele.exposure == "G") |
              (expo_formatted\$effect_allele.exposure == "G" & expo_formatted\$other_allele.exposure == "C")), ]
        
        # Add exposure metadata
        expo_formatted\$id.exposure <- expo_formatted\$gene_id
        expo_formatted\$exposure <- expo_formatted\$gene_id
        
        # Save
        saveRDS(expo_formatted, sprintf("exposure_by_chr_%s.rds", i))
        cat(sprintf("Chromosome %s: %d SNPs\\n", i, nrow(expo_formatted)))
    }
    
    cat("\\n✓ Exposure preprocessing complete!\\n")
    """
}
