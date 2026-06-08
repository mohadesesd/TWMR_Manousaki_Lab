process RESULTS_PLOTTING {
    container 'rocker/tidyverse:latest'
    publishDir "${params.outdir}/${params.tissue}/plots", mode: 'copy'
    
    input:
    path 'mr_results_*.rds'
    path coloc_results
    path gene_conversion_table
    val tissue
    val outdir
    
    output:
    path "*.pdf", emit: plots
    path "*.png", emit: png_plots
    path "*_summary_stats.txt", emit: summary_stats
    path "significant_results_with_names.rds", optional: true, emit: named_results
    
    script:
    """
    #!/usr/bin/env Rscript
    
    library(ggplot2)
    library(dplyr)
    library(gridExtra)
    library(data.table)
    library(tidyverse)
    
    # Load all MR results
    mr_files <- list.files(".", pattern = "mr_results_.*\\\\.rds\$")
    mr_combined <- data.frame()
    for (file in mr_files) {
        df <- readRDS(file)
        if (nrow(df) > 0) mr_combined <- rbind(mr_combined, df)
    }
    
    # Load colocalization results
    coloc <- readRDS("${coloc_results}")
    
    # Load gene name conversion table (columns: id, name)
    gtf_reduced <- readRDS("${gene_conversion_table}")
    
    # Helper: map id.exposure (e.g. ENSG00000172260.15) -> gene symbol
    add_gene_names <- function(df) {
        ens <- sub(".*?(ENSG[0-9]+).*", "\\\\1", df\$id.exposure)
        nm <- gtf_reduced\$name[match(ens, gtf_reduced\$id)]
        nm[is.na(nm)] <- df\$id.exposure[is.na(nm)]
        df\$geneName <- nm
        df
    }
    
    # Create summary statistics
    summary_text <- sprintf("
    ===== MR Analysis Summary for ${params.tissue} =====
    
    Total Number of Genes Tested: %d
    Number of Significant Associations: %d
    MR Methods Used: %s
    
    ===== Colocalization Summary =====
    
    Genes Tested: %d
    Genes with H4 > 0.8: %d
    ", 
    length(unique(mr_combined\$id.exposure)),
    nrow(mr_combined),
    paste(unique(mr_combined\$method), collapse = ", "),
    nrow(coloc),
    sum(coloc\$H4 > 0.8, na.rm = TRUE))
    writeLines(summary_text, "${params.tissue}_summary_stats.txt")
    
    # ---- Plot 1: MR effect sizes ----
    if (nrow(mr_combined) > 0) {
        p1 <- ggplot(mr_combined, aes(x = b, y = -log10(pval))) +
            geom_point(alpha = 0.6) +
            geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
            labs(title = sprintf("MR Effect Sizes - %s", "${params.tissue}"),
                 x = "Effect Size (Beta)", y = "-log10(P-value)") +
            theme_minimal() +
            theme(plot.title = element_text(hjust = 0.5, face = "bold"))
        ggsave("01_mr_manhattan_plot.pdf", p1, width = 10, height = 6)
        ggsave("01_mr_manhattan_plot.png", p1, width = 10, height = 6, dpi = 300)
    }
    
    # ---- Plot 2: coloc H4 distribution ----
    if (nrow(coloc) > 0) {
        p2 <- ggplot(coloc, aes(x = H4)) +
            geom_histogram(binwidth = 0.05, fill = "steelblue", alpha = 0.7) +
            geom_vline(xintercept = 0.8, linetype = "dashed", color = "red") +
            labs(title = sprintf("Colocalization H4 Distribution - %s", "${params.tissue}"),
                 x = "H4 Posterior Probability", y = "Count") +
            theme_minimal() +
            theme(plot.title = element_text(hjust = 0.5, face = "bold"))
        ggsave("02_coloc_h4_distribution.pdf", p2, width = 10, height = 6)
        ggsave("02_coloc_h4_distribution.png", p2, width = 10, height = 6, dpi = 300)
    }
    
    # ---- Plot 3: method comparison ----
    if (nrow(mr_combined) > 0 && length(unique(mr_combined\$method)) > 1) {
        p3 <- ggplot(mr_combined, aes(x = method, y = -log10(pval))) +
            geom_boxplot(fill = "lightblue", alpha = 0.7) +
            geom_jitter(width = 0.2, alpha = 0.4) +
            geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
            labs(title = sprintf("MR Method Comparison - %s", "${params.tissue}"),
                 x = "Method", y = "-log10(P-value)") +
            theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                  plot.title = element_text(hjust = 0.5, face = "bold"))
        ggsave("03_mr_method_comparison.pdf", p3, width = 10, height = 6)
        ggsave("03_mr_method_comparison.png", p3, width = 10, height = 6, dpi = 300)
    }
    
    # ---- Plot 4: coloc hypothesis summary ----
    if (nrow(coloc) > 0) {
        h_cols <- c("H0", "H1", "H2", "H3", "H4")
        h_data <- coloc %>%
            pivot_longer(all_of(h_cols), names_to = "Hypothesis", values_to = "Probability") %>%
            mutate(Hypothesis = factor(Hypothesis, levels = h_cols))
        p4 <- ggplot(h_data, aes(x = Hypothesis, y = Probability, fill = Hypothesis)) +
            geom_boxplot(alpha = 0.7) +
            scale_fill_brewer(palette = "Set2") +
            labs(title = sprintf("Colocalization Hypotheses - %s", "${params.tissue}"),
                 x = "Hypothesis", y = "Posterior Probability") +
            theme_minimal() +
            theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")
        ggsave("04_coloc_hypothesis_summary.pdf", p4, width = 10, height = 6)
        ggsave("04_coloc_hypothesis_summary.png", p4, width = 10, height = 6, dpi = 300)
    }
    
    # ---- Plot 5: FOREST PLOT of significant genes (with gene names) ----
    if (nrow(mr_combined) > 0) {
        # Genome-wide FDR to define significance for the forest plot
        mr_combined\$pval_fdr <- p.adjust(mr_combined\$pval, method = "BH")
        signif <- mr_combined %>% filter(pval_fdr < 0.05)
        
        cat(sprintf("Forest plot: %d genes at FDR < 0.05\\n", nrow(signif)))
        
        if (nrow(signif) > 0) {
            # OR and 95% CI (same formulas as your original script)
            signif\$OR <- exp(signif\$b)
            signif\$CI_lower <- signif\$b - signif\$se * qnorm(0.975)
            signif\$CI_upper <- signif\$b + signif\$se * qnorm(0.975)
            
            # Map to gene names
            signif <- add_gene_names(signif)
            
            # Order by effect size; make labels unique so duplicate symbols don't collapse
            signif <- signif[order(signif\$b), ]
            signif\$geneLabel <- make.unique(as.character(signif\$geneName))
            signif\$geneLabel <- factor(signif\$geneLabel, levels = signif\$geneLabel)
            
            p_forest <- ggplot(signif, aes(x = b, y = geneLabel, xmin = CI_lower, xmax = CI_upper)) +
                geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
                geom_errorbarh(height = 0.25, color = "steelblue") +
                geom_point(size = 2, color = "darkblue") +
                labs(title = sprintf("Significant MR Genes (FDR < 0.05) - %s", "${params.tissue}"),
                     x = "Effect (Beta)", y = "Gene") +
                theme_minimal() +
                theme(plot.title = element_text(hjust = 0.5, face = "bold"))
            
            h <- max(4, nrow(signif) * 0.35)
            ggsave("05_forest_plot_significant.pdf", p_forest, width = 10, height = h, limitsize = FALSE)
            ggsave("05_forest_plot_significant.png", p_forest, width = 10, height = h, dpi = 300, limitsize = FALSE)
            
            saveRDS(signif, "significant_results_with_names.rds")
            write.csv(signif, "significant_results_with_names.csv", row.names = FALSE)
        } else {
            cat("No FDR-significant genes for forest plot\\n")
        }
    }
    
    cat("Plotting complete\\n")
    """
}
