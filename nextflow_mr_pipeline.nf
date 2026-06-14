#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.tissue = "Whole_Blood"
params.exposure_file = null
params.outcome_file = null
params.outdir = "results"
params.sample_size = 670
params.fstat_threshold = 10
params.p_threshold = 0.05
params.lddlink_token = "23a23730fd0b"
params.ld_threshold = 0.8
params.genome_build = "grch38"   // LDlink genome build: grch37 | grch38 | grch38_high_coverage
params.population = "ALL"        // LDlink reference population (e.g. EUR, ALL)
params.coloc_h4_threshold = 0.8
params.window_size = 500000
params.chromosomes = 1..22
params.gene_conversion_table = "../GeneCode_ConversionTable.rds"

include { EXPOSURE_PREPROCESSING } from './modules/exposure_preprocessing.nf'
include { FIND_MISSING_SNPS; PROXY_BATCH_LDLINK; PROCESS_PROXY_RESULTS; MERGE_PROXIED_DATA } from './modules/proxy_batch.nf'
include { OUTCOME_PREPROCESSING } from './modules/outcome_preprocessing.nf'
include { MR_ANALYSIS } from './modules/mr_analysis.nf'
include { COLOCALIZATION } from './modules/colocalization.nf'
include { RESULTS_PLOTTING } from './modules/results_plotting.nf'

workflow {
    main:
        if (!params.exposure_file || !params.outcome_file) {
            error "ERROR: --exposure_file and --outcome_file are required"
        }
        
        log.info "Starting MR Pipeline for ${params.tissue}..."
        
        chromosomes_ch = Channel.from(params.chromosomes)
        
        exposure_preprocessed = EXPOSURE_PREPROCESSING(
            params.exposure_file, params.tissue, params.sample_size,
            params.fstat_threshold, params.outdir
        )
        
        missing_snps = FIND_MISSING_SNPS(
            exposure_preprocessed.exposure_by_chr.collect(),
            params.outcome_file, params.tissue
        )
        
        proxy_results = PROXY_BATCH_LDLINK(
            missing_snps.missing_snps, params.lddlink_token, params.genome_build, params.population
        )
        
        proxies_processed = PROCESS_PROXY_RESULTS(
            proxy_results.proxy_results, params.outcome_file, params.ld_threshold
        )
        
        proxied_outcome_data = MERGE_PROXIED_DATA(
            proxies_processed.proxies_processed,
            missing_snps.missing_snps,
            params.outcome_file,
            params.ld_threshold
        )
        
        outcome_preprocessed = OUTCOME_PREPROCESSING(
            params.outcome_file, proxied_outcome_data.outcome_with_proxies,
            params.tissue, chromosomes_ch, params.outdir
        )
        
        log.info "Running MR Analysis..."
        exposure_files = exposure_preprocessed.exposure_by_chr
            .flatten()
            .map { f -> def m = (f.getName() =~ /exposure_by_chr_(\d+)\.rds/); tuple(m[0][1] as Integer, f) }
        outcome_files = outcome_preprocessed.outcome_by_chr
            .flatten()
            .map { f -> def m = (f.getName() =~ /outcome_chr_(\d+)\.rds/); tuple(m[0][1] as Integer, f) }
        mr_input = exposure_files.join(outcome_files)
        
        mr_results = MR_ANALYSIS(
            mr_input, params.tissue, params.p_threshold, params.outdir
        )
        
        log.info "Running colocalization..."
        coloc_results = COLOCALIZATION(
            mr_results.harmonized_data.collect(),
            mr_results.filtered_results.collect(),
            exposure_preprocessed.exposure_raw,
            outcome_preprocessed.outcome_raw,
            params.tissue,
            params.window_size,
            params.coloc_h4_threshold,
            params.outdir
        )
        
        log.info "Generating plots..."
        final_plots = RESULTS_PLOTTING(
            mr_results.all_results.collect(),
            coloc_results.coloc_table,
            file(params.gene_conversion_table),
            params.tissue,
            params.outdir
        )
        
        log.info "Pipeline completed!"
}
