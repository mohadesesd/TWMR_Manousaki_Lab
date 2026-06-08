#!/usr/bin/env Rscript

# Test LDlink API connectivity and proxy batch function

library(LDlinkR)
library(data.table)

cat(strrep("=", 60), "\n")
cat("LDlink API Diagnostic Test\n")
cat(strrep("=", 60), "\n\n")

# Your token
token <- "23a23730fd0b"
cat("Token: ", token, "\n")

# Test 1: Simple API call with a single SNP
cat("TEST 1: Single SNP query (rs7412)\n")
cat(strrep("-", 40), "\n")

tryCatch({
  cat("Attempting to query rs7412 via LDproxy()...\n")
  result <- LDproxy(snp = "rs7412", pop = "ALL", token = token, r2d = "r2")
  cat("✓ Success! LDproxy works\n")
  print(head(result))
}, error = function(e) {
  cat("✗ Error:\n")
  cat("Message:", e$message, "\n")
})

cat("\n\n")

# Test 2: Batch API call
cat("TEST 2: Batch query (LDproxy_batch)\n")
cat(strrep("-", 40), "\n")

test_snps <- c("rs7412", "rs429358", "rs7259620")
cat("Testing with SNPs:", paste(test_snps, collapse=", "), "\n\n")

tryCatch({
  cat("Attempting LDproxy_batch()...\n")
  result <- LDproxy_batch(
    snps = test_snps,
    token = token,
    genome_build = "grch38",
    pop = "ALL",
    append = TRUE
  )
  cat("✓ Success! LDproxy_batch works\n")
  cat("Result type:", class(result), "\n")
  
  # Check what files were created
  cat("\nFiles created:\n")
  files <- list.files(".", pattern = "\\.txt$")
  if (length(files) > 0) {
    for (f in files) {
      size <- file.size(f)
      cat("  -", f, "(", size, "bytes)\n")
    }
  } else {
    cat("  (No .txt files created)\n")
  }
  
}, error = function(e) {
  cat("✗ Error in LDproxy_batch:\n")
  cat("Message:", e$message, "\n\n")
  
  cat("Troubleshooting:\n")
  cat("1. Check token is valid (32 chars):", nchar(token) == 32, "\n")
  cat("2. Check internet connection\n")
  cat("3. Check if LDlink server is up: https://ldlink.nci.nih.gov/\n")
  cat("4. Try with single SNP query above\n")
})

cat("\n")
cat(strrep("=", 60), "\n")
cat("Diagnostic complete.\n")
cat(strrep("=", 60), "\n")
