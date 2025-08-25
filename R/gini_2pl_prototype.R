# Test script for 2PL Gini Index Prototype
# Save this as: test_prototype.R

# Clear environment
rm(list = ls())

# Source the prototype
cat("Loading prototype functions...\n")
source("R/01_gini_prototype.R")  # Adjust path as needed

# Test 1: Simple simulation
cat("\n", strrep("=", 50), "\n")
cat("TEST 1: Basic functionality\n")
cat(strrep("=", 50), "\n")

result1 <- tryCatch({
  run_single_simulation(
    n_items = 10,      # Fewer items for quick test
    n_ref = 200,       # Smaller sample
    n_foc = 200,       # Balanced groups
    prop_dif = 0.2,    # 20% DIF
    dif_size_b = 0.4,
    dif_size_a = 0.2,
    weight_a = 0.5,
    seed = 123
  )
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  NULL
})

if (!is.null(result1)) {
  cat("✓ Basic simulation successful\n")
  cat("  True DIF items:", result1$true_dif, "\n")
  cat("  Gini anchors:", result1$anchors$gini, "\n")
  cat("  Gini value:", round(result1$gini_details$gini_value, 3), "\n")
} else {
  cat("✗ Basic simulation failed\n")
}

# Test 2: Different sample sizes
cat("\n", strrep("=", 50), "\n")
cat("TEST 2: Unbalanced groups\n")
cat(strrep("=", 50), "\n")

result2 <- tryCatch({
  run_single_simulation(
    n_items = 15,
    n_ref = 500,
    n_foc = 150,  # Unbalanced
    prop_dif = 0.3,
    seed = 456
  )
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  NULL
})

if (!is.null(result2)) {
  cat("✓ Unbalanced groups successful\n")
  
  # Compare methods
  cat("\nMethod comparison (F1 scores):\n")
  cat("  Gini:", round(result2$performance$gini$f1_score, 3), "\n")
  cat("  CLF:", round(result2$performance$clf$f1_score, 3), "\n")
  cat("  Random:", round(result2$performance$random$f1_score, 3), "\n")
} else {
  cat("✗ Unbalanced groups failed\n")
}

# Test 3: Extreme conditions
cat("\n", strrep("=", 50), "\n")
cat("TEST 3: Extreme DIF\n")
cat(strrep("=", 50), "\n")

result3 <- tryCatch({
  run_single_simulation(
    n_items = 20,
    n_ref = 300,
    n_foc = 300,
    prop_dif = 0.5,  # 50% DIF items
    dif_size_b = 0.8,  # Large DIF
    dif_size_a = 0.5,
    seed = 789
  )
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  NULL
})

if (!is.null(result3)) {
  cat("✓ Extreme DIF successful\n")
  cat("  Proportion correct anchors (Gini):", 
      round(result3$performance$gini$precision, 3), "\n")
} else {
  cat("✗ Extreme DIF failed\n")
}

# Summary
cat("\n", strrep("=", 50), "\n")
cat("TEST SUMMARY\n")
cat(strrep("=", 50), "\n")

tests_passed <- sum(!is.null(result1), !is.null(result2), !is.null(result3))
cat("Tests passed:", tests_passed, "/ 3\n")

if (tests_passed == 3) {
  cat("\n✓ All tests passed! Prototype is working correctly.\n")
  
  # Save a test result for inspection
  saveRDS(result1, "test_result.rds")
  cat("Test result saved to 'test_result.rds'\n")
  
  # Create simple comparison table
  cat("\n", strrep("=", 50), "\n")
  cat("PERFORMANCE COMPARISON\n")
  cat(strrep("=", 50), "\n")
  
  methods <- c("gini", "clf", "random")
  metrics <- data.frame(
    Method = c("Gini", "CLF", "Random"),
    Sensitivity = sapply(methods, function(m) 
      round(result1$performance[[m]]$sensitivity, 3)),
    Specificity = sapply(methods, function(m) 
      round(result1$performance[[m]]$specificity, 3)),
    F1_Score = sapply(methods, function(m) 
      round(result1$performance[[m]]$f1_score, 3))
  )
  print(metrics, row.names = FALSE)
  
} else {
  cat("\n✗ Some tests failed. Check error messages above.\n")
}

cat("\n", strrep("=", 50), "\n")
cat("Testing complete!\n")
