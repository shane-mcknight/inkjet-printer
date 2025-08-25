################################################################################
# REFINED 2PL GINI METHOD WITH SINGLE ANCHOR SELECTION
# Literature-based comparison methods
# 
# CLF implementation based on:
# Asparouhov, T., & Muthén, B. (2014). Multiple-group factor analysis alignment.
# Structural Equation Modeling: A Multidisciplinary Journal, 21(4), 495-508.
################################################################################

library(mirt)
library(ineq)

# ==============================================================================
# STEP 1: Generate Simple 2PL Data with DIF
# ==============================================================================
generate_data <- function(n_ref = 1500, n_foc = 1000, n_items = 20, prop_dif = 0.3, seed = 123) {
  set.seed(seed)
  
  # True parameters
  a_true <- runif(n_items, 0.8, 1.5)  # Discrimination
  b_true <- rnorm(n_items, 0, 1)      # Difficulty
  
  # Which items have DIF?
  n_dif_items <- round(n_items * prop_dif)
  dif_items <- sample(1:n_items, n_dif_items)
  
  # Create DIF (moderate magnitudes for realistic scenario)
  a_foc <- a_true
  b_foc <- b_true
  a_foc[dif_items] <- a_true[dif_items] + rnorm(n_dif_items, 0.4, 0.1)  # DIF in discrimination
  b_foc[dif_items] <- b_true[dif_items] + rnorm(n_dif_items, 0.6, 0.15) # DIF in difficulty
  
  # Generate abilities
  theta_ref <- rnorm(n_ref, 0, 1)
  theta_foc <- rnorm(n_foc, 0, 1)  # Could add mean difference here
  
  # Generate responses (2PL model)
  # P(X=1) = 1 / (1 + exp(-a*(theta - b)))
  
  # Reference group
  data_ref <- matrix(0, n_ref, n_items)
  colnames(data_ref) <- paste0("Item", 1:n_items)
  for (i in 1:n_ref) {
    for (j in 1:n_items) {
      prob <- 1 / (1 + exp(-a_true[j] * (theta_ref[i] - b_true[j])))
      data_ref[i, j] <- rbinom(1, 1, prob)
    }
  }
  
  # Focal group  
  data_foc <- matrix(0, n_foc, n_items)
  colnames(data_foc) <- paste0("Item", 1:n_items)
  for (i in 1:n_foc) {
    for (j in 1:n_items) {
      prob <- 1 / (1 + exp(-a_foc[j] * (theta_foc[i] - b_foc[j])))
      data_foc[i, j] <- rbinom(1, 1, prob)
    }
  }
  
  # Combine
  data_all <- rbind(data_ref, data_foc)
  colnames(data_all) <- paste0("Item", 1:n_items)
  group <- c(rep("ref", n_ref), rep("foc", n_foc))
  
  cat("Data generated:\n")
  cat("  Reference group:", n_ref, "subjects\n")
  cat("  Focal group:", n_foc, "subjects\n")
  cat("  Items:", n_items, "\n")
  cat("  True DIF items:", dif_items, "\n\n")
  
  return(list(
    data = data_all,
    group = group,
    dif_items = dif_items,
    n_items = n_items,
    true_params = list(
      a_ref = a_true,
      b_ref = b_true,
      a_foc = a_foc,
      b_foc = b_foc
    )
  ))
}

# ==============================================================================
# STEP 2: Estimate Parameters with Standard Errors
# ==============================================================================
estimate_params <- function(data, group) {
  
  cat("Estimating 2PL parameters...\n")
  
  # Separate estimation for each group
  data_ref <- data[group == "ref", ]
  data_foc <- data[group == "foc", ]
  
  # Ensure column names are present
  if (is.null(colnames(data_ref))) {
    colnames(data_ref) <- paste0("Item", 1:ncol(data_ref))
  }
  if (is.null(colnames(data_foc))) {
    colnames(data_foc) <- paste0("Item", 1:ncol(data_foc))
  }
  
  # Fit 2PL for reference group
  mod_ref <- mirt(data_ref, model = 1, itemtype = "2PL", verbose = FALSE, 
                  technical = list(NCYCLES = 1000))
  params_ref <- coef(mod_ref, simplify = TRUE)$items
  a_ref <- params_ref[, "a1"]
  b_ref <- -params_ref[, "d"] / params_ref[, "a1"]
  
  # Get standard errors for reference group
  se_ref <- coef(mod_ref, printSE = TRUE)
  se_a_ref <- numeric(length(a_ref))
  se_b_ref <- numeric(length(b_ref))
  for(i in 1:length(a_ref)) {
    item_se <- se_ref[[i]]
    if(!is.null(item_se) && nrow(item_se) >= 2) {
      se_a_ref[i] <- item_se[1, 2]  # SE for a parameter
      se_d <- item_se[2, 2]  # SE for d parameter
      # Delta method for SE of b = -d/a
      se_b_ref[i] <- abs(se_d / a_ref[i])  # Simplified delta method
    }
  }
  
  # Fit 2PL for focal group
  mod_foc <- mirt(data_foc, model = 1, itemtype = "2PL", verbose = FALSE,
                  technical = list(NCYCLES = 1000))
  params_foc <- coef(mod_foc, simplify = TRUE)$items
  a_foc <- params_foc[, "a1"]
  b_foc <- -params_foc[, "d"] / params_foc[, "a1"]
  
  # Get standard errors for focal group
  se_foc <- coef(mod_foc, printSE = TRUE)
  se_a_foc <- numeric(length(a_foc))
  se_b_foc <- numeric(length(b_foc))
  for(i in 1:length(a_foc)) {
    item_se <- se_foc[[i]]
    if(!is.null(item_se) && nrow(item_se) >= 2) {
      se_a_foc[i] <- item_se[1, 2]
      se_d <- item_se[2, 2]
      se_b_foc[i] <- abs(se_d / a_foc[i])
    }
  }
  
  cat("  Parameters estimated successfully\n")
  cat("  Reference convergence:", mod_ref@OptimInfo$converged, "\n")
  cat("  Focal convergence:", mod_foc@OptimInfo$converged, "\n\n")
  
  return(list(
    a_ref = a_ref,
    b_ref = b_ref,
    a_foc = a_foc,
    b_foc = b_foc,
    se_a_ref = se_a_ref,
    se_b_ref = se_b_ref,
    se_a_foc = se_a_foc,
    se_b_foc = se_b_foc,
    n_items = ncol(data),
    mod_ref = mod_ref,
    mod_foc = mod_foc
  ))
}

# ==============================================================================
# STEP 3: GINI METHOD - Calculate Gini Index
# ==============================================================================
calculate_gini <- function(params, shift_a = 0, shift_b = 0, weight_a = 0.5) {
  
  # Apply shifts to reference group
  a_ref_shifted <- params$a_ref + shift_a
  b_ref_shifted <- params$b_ref + shift_b
  
  # Calculate differences
  diff_a <- abs(params$a_foc - a_ref_shifted)
  diff_b <- abs(params$b_foc - b_ref_shifted)
  
  # Standardize (avoid division by zero)
  sd_a <- sd(diff_a)
  sd_b <- sd(diff_b)
  if(sd_a < 1e-10) sd_a <- 1
  if(sd_b < 1e-10) sd_b <- 1
  
  diff_a_std <- diff_a / sd_a
  diff_b_std <- diff_b / sd_b
  
  # Combined distance (weighted)
  distances <- weight_a * diff_a_std + (1 - weight_a) * diff_b_std
  
  # Calculate Gini coefficient
  gini <- ineq::Gini(distances)
  
  return(list(gini = gini, distances = distances))
}

# Find Optimal Shifts for Gini
find_optimal_shifts_gini <- function(params, weight_a = 0.5) {
  
  cat("GINI METHOD: Finding optimal shifts...\n")
  
  # Objective function (negative because optim minimizes)
  objective <- function(shifts) {
    -calculate_gini(params, shifts[1], shifts[2], weight_a)$gini
  }
  
  # Try multiple starting values
  best_result <- NULL
  best_value <- Inf
  
  start_values <- list(
    c(0, 0),
    c(0.1, 0.1),
    c(-0.1, -0.1),
    c(0.1, -0.1),
    c(-0.1, 0.1)
  )
  
  for(start in start_values) {
    result <- optim(
      par = start,
      fn = objective,
      method = "Nelder-Mead",
      control = list(maxit = 1000)
    )
    
    if(result$value < best_value) {
      best_value <- result$value
      best_result <- result
    }
  }
  
  # Get final Gini and distances
  final <- calculate_gini(params, best_result$par[1], best_result$par[2], weight_a)
  
  cat("  Optimal shift_a:", round(best_result$par[1], 3), "\n")
  cat("  Optimal shift_b:", round(best_result$par[2], 3), "\n")
  cat("  Gini value:", round(final$gini, 3), "\n")
  
  return(list(
    shifts = best_result$par,
    gini = final$gini,
    distances = final$distances
  ))
}

# ==============================================================================
# STEP 4: COMPONENT LOSS FUNCTION (CLF) METHOD - Asparouhov & Muthén (2014)
# ==============================================================================
clf_method <- function(params, weight_a = 0.5, epsilon = 0.001) {
  
  cat("CLF METHOD: Calculating component loss function (Asparouhov & Muthén, 2014)...\n")
  
  # This implements the Component Loss Function from:
  # Asparouhov, T., & Muthén, B. (2014). Multiple-group factor analysis alignment.
  # Structural Equation Modeling, 21(4), 495-508.
  # 
  # The CLF is: f(x) = sqrt(x^2 + epsilon) - sqrt(epsilon)
  # This function is more robust to outliers than squared loss
  # 
  # Adapted here for anchor selection in 2PL models:
  # - Each item is tested as a potential anchor
  # - Scales are aligned based on that anchor
  # - Total loss is calculated for all other items
  # - Item with minimum loss is selected
  
  n_items <- params$n_items
  clf_values <- numeric(n_items)
  
  # Define the component loss function f(x) = sqrt(x^2 + epsilon) - sqrt(epsilon)
  # This is the specific CLF from Asparouhov & Muthén (2014)
  clf_function <- function(x, eps = epsilon) {
    sqrt(x^2 + eps) - sqrt(eps)
  }
  
  # For each item as potential anchor
  for(i in 1:n_items) {
    # Use item i as the anchor for alignment
    # The shift is based on the difference for this anchor item
    shift_a <- params$a_foc[i] - params$a_ref[i]
    shift_b <- params$b_foc[i] - params$b_ref[i]
    
    # Apply these shifts to align all reference parameters
    a_ref_aligned <- params$a_ref + shift_a
    b_ref_aligned <- params$b_ref + shift_b
    
    # Calculate differences for all OTHER items (exclude anchor)
    other_items <- setdiff(1:n_items, i)
    
    diff_a <- params$a_foc[other_items] - a_ref_aligned[other_items]
    diff_b <- params$b_foc[other_items] - b_ref_aligned[other_items]
    
    # Apply the component loss function to each difference
    # Note: Asparouhov & Muthén apply CLF to raw differences, not standardized
    loss_a <- sapply(diff_a, clf_function)
    loss_b <- sapply(diff_b, clf_function)
    
    # Total loss is the weighted sum
    # Weight by parameter importance (discrimination vs difficulty)
    clf_values[i] <- sum(weight_a * loss_a) + sum((1 - weight_a) * loss_b)
  }
  
  # Item with minimum CLF value is the best anchor
  # (minimizes the loss when used for alignment)
  anchor_item <- which.min(clf_values)
  
  cat("  CLF values range:", round(min(clf_values), 3), "-", round(max(clf_values), 3), "\n")
  cat("  Selected anchor:", anchor_item, "\n")
  cat("  Epsilon used:", epsilon, "\n")
  
  return(list(
    clf_values = clf_values,
    anchor = anchor_item,
    epsilon = epsilon
  ))
}

# ==============================================================================
# STEP 5: ITERATIVE FORWARD SELECTION METHOD
# ==============================================================================
iterative_forward_method <- function(params, weight_a = 0.5) {
  
  cat("ITERATIVE FORWARD METHOD: Selecting anchor...\n")
  
  n_items <- params$n_items
  criterion_values <- numeric(n_items)
  
  for(item in 1:n_items) {
    # Use this item as anchor - align scales based on it
    shift_a <- params$a_foc[item] - params$a_ref[item]
    shift_b <- params$b_foc[item] - params$b_ref[item]
    
    # Apply alignment to all items
    a_ref_aligned <- params$a_ref + shift_a
    b_ref_aligned <- params$b_ref + shift_b
    
    # Calculate total squared differences for all items after alignment
    all_diff_a <- params$a_foc - a_ref_aligned
    all_diff_b <- params$b_foc - b_ref_aligned
    
    # Weight by standard errors if available
    if(!is.null(params$se_a_ref) && !is.null(params$se_a_foc) && 
       all(params$se_a_ref > 0) && all(params$se_a_foc > 0)) {
      # Use pooled SEs as weights (items with larger SEs get less weight)
      se_pooled_a <- sqrt(params$se_a_ref^2 + params$se_a_foc^2)
      se_pooled_b <- sqrt(params$se_b_ref^2 + params$se_b_foc^2)
      
      # Avoid division by zero
      se_pooled_a[se_pooled_a < 0.01] <- 0.01
      se_pooled_b[se_pooled_b < 0.01] <- 0.01
      
      # Weight differences by inverse of SEs
      weighted_diff_a <- all_diff_a / se_pooled_a
      weighted_diff_b <- all_diff_b / se_pooled_b
    } else {
      # No weighting if SEs not available
      weighted_diff_a <- all_diff_a
      weighted_diff_b <- all_diff_b
    }
    
    # Standardize
    sd_a <- sd(weighted_diff_a)
    sd_b <- sd(weighted_diff_b)
    if(sd_a > 0) weighted_diff_a <- weighted_diff_a / sd_a
    if(sd_b > 0) weighted_diff_b <- weighted_diff_b / sd_b
    
    # Criterion is weighted sum of squared standardized differences
    criterion_values[item] <- sum(weight_a * weighted_diff_a^2 + 
                                    (1 - weight_a) * weighted_diff_b^2)
  }
  
  # Select item with minimum criterion
  best_item <- which.min(criterion_values)
  
  cat("  Criterion range:", round(min(criterion_values), 3), "-", 
      round(max(criterion_values), 3), "\n")
  cat("  Selected anchor:", best_item, "\n")
  
  return(list(
    anchor = best_item,
    criterion_values = criterion_values
  ))
}

# ==============================================================================
# STEP 6: RANDOM SELECTION METHOD (Baseline)
# ==============================================================================
random_method <- function(n_items, seed = NULL) {
  
  cat("RANDOM METHOD: Selecting anchor...\n")
  
  if(!is.null(seed)) set.seed(seed)
  anchor <- sample(1:n_items, 1)
  
  cat("  Selected anchor:", anchor, "\n")
  
  return(list(anchor = anchor))
}

# ==============================================================================
# STEP 7: ALL ITEMS METHOD (Reference - uses all items)
# ==============================================================================
all_items_method <- function(n_items) {
  
  cat("ALL ITEMS METHOD: Using all items as anchors...\n")
  cat("  Note: This method doesn't select a single anchor\n")
  
  return(list(anchors = 1:n_items))
}

# ==============================================================================
# STEP 8: EVALUATE SINGLE ANCHOR PERFORMANCE
# ==============================================================================
evaluate_single_anchor <- function(selected_anchor, true_dif_items, n_items) {
  
  # Is the selected anchor truly DIF-free?
  is_correct <- !(selected_anchor %in% true_dif_items)
  
  # True anchors are all non-DIF items
  true_anchors <- setdiff(1:n_items, true_dif_items)
  
  # Probability of correct selection by chance
  chance_prob <- length(true_anchors) / n_items
  
  cat("  Evaluation:\n")
  cat("    Anchor is DIF-free:", is_correct, "\n")
  cat("    True DIF-free items:", length(true_anchors), "out of", n_items, "\n")
  cat("    Chance probability:", round(chance_prob, 3), "\n")
  
  return(list(
    is_correct = is_correct,
    selected = selected_anchor,
    true_anchors = true_anchors,
    chance_prob = chance_prob
  ))
}

# ==============================================================================
# STEP 9: RUN COMPARISON ACROSS METHODS
# ==============================================================================
compare_methods <- function(sim_data, params, weight_a = 0.5, seed = 123) {
  
  cat("\n", strrep("=", 60), "\n")
  cat("COMPARING ANCHOR SELECTION METHODS\n")
  cat(strrep("=", 60), "\n\n")
  
  results <- list()
  
  # Method 1: Gini
  gini_result <- find_optimal_shifts_gini(params, weight_a)
  gini_anchor <- which.min(gini_result$distances)  # Single best anchor
  results$gini <- evaluate_single_anchor(gini_anchor, sim_data$dif_items, sim_data$n_items)
  results$gini$method <- "Gini"
  results$gini$distances <- gini_result$distances
  cat("\n")
  
  # Method 2: CLF (Asparouhov & Muthén, 2014)
  clf_result <- clf_method(params, weight_a, epsilon = 0.001)
  results$clf <- evaluate_single_anchor(clf_result$anchor, sim_data$dif_items, sim_data$n_items)
  results$clf$method <- "CLF"
  results$clf$clf_values <- clf_result$clf_values
  cat("\n")
  
  # Method 3: Iterative Forward
  iter_result <- iterative_forward_method(params, weight_a)
  results$iterative <- evaluate_single_anchor(iter_result$anchor, sim_data$dif_items, sim_data$n_items)
  results$iterative$method <- "Iterative"
  results$iterative$criterion_values <- iter_result$criterion_values
  cat("\n")
  
  # Method 4: Random
  random_result <- random_method(sim_data$n_items, seed)
  results$random <- evaluate_single_anchor(random_result$anchor, sim_data$dif_items, sim_data$n_items)
  results$random$method <- "Random"
  cat("\n")
  
  # Method 5: All items (for reference)
  all_result <- all_items_method(sim_data$n_items)
  results$all <- list(
    method = "All Items",
    anchors = all_result$anchors,
    is_correct = NA
  )
  
  return(results)
}

# ==============================================================================
# STEP 10: VISUALIZE RESULTS
# ==============================================================================
visualize_results <- function(results, true_dif_items, n_items) {
  
  cat("\n", strrep("=", 60), "\n")
  cat("RESULTS VISUALIZATION\n")
  cat(strrep("=", 60), "\n\n")
  
  # Create item status table
  item_table <- data.frame(
    Item = 1:n_items,
    True_Status = ifelse(1:n_items %in% true_dif_items, "DIF", "Anchor"),
    Gini_Distance = round(results$gini$distances, 3),
    CLF_Value = round(results$clf$clf_values, 3),
    Iter_Criterion = round(results$iterative$criterion_values, 3),
    Selected_Gini = ifelse(1:n_items == results$gini$selected, "<<<", ""),
    Selected_CLF = ifelse(1:n_items == results$clf$selected, "<<<", ""),
    Selected_Iter = ifelse(1:n_items == results$iterative$selected, "<<<", "")
  )
  
  # Sort by Gini distance
  item_table <- item_table[order(item_table$Gini_Distance), ]
  
  cat("Item Analysis Table (sorted by Gini distance):\n")
  cat("Top 10 items shown:\n\n")
  print(head(item_table, 10), row.names = FALSE)
  
  # Summary comparison
  cat("\n", strrep("-", 60), "\n")
  cat("SUMMARY COMPARISON:\n")
  cat(strrep("-", 60), "\n")
  
  summary_table <- data.frame(
    Method = c("Gini", "CLF", "Iterative", "Random"),
    Selected_Item = c(results$gini$selected, 
                      results$clf$selected,
                      results$iterative$selected,
                      results$random$selected),
    Is_Correct = c(results$gini$is_correct,
                   results$clf$is_correct,
                   results$iterative$is_correct,
                   results$random$is_correct)
  )
  
  print(summary_table, row.names = FALSE)
  
  # Agreement between methods
  cat("\n", strrep("-", 60), "\n")
  cat("METHOD AGREEMENT:\n")
  
  selected_items <- c(results$gini$selected, 
                      results$clf$selected,
                      results$iterative$selected)
  
  # Count how many methods agree
  agreement_table <- table(selected_items)
  
  if(max(agreement_table) == 3) {
    cat("All deterministic methods agree on item", names(agreement_table)[which.max(agreement_table)], "\n")
  } else {
    cat("Methods disagree:\n")
    for(method in c("Gini", "CLF", "Iterative")) {
      cat(sprintf("  %-10s selected: %2d\n", method, results[[tolower(method)]]$selected))
    }
    cat("\nConsensus items (selected by multiple methods):\n")
    consensus <- agreement_table[agreement_table > 1]
    if(length(consensus) > 0) {
      for(i in 1:length(consensus)) {
        cat("  Item", names(consensus)[i], "selected by", consensus[i], "methods\n")
      }
    } else {
      cat("  No consensus among methods\n")
    }
  }
}

# ==============================================================================
# STEP 11: MULTIPLE REPLICATIONS FOR SIMULATION STUDY
# ==============================================================================
run_simulation_study <- function(n_reps = 100, n_ref = 1500, n_foc = 1000, 
                                 n_items = 20, prop_dif = 0.3, weight_a = 0.5) {
  
  cat("\n", strrep("=", 60), "\n")
  cat("SIMULATION STUDY:", n_reps, "replications\n")
  cat(strrep("=", 60), "\n")
  cat("Settings: n_ref =", n_ref, ", n_foc =", n_foc, 
      ", n_items =", n_items, "\n")
  cat("         prop_dif =", prop_dif, ", weight_a =", weight_a, "\n\n")
  
  # Storage for results
  accuracy <- data.frame(
    gini = numeric(n_reps),
    clf = numeric(n_reps),
    iterative = numeric(n_reps),
    random = numeric(n_reps)
  )
  
  # Track convergence
  convergence_count <- 0
  
  pb <- txtProgressBar(min = 0, max = n_reps, style = 3)
  
  for(rep in 1:n_reps) {
    setTxtProgressBar(pb, rep)
    
    # Generate data
    sim_data <- generate_data(n_ref, n_foc, n_items, prop_dif, seed = rep * 100)
    
    # Estimate parameters
    params <- estimate_params(sim_data$data, sim_data$group)
    
    # Check convergence
    if(!is.null(params$mod_ref) && !is.null(params$mod_foc)) {
      if(params$mod_ref@OptimInfo$converged && params$mod_foc@OptimInfo$converged) {
        convergence_count <- convergence_count + 1
      }
    }
    
    # Run methods (suppress output)
    invisible(capture.output({
      results <- compare_methods(sim_data, params, weight_a, seed = rep * 100)
    }))
    
    # Store accuracy
    accuracy$gini[rep] <- results$gini$is_correct
    accuracy$clf[rep] <- results$clf$is_correct
    accuracy$iterative[rep] <- results$iterative$is_correct
    accuracy$random[rep] <- results$random$is_correct
  }
  
  close(pb)
  
  # Calculate accuracy rates
  cat("\n\n", strrep("-", 60), "\n")
  cat("SIMULATION RESULTS (", n_reps, "replications):\n")
  cat(strrep("-", 60), "\n")
  
  accuracy_rates <- colMeans(accuracy)
  
  cat("Convergence rate:", round(convergence_count/n_reps, 3), "\n\n")
  
  cat("Accuracy Rates (probability of selecting a true anchor):\n")
  cat("  Gini Method:      ", sprintf("%.3f", accuracy_rates["gini"]), "\n")
  cat("  CLF Method:       ", sprintf("%.3f", accuracy_rates["clf"]), "\n")
  cat("  Iterative Method: ", sprintf("%.3f", accuracy_rates["iterative"]), "\n")
  cat("  Random Method:    ", sprintf("%.3f", accuracy_rates["random"]), "\n")
  cat("  Chance Level:     ", sprintf("%.3f", 1 - prop_dif), "\n")
  
  # Statistical tests vs random
  cat("\n", strrep("-", 60), "\n")
  cat("Statistical Tests (vs Random):\n")
  
  for(method in c("gini", "clf", "iterative")) {
    test <- prop.test(c(sum(accuracy[[method]]), sum(accuracy$random)), 
                      c(n_reps, n_reps))
    sig <- ifelse(test$p.value < 0.001, "***",
                  ifelse(test$p.value < 0.01, "**",
                         ifelse(test$p.value < 0.05, "*", "")))
    cat(sprintf("  %-10s vs random: p = %.4f %s\n", method, test$p.value, sig))
  }
  
  # Pairwise comparisons between methods
  cat("\n", strrep("-", 60), "\n")
  cat("Pairwise Comparisons (McNemar tests):\n")
  
  methods <- c("gini", "clf", "iterative")
  for(i in 1:(length(methods)-1)) {
    for(j in (i+1):length(methods)) {
      # McNemar test for paired binary data
      contingency <- table(accuracy[[methods[i]]], accuracy[[methods[j]]])
      if(nrow(contingency) == 2 && ncol(contingency) == 2) {
        test <- mcnemar.test(contingency)
        sig <- ifelse(test$p.value < 0.05, "*", "")
        cat(sprintf("  %-10s vs %-10s: p = %.4f %s\n", 
                    methods[i], methods[j], test$p.value, sig))
      }
    }
  }
  
  return(list(
    accuracy = accuracy,
    accuracy_rates = accuracy_rates,
    convergence_rate = convergence_count/n_reps,
    settings = list(n_reps = n_reps, n_ref = n_ref, n_foc = n_foc,
                    n_items = n_items, prop_dif = prop_dif, weight_a = weight_a)
  ))
}

# ==============================================================================
# MAIN EXECUTION FUNCTIONS
# ==============================================================================

# Single demonstration
run_single_demo <- function() {
  cat("=" , strrep("=", 60), "\n")
  cat("SINGLE ANCHOR SELECTION - METHOD COMPARISON\n")
  cat(strrep("=", 61), "\n\n")
  
  # Step 1: Generate data
  sim_data <- generate_data(n_ref = 1500, n_foc = 1000, n_items = 20, prop_dif = 0.3)
  
  # Step 2: Estimate parameters
  params <- estimate_params(sim_data$data, sim_data$group)
  
  # Step 3: Compare methods
  results <- compare_methods(sim_data, params, weight_a = 0.5)
  
  # Step 4: Visualize results
  visualize_results(results, sim_data$dif_items, sim_data$n_items)
  
  return(results)
}

# Run simulation study with different conditions
run_condition_comparison <- function() {
  cat("\n", strrep("=", 70), "\n")
  cat("COMPARING DIFFERENT CONDITIONS\n")
  cat(strrep("=", 70), "\n\n")
  
  # Test different conditions
  conditions <- expand.grid(
    prop_dif = c(0.2, 0.3, 0.4),
    weight_a = c(0.3, 0.5, 0.7)
  )
  
  all_results <- list()
  
  for(i in 1:nrow(conditions)) {
    cat("\n", strrep("=", 70), "\n")
    cat("CONDITION", i, ": prop_dif =", conditions$prop_dif[i], 
        ", weight_a =", conditions$weight_a[i], "\n")
    cat(strrep("=", 70), "\n")
    
    result <- run_simulation_study(
      n_reps = 100,
      prop_dif = conditions$prop_dif[i],
      weight_a = conditions$weight_a[i]
    )
    
    all_results[[i]] <- result
  }
  
  # Summary table
  cat("\n\n", strrep("=", 70), "\n")
  cat("SUMMARY ACROSS ALL CONDITIONS\n")
  cat(strrep("=", 70), "\n\n")
  
  summary_df <- data.frame()
  for(i in 1:length(all_results)) {
    row <- data.frame(
      prop_dif = conditions$prop_dif[i],
      weight_a = conditions$weight_a[i],
      gini = all_results[[i]]$accuracy_rates["gini"],
      clf = all_results[[i]]$accuracy_rates["clf"],
      iterative = all_results[[i]]$accuracy_rates["iterative"],
      random = all_results[[i]]$accuracy_rates["random"]
    )
    summary_df <- rbind(summary_df, row)
  }
  
  print(summary_df, row.names = FALSE)
  
  return(all_results)
}

# ==============================================================================
# EXECUTE DEMONSTRATION
# ==============================================================================

# Run single demonstration
cat("\n", strrep("=", 70), "\n")
cat("REFINED 2PL GINI METHOD - SINGLE ANCHOR SELECTION\n")
cat("CLF based on Asparouhov & Muthén (2014)\n")
cat(strrep("=", 70), "\n\n")

demo_results <- run_single_demo()

# Run simulation study
cat("\n\n")
simulation_results <- run_simulation_study(n_reps = 100)

# Optional: Run full condition comparison
condition_results <- run_condition_comparison()

cat("\n", strrep("=", 70), "\n")
cat("ANALYSIS COMPLETE\n")
cat(strrep("=", 70), "\n")