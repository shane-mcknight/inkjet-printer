################################################################################
# 2PL Gini Index Anchor Selection - Working Prototype
# Focus: Core implementation of Gini method extended to 2PL
################################################################################

# Load essential packages
library(mirt)
library(tidyverse)
library(ineq)

# ==============================================================================
# PART 1: DATA GENERATION
# ==============================================================================

#' Generate 2PL data with known DIF
#' @param n_items Number of items
#' @param n_ref Reference group sample size
#' @param n_foc Focal group sample size
#' @param prop_dif Proportion of DIF items
#' @param dif_size_b Magnitude of DIF in difficulty
#' @param dif_size_a Magnitude of DIF in discrimination
#' @param seed Random seed
generate_2pl_data <- function(n_items = 20, 
                              n_ref = 500, 
                              n_foc = 300,
                              prop_dif = 0.2,
                              dif_size_b = 0.5,
                              dif_size_a = 0.3,
                              seed = 123) {
  set.seed(seed)
  
  # Generate true parameters
  a_true <- runif(n_items, 0.5, 2.0)  # Discrimination parameters
  b_true <- rnorm(n_items, 0, 1)      # Difficulty parameters
  
  # Determine which items have DIF
  n_dif <- round(n_items * prop_dif)
  dif_items <- sample(1:n_items, n_dif)
  
  # Create focal group parameters with DIF
  a_foc <- a_true
  b_foc <- b_true
  a_foc[dif_items] <- a_true[dif_items] + dif_size_a
  b_foc[dif_items] <- b_true[dif_items] + dif_size_b
  
  # Generate ability parameters
  theta_ref <- rnorm(n_ref, 0, 1)
  theta_foc <- rnorm(n_foc, 0, 1)
  
  # Generate response data using 2PL model
  generate_responses <- function(theta, a, b) {
    n_person <- length(theta)
    n_item <- length(a)
    responses <- matrix(0, n_person, n_item)
    
    for (i in 1:n_person) {
      for (j in 1:n_item) {
        prob <- 1 / (1 + exp(-a[j] * (theta[i] - b[j])))
        responses[i, j] <- rbinom(1, 1, prob)
      }
    }
    return(responses)
  }
  
  # Generate data for both groups
  data_ref <- generate_responses(theta_ref, a_true, b_true)
  data_foc <- generate_responses(theta_foc, a_foc, b_foc)
  
  # Combine data
  data_combined <- rbind(data_ref, data_foc)
  group <- c(rep("reference", n_ref), rep("focal", n_foc))
  
  # Return everything
  return(list(
    data = data_combined,
    group = group,
    dif_items = dif_items,
    true_params = list(
      a_ref = a_true, b_ref = b_true,
      a_foc = a_foc, b_foc = b_foc
    ),
    n_ref = n_ref,
    n_foc = n_foc
  ))
}

# ==============================================================================
# PART 2: PARAMETER ESTIMATION
# ==============================================================================

#' Estimate 2PL parameters for both groups
#' @param data Response matrix
#' @param group Group membership vector
estimate_2pl_params <- function(data, group) {
  
  # Fit multigroup 2PL model
  model <- multipleGroup(data, 
                         model = 1,  # 2PL model
                         group = group,
                         SE = TRUE,
                         invariance = c('free_means', 'free_var'))
  
  # Extract parameters for reference group
  coef_ref <- coef(model, group = 1, simplify = TRUE)
  a_ref <- coef_ref$items[, "a1"]
  b_ref <- -coef_ref$items[, "d"] / coef_ref$items[, "a1"]  # Convert to difficulty
  
  # Extract parameters for focal group
  coef_foc <- coef(model, group = 2, simplify = TRUE)
  a_foc <- coef_foc$items[, "a1"]
  b_foc <- -coef_foc$items[, "d"] / coef_foc$items[, "a1"]
  
  # Get standard errors (using delta method for b)
  se_ref <- coef(model, group = 1, printSE = TRUE)
  se_foc <- coef(model, group = 2, printSE = TRUE)
  
  return(list(
    a_ref = a_ref,
    b_ref = b_ref,
    a_foc = a_foc,
    b_foc = b_foc,
    model = model,
    se_ref = se_ref,
    se_foc = se_foc
  ))
}

# ==============================================================================
# PART 3: GINI INDEX CALCULATION
# ==============================================================================

#' Calculate Gini index for 2PL anchor selection
#' @param params Estimated parameters
#' @param shift_a Shift for discrimination parameters
#' @param shift_b Shift for difficulty parameters
#' @param weight_a Weight for discrimination (0-1)
calculate_gini_2pl <- function(params, shift_a = 0, shift_b = 0, weight_a = 0.5) {
  
  # Apply shifts
  a_ref_shifted <- params$a_ref + shift_a
  b_ref_shifted <- params$b_ref + shift_b
  
  # Calculate differences
  diff_a <- params$a_foc - a_ref_shifted
  diff_b <- params$b_foc - b_ref_shifted
  
  # Standardize differences
  std_diff_a <- abs(diff_a) / sd(diff_a, na.rm = TRUE)
  std_diff_b <- abs(diff_b) / sd(diff_b, na.rm = TRUE)
  
  # Combined distance metric
  distances <- weight_a * std_diff_a + (1 - weight_a) * std_diff_b
  
  # Calculate Gini index
  gini_value <- ineq::Gini(distances)
  
  return(list(
    gini = gini_value,
    distances = distances
  ))
}

#' Optimize shifts to maximize Gini index
#' @param params Estimated parameters
#' @param weight_a Weight for discrimination
#' @param n_starts Number of random starting values
optimize_gini_2pl <- function(params, weight_a = 0.5, n_starts = 10) {
  
  # Objective function (we want to maximize, so negate)
  objective <- function(shifts) {
    -calculate_gini_2pl(params, shifts[1], shifts[2], weight_a)$gini
  }
  
  # Try multiple starting values
  best_result <- NULL
  best_value <- Inf
  
  for (i in 1:n_starts) {
    # Random starting values
    start_vals <- rnorm(2, 0, 0.1)
    
    # Optimize
    result <- tryCatch({
      optim(start_vals, objective, method = "L-BFGS-B",
            lower = c(-2, -2), upper = c(2, 2))
    }, error = function(e) NULL)
    
    if (!is.null(result) && result$value < best_value) {
      best_value <- result$value
      best_result <- result
    }
  }
  
  # Get final Gini and distances with optimal shifts
  final_gini <- calculate_gini_2pl(params, 
                                   best_result$par[1], 
                                   best_result$par[2], 
                                   weight_a)
  
  # Identify anchor items (lowest 25% of distances)
  n_items <- length(params$a_ref)
  n_anchor <- ceiling(n_items * 0.25)
  anchor_items <- order(final_gini$distances)[1:n_anchor]
  
  return(list(
    optimal_shifts = best_result$par,
    gini_value = -best_result$value,
    distances = final_gini$distances,
    anchor_items = anchor_items,
    convergence = best_result$convergence
  ))
}

# ==============================================================================
# PART 4: DIF TESTING
# ==============================================================================

#' Perform Wald test for DIF
#' @param model mirt model object
#' @param anchor_items Vector of anchor item indices
wald_dif_test <- function(model, anchor_items) {
  n_items <- extract.mirt(model, "nitems")
  test_items <- setdiff(1:n_items, anchor_items)
  
  # Build constraint matrix for testing DIF
  results <- data.frame(
    item = integer(),
    chi_sq = numeric(),
    df = integer(),
    p_value = numeric()
  )
  
  for (item in test_items) {
    # Test equality of both a and d parameters
    constrain <- paste0("a1_g1 - a1_g2 = 0\n",
                        "d_g1 - d_g2 = 0")
    
    test <- tryCatch({
      DIF(model, which.par = c(paste0("a1_", item), paste0("d_", item)),
          groups = c("reference", "focal"))
    }, error = function(e) NULL)
    
    if (!is.null(test)) {
      results <- rbind(results, data.frame(
        item = item,
        chi_sq = test$X2,
        df = test$df,
        p_value = test$p
      ))
    }
  }
  
  # Apply multiple testing correction
  results$p_adjusted <- p.adjust(results$p_value, method = "BH")
  results$dif_detected <- results$p_adjusted < 0.05
  
  return(results)
}

# ==============================================================================
# PART 5: COMPARISON METHODS
# ==============================================================================

#' All items as anchors (traditional approach)
all_items_anchor <- function(params) {
  n_items <- length(params$a_ref)
  return(1:n_items)
}

#' Random anchor selection
random_anchor <- function(params, prop_anchor = 0.25, seed = 123) {
  set.seed(seed)
  n_items <- length(params$a_ref)
  n_anchor <- ceiling(n_items * prop_anchor)
  return(sample(1:n_items, n_anchor))
}

#' Component Loss Function (CLF) method
clf_anchor <- function(params, n_anchor = NULL) {
  n_items <- length(params$a_ref)
  if (is.null(n_anchor)) n_anchor <- ceiling(n_items * 0.25)
  
  # Calculate CLF criterion for each item
  diff_a <- params$a_foc - params$a_ref
  diff_b <- params$b_foc - params$b_ref
  
  # CLF uses squared differences
  clf_values <- diff_a^2 + diff_b^2
  
  # Select items with smallest CLF values
  anchor_items <- order(clf_values)[1:n_anchor]
  return(anchor_items)
}

# ==============================================================================
# PART 6: PERFORMANCE EVALUATION
# ==============================================================================

#' Evaluate anchor selection performance
#' @param selected_anchors Vector of selected anchor items
#' @param true_dif_items Vector of true DIF items
#' @param n_items Total number of items
evaluate_performance <- function(selected_anchors, true_dif_items, n_items) {
  
  # True anchor items (non-DIF items)
  true_anchors <- setdiff(1:n_items, true_dif_items)
  
  # Calculate performance metrics
  tp <- length(intersect(selected_anchors, true_anchors))  # True positives
  fp <- length(intersect(selected_anchors, true_dif_items))  # False positives
  tn <- length(intersect(setdiff(1:n_items, selected_anchors), true_dif_items))  # True negatives
  fn <- length(intersect(setdiff(1:n_items, selected_anchors), true_anchors))  # False negatives
  
  # Calculate metrics
  sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA
  specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA
  precision <- if ((tp + fp) > 0) tp / (tp + fp) else NA
  f1 <- if (!is.na(precision) && !is.na(sensitivity) && (precision + sensitivity) > 0) {
    2 * (precision * sensitivity) / (precision + sensitivity)
  } else NA
  
  return(list(
    tp = tp, fp = fp, tn = tn, fn = fn,
    sensitivity = sensitivity,
    specificity = specificity,
    precision = precision,
    f1_score = f1
  ))
}

# ==============================================================================
# PART 7: MAIN SIMULATION FUNCTION
# ==============================================================================

#' Run single simulation iteration
run_single_simulation <- function(n_items = 20,
                                  n_ref = 500,
                                  n_foc = 300,
                                  prop_dif = 0.2,
                                  dif_size_b = 0.5,
                                  dif_size_a = 0.3,
                                  weight_a = 0.5,
                                  seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Generate data
  cat("Generating data...\n")
  sim_data <- generate_2pl_data(n_items, n_ref, n_foc, prop_dif, 
                                dif_size_b, dif_size_a, seed)
  
  # Estimate parameters
  cat("Estimating 2PL parameters...\n")
  params <- estimate_2pl_params(sim_data$data, sim_data$group)
  
  # Apply different anchor selection methods
  cat("Applying anchor selection methods...\n")
  
  # 1. Gini method
  gini_result <- optimize_gini_2pl(params, weight_a)
  
  # 2. All items
  all_anchors <- all_items_anchor(params)
  
  # 3. Random selection
  random_anchors <- random_anchor(params)
  
  # 4. CLF method
  clf_anchors <- clf_anchor(params)
  
  # Evaluate performance
  cat("Evaluating performance...\n")
  perf_gini <- evaluate_performance(gini_result$anchor_items, 
                                    sim_data$dif_items, n_items)
  perf_all <- evaluate_performance(all_anchors, 
                                   sim_data$dif_items, n_items)
  perf_random <- evaluate_performance(random_anchors, 
                                      sim_data$dif_items, n_items)
  perf_clf <- evaluate_performance(clf_anchors, 
                                   sim_data$dif_items, n_items)
  
  # Compile results
  results <- list(
    conditions = list(
      n_items = n_items,
      n_ref = n_ref,
      n_foc = n_foc,
      prop_dif = prop_dif,
      dif_size_b = dif_size_b,
      dif_size_a = dif_size_a,
      weight_a = weight_a
    ),
    true_dif = sim_data$dif_items,
    anchors = list(
      gini = gini_result$anchor_items,
      all = all_anchors,
      random = random_anchors,
      clf = clf_anchors
    ),
    performance = list(
      gini = perf_gini,
      all = perf_all,
      random = perf_random,
      clf = perf_clf
    ),
    gini_details = list(
      optimal_shifts = gini_result$optimal_shifts,
      gini_value = gini_result$gini_value,
      distances = gini_result$distances
    )
  )
  
  return(results)
}

# ==============================================================================
# PART 8: VISUALIZATION FUNCTIONS
# ==============================================================================

#' Plot parameter space with anchor visualization
plot_parameter_space <- function(results, params) {
  
  # Create data frame for plotting
  plot_data <- data.frame(
    item = 1:length(params$a_ref),
    a_ref = params$a_ref,
    b_ref = params$b_ref,
    a_foc = params$a_foc,
    b_foc = params$b_foc,
    is_dif = 1:length(params$a_ref) %in% results$true_dif,
    is_anchor_gini = 1:length(params$a_ref) %in% results$anchors$gini
  )
  
  # Create plot
  p <- ggplot(plot_data) +
    # Reference group
    geom_point(aes(x = b_ref, y = a_ref, color = "Reference"), 
               size = 3, alpha = 0.7) +
    # Focal group
    geom_point(aes(x = b_foc, y = a_foc, color = "Focal"), 
               size = 3, alpha = 0.7) +
    # Connect same items
    geom_segment(aes(x = b_ref, y = a_ref, xend = b_foc, yend = a_foc),
                 alpha = 0.3) +
    # Highlight DIF items
    geom_point(data = filter(plot_data, is_dif),
               aes(x = b_ref, y = a_ref), 
               shape = 21, size = 5, stroke = 2, color = "red") +
    # Highlight selected anchors
    geom_point(data = filter(plot_data, is_anchor_gini),
               aes(x = b_ref, y = a_ref), 
               shape = 24, size = 4, fill = "green", color = "darkgreen") +
    scale_color_manual(values = c("Reference" = "blue", "Focal" = "orange")) +
    labs(title = "2PL Parameter Space with Anchor Selection",
         x = "Difficulty (b)",
         y = "Discrimination (a)",
         color = "Group") +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  return(p)
}

#' Plot Lorenz curve for Gini index
plot_lorenz_curve <- function(distances) {
  
  # Sort distances
  sorted_dist <- sort(distances)
  n <- length(sorted_dist)
  
  # Calculate cumulative proportions
  cum_prop <- cumsum(sorted_dist) / sum(sorted_dist)
  item_prop <- (1:n) / n
  
  # Create data frame
  lorenz_data <- data.frame(
    item_prop = c(0, item_prop),
    cum_prop = c(0, cum_prop)
  )
  
  # Calculate Gini
  gini_val <- ineq::Gini(distances)
  
  # Create plot
  p <- ggplot(lorenz_data, aes(x = item_prop, y = cum_prop)) +
    geom_line(size = 1.5, color = "blue") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    geom_ribbon(aes(ymin = item_prop, ymax = cum_prop), 
                fill = "lightblue", alpha = 0.3) +
    annotate("text", x = 0.7, y = 0.3, 
             label = paste("Gini =", round(gini_val, 3)),
             size = 5) +
    labs(title = "Lorenz Curve for Item Distances",
         x = "Cumulative Proportion of Items",
         y = "Cumulative Proportion of Distance") +
    theme_minimal()
  
  return(p)
}

# ==============================================================================
# PART 9: EXAMPLE USAGE
# ==============================================================================

# Run a single simulation
cat("Running example simulation...\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

example_result <- run_single_simulation(
  n_items = 20,
  n_ref = 500,
  n_foc = 300,
  prop_dif = 0.3,
  dif_size_b = 0.5,
  dif_size_a = 0.3,
  weight_a = 0.5,
  seed = 42
)

# Print results
cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("RESULTS SUMMARY\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

cat("\nTrue DIF items:", example_result$true_dif, "\n")

cat("\nSelected anchors:\n")
cat("  Gini method:", example_result$anchors$gini, "\n")
cat("  CLF method:", example_result$anchors$clf, "\n")
cat("  Random:", example_result$anchors$random, "\n")

cat("\nPerformance metrics:\n")
methods <- c("Gini", "CLF", "Random", "All")
for (method in tolower(methods)) {
  perf <- example_result$performance[[method]]
  cat(sprintf("\n%s method:\n", method))
  cat(sprintf("  Sensitivity: %.3f\n", perf$sensitivity))
  cat(sprintf("  Specificity: %.3f\n", perf$specificity))
  cat(sprintf("  Precision: %.3f\n", perf$precision))
  cat(sprintf("  F1 Score: %.3f\n", perf$f1_score))
}

cat("\nGini optimization details:\n")
cat("  Optimal shifts (a, b):", 
    round(example_result$gini_details$optimal_shifts, 4), "\n")
cat("  Gini value:", round(example_result$gini_details$gini_value, 4), "\n")

# Create visualizations
cat("\nCreating visualizations...\n")

# Regenerate data and params for plotting
sim_data <- generate_2pl_data(seed = 42)
params <- estimate_2pl_params(sim_data$data, sim_data$group)

# Parameter space plot
p1 <- plot_parameter_space(example_result, params)

# Lorenz curve
p2 <- plot_lorenz_curve(example_result$gini_details$distances)

# Display plots (if in interactive session)
if (interactive()) {
  print(p1)
  print(p2)
}

cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("Prototype simulation complete!\n")
cat(paste(rep("=", 50), collapse = ""), "\n")