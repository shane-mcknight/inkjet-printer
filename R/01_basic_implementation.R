################################################################################
# BASIC 2PL GINI METHOD - Minimal Working Implementation
# Large samples: 1500 reference, 1000 focal
################################################################################

library(mirt)
library(ineq)

# ==============================================================================
# STEP 1: Generate Simple 2PL Data with DIF
# ==============================================================================

generate_data <- function(n_ref = 1500, n_foc = 1000, n_items = 20, prop_dif = 0.3) {
  
  set.seed(123)
  
  # True parameters
  a_true <- runif(n_items, 0.8, 1.5)  # Discrimination
  b_true <- rnorm(n_items, 0, 1)      # Difficulty
  
  # Which items have DIF?
  n_dif_items <- round(n_items * prop_dif)
  dif_items <- sample(1:n_items, n_dif_items)
  
  # Create DIF (bigger magnitudes)
  a_foc <- a_true
  b_foc <- b_true
  a_foc[dif_items] <- a_true[dif_items] + 0.5  # Large DIF in discrimination
  b_foc[dif_items] <- b_true[dif_items] + 0.8  # Large DIF in difficulty
  
  # Generate abilities
  theta_ref <- rnorm(n_ref, 0, 1)
  theta_foc <- rnorm(n_foc, 0, 1)
  
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
  
  # Add column names (required by mirt)
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
    n_items = n_items
  ))
}

# ==============================================================================
# STEP 2: Estimate Parameters (Simple Approach)
# ==============================================================================

estimate_params <- function(data, group) {
  
  cat("Estimating 2PL parameters...\n")
  
  # Separate estimation for each group (simplest approach)
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
  mod_ref <- mirt(data_ref, model = 1, itemtype = "2PL", verbose = FALSE)
  params_ref <- coef(mod_ref, simplify = TRUE)$items
  a_ref <- params_ref[, "a1"]
  b_ref <- -params_ref[, "d"] / params_ref[, "a1"]
  
  # Fit 2PL for focal group
  mod_foc <- mirt(data_foc, model = 1, itemtype = "2PL", verbose = FALSE)
  params_foc <- coef(mod_foc, simplify = TRUE)$items
  a_foc <- params_foc[, "a1"]
  b_foc <- -params_foc[, "d"] / params_foc[, "a1"]
  
  cat("  Parameters estimated successfully\n\n")
  
  return(list(
    a_ref = a_ref,
    b_ref = b_ref,
    a_foc = a_foc,
    b_foc = b_foc
  ))
}

# ==============================================================================
# STEP 3: Calculate Gini Index
# ==============================================================================

calculate_gini <- function(params, shift_a = 0, shift_b = 0) {
  
  # Apply shifts to reference group
  a_ref_shifted <- params$a_ref + shift_a
  b_ref_shifted <- params$b_ref + shift_b
  
  # Calculate differences
  diff_a <- abs(params$a_foc - a_ref_shifted)
  diff_b <- abs(params$b_foc - b_ref_shifted)
  
  # Standardize
  diff_a_std <- diff_a / sd(diff_a)
  diff_b_std <- diff_b / sd(diff_b)
  
  # Combined distance (equal weight)
  distances <- (diff_a_std + diff_b_std) / 2
  
  # Calculate Gini
  gini <- ineq::Gini(distances)
  
  return(list(gini = gini, distances = distances))
}

# ==============================================================================
# STEP 4: Find Optimal Shifts
# ==============================================================================

find_optimal_shifts <- function(params) {
  
  cat("Finding optimal shifts to maximize Gini...\n")
  
  # Objective function (negative because optim minimizes)
  objective <- function(shifts) {
    -calculate_gini(params, shifts[1], shifts[2])$gini
  }
  
  # Optimize
  result <- optim(
    par = c(0, 0),  # Starting values
    fn = objective,
    method = "Nelder-Mead"
  )
  
  # Get final Gini and distances
  final <- calculate_gini(params, result$par[1], result$par[2])
  
  cat("  Optimal shifts found\n")
  cat("  Shift_a:", round(result$par[1], 3), "\n")
  cat("  Shift_b:", round(result$par[2], 3), "\n")
  cat("  Gini value:", round(final$gini, 3), "\n\n")
  
  return(list(
    shifts = result$par,
    gini = final$gini,
    distances = final$distances
  ))
}

# ==============================================================================
# STEP 5: Select Anchor Items
# ==============================================================================

select_anchors <- function(distances, prop_anchor = 0.25) {
  
  n_items <- length(distances)
  n_anchor <- ceiling(n_items * prop_anchor)
  
  # Items with smallest distances are anchors
  anchor_items <- order(distances)[1:n_anchor]
  
  cat("Anchor selection:\n")
  cat("  Number of anchors:", n_anchor, "\n")
  cat("  Selected items:", anchor_items, "\n\n")
  
  return(anchor_items)
}

# ==============================================================================
# STEP 6: Evaluate Performance
# ==============================================================================

evaluate <- function(selected_anchors, true_dif_items, n_items) {
  
  true_anchors <- setdiff(1:n_items, true_dif_items)
  
  # How many selected anchors are truly DIF-free?
  correct <- length(intersect(selected_anchors, true_anchors))
  
  # Metrics
  precision <- correct / length(selected_anchors)
  recall <- correct / length(true_anchors)
  
  cat("Performance:\n")
  cat("  Correct anchors:", correct, "out of", length(selected_anchors), "\n")
  cat("  Precision:", round(precision, 3), "\n")
  cat("  Recall:", round(recall, 3), "\n\n")
  
  return(list(precision = precision, recall = recall))
}

# ==============================================================================
# MAIN ANALYSIS
# ==============================================================================

cat("=" , strrep("=", 60), "\n")
cat("BASIC 2PL GINI METHOD - DEMONSTRATION\n")
cat(strrep("=", 61), "\n\n")

# Step 1: Generate data
sim_data <- generate_data(n_ref = 1500, n_foc = 1000, n_items = 20, prop_dif = 0.3)

# Step 2: Estimate parameters
params <- estimate_params(sim_data$data, sim_data$group)

# Step 3-4: Find optimal shifts and calculate Gini
gini_result <- find_optimal_shifts(params)

# Step 5: Select anchors
anchors <- select_anchors(gini_result$distances)

# Step 6: Evaluate
performance <- evaluate(anchors, sim_data$dif_items, sim_data$n_items)

# Summary
cat(strrep("=", 61), "\n")
cat("SUMMARY\n")
cat(strrep("=", 61), "\n")
cat("True DIF items:", sim_data$dif_items, "\n")
cat("Selected anchors:", anchors, "\n")
cat("Overlap with true anchors:", 
    length(intersect(anchors, setdiff(1:sim_data$n_items, sim_data$dif_items))), "\n")

# Simple visualization of distances
cat("\nItem distances (sorted):\n")
dist_df <- data.frame(
  Item = 1:length(gini_result$distances),
  Distance = round(gini_result$distances, 3),
  Is_DIF = 1:length(gini_result$distances) %in% sim_data$dif_items,
  Is_Anchor = 1:length(gini_result$distances) %in% anchors
)
dist_df <- dist_df[order(dist_df$Distance), ]
print(head(dist_df, 10))

cat("\n", strrep("=", 61), "\n")
cat("Basic implementation complete!\n")