library(matrixStats)
clonal_synergy_core <- function(R1, R2, RC, type = "bliss", sizes = NULL, t = 1, scaling = TRUE) {
  # Ensure sizes
  if (is.null(sizes) || length(sizes) == 0) sizes <- rep(1, length(R1))
  
  # Optional scaling (no centering)
  if (scaling) {
    R1 <- scale(R1, center = FALSE)
    R2 <- scale(R2, center = FALSE)
    RC <- scale(RC, center = FALSE)
  }
  
  # Clone-level expected combo (per clone)
  if (type == "bliss") {
    R_comb <- R1 + R2 - RC
  } else if (type == "HSA") {
    R_comb <- pmin(R1, R2)
  } else if (type == "LSA") {
    R_comb <- pmax(R1, R2)
  } else if (type == "x_only") {
    R_comb <- R1
  } else if (type == "y_only") {
    R_comb <- R2
  } else {
    stop("Unknown type")
  }
  
  # Population totals
  cellsC  <- sum(sizes * exp(RC * t))
  cellsR1 <- sum(sizes * exp(R1 * t))
  cellsR2 <- sum(sizes * exp(R2 * t))
  
  comb_total <- log(sum(sizes * exp(R_comb * t)) / sum(sizes)) / t
  RC_total   <- log(cellsC / sum(sizes)) / t
  R1_total   <- log(cellsR1 / sum(sizes)) / t
  R2_total   <- log(cellsR2 / sum(sizes)) / t
  
  # Population-level expected combo
  if (type == "bliss") {
    R_pop_combo <- R1_total + R2_total - RC_total
  } else if (type == "HSA") {
    R_pop_combo <- pmin(R1_total, R2_total)
  } else if (type == "LSA") {
    R_pop_combo <- pmax(R1_total, R2_total)
  } else if (type == "x_only") {
    R_pop_combo <- R1_total
  } else if (type == "y_only") {
    R_pop_combo <- R2_total
  } else {
    stop("Unknown type")
  }
  
  # Clonal-synergy quantities
  clonal_synergy        <- (R_pop_combo - comb_total)
  clonal_synergy_ratio  <- clonal_synergy / abs(RC_total)
  
  list(
    expected_r_pop     = R_pop_combo,
    expected_r_clonal  = comb_total,
    RC_total           = RC_total,
    R1_total           = R1_total,
    R2_total           = R2_total,
    clonal_synergy     = clonal_synergy,
    clonal_synergy_ratio = clonal_synergy_ratio
  )
}

clonal_synergy <- function(
    R1, R2, RC,
    type = "bliss", sizes = NULL, t = 1, scaling = FALSE,
    calculate_permutation_results = T, n_perm = 1000,
    synergy_type = "clonal_synergy_ratio", eps = 1e-12,
    calculate_relative_to_Baseline = "relative" # or "relative"
) {
  results <- clonal_synergy_core(R1, R2, RC, type = type, sizes = sizes, t = t, scaling = scaling)
  
  if ( type!="bliss" & type!="HSA" ) {
    # set calculate_relative_to_Baseline to FALSE
    calculate_relative_to_Baseline = FALSE
    cat("Note: calculate_relative_to_Baseline is only supported for 'bliss' and 'HSA' types. Setting to FALSE.\n")
  }
  if (calculate_permutation_results) {
    perm_results <- replicate(n_perm, {
      # Permute clone-wise differential effects relative to control and re-add control
      R1_yhat <- sample(R1 - RC)
      R2_yhat <- sample(R2 - RC)
      R1perm  <- RC + R1_yhat
      R2perm  <- RC + R2_yhat
      clonal_synergy_core(R1perm, R2perm, RC, type = type, sizes = sizes, t = t, scaling = scaling)[[synergy_type]]
    })
    
    if (calculate_relative_to_Baseline == "relative") {
      m <- mean(perm_results)
      perm_results <- perm_results - m
      results[[synergy_type]] <- results[[synergy_type]] - m
    }
    
    sd_val <- stats::sd(perm_results)
    m      <- mean(perm_results)
    med    <- stats::median(perm_results)
    madv   <- stats::mad(perm_results, constant = 1.4826)
    
    z_sd   <- (results[[synergy_type]] - m)   / max(sd_val, eps)
    z_mad  <- (results[[synergy_type]] - med) / max(madv,   eps)
    
    p_upper <- (sum(perm_results >= results[[synergy_type]]) + 1) / (n_perm + 1)
    p_lower <- (sum(perm_results <= results[[synergy_type]]) + 1) / (n_perm + 1)
    
    o <- list(
      expected_r_pop     = results$expected_r_pop,
      expected_r_clonal  = results$expected_r_clonal,
      RC_total           = results$RC_total,
      R1_total           = results$R1_total,
      R2_total           = results$R2_total,
      clonal_synergy     = results[[synergy_type]],
      p_synergy          = p_upper,
      p_antagonism       = p_lower,
      z                  = z_sd,
     # z_mad              = z_mad,
      perm_results       = perm_results,
      sd                 = sd_val,
      synergy_type    = synergy_type
    )
   # names(o)[which(names(o) == "clonal_synergy")] <- permutation_var
    o
  } else {
    o <- list(
      expected_r_pop     = results$expected_r_pop,
      expected_r_clonal  = results$expected_r_clonal,
      RC_total           = results$RC_total,
      R1_total           = results$R1_total,
      R2_total           = results$R2_total,
      clonal_synergy     = results[[synergy_type]],
      p_synergy          = NA,
      p_antagonism       = NA,
      z                  = NA,
      #z_mad              = NA,
      perm_results       = NA,
      sd                 = NA,
      synergy_type    = synergy_type
    )
    #names(o)[which(names(o) == "clonal_synergy")] <- permutation_var
    o
  }
}
.stouffer  <- function(z) {sum(z) / sqrt(length(z) ) }

clonal_interaction_z=function(
    R1, R2, RC, sizes, t,synergy_type="clonal_synergy_ratio",calculate_relative_to_Baseline="relative",scaling= FALSE,
    n_perm = 1000, eps = 1e-12, dev.output=F
) {

  bliss=clonal_synergy(
    R1, R2, RC,
    type = "bliss", sizes = sizes, t = t, scaling = scaling,
    calculate_permutation_results = T, n_perm = n_perm,
    synergy_type = synergy_type, eps = eps,
    calculate_relative_to_Baseline = calculate_relative_to_Baseline # or "relative"
  )
  
  HSA=clonal_synergy(
    R1, R2, RC,
    type = "HSA", sizes = sizes, t = t, scaling = scaling,
    calculate_permutation_results = T, n_perm = n_perm,
    synergy_type = synergy_type, eps = eps,
    calculate_relative_to_Baseline = calculate_relative_to_Baseline # or "relative"
  )
  if (dev.output){
    return(c(bliss=bliss, HSA=HSA))
  } else {
    return(c(bliss_z=bliss$z,HSA_z=HSA$z,zS=.stouffer(c(bliss$z,HSA$z)), HSA =HSA$clonal_synergy, bliss=bliss$clonal_synergy,cS=(bliss$clonal_synergy+HSA$clonal_synergy)/2  )  )
  }

}


clonal_interaction_z=function(
    R1, R2, RC, sizes, t,synergy_type="clonal_synergy_ratio",calculate_relative_to_Baseline="relative",scaling= FALSE,
    n_perm = 100, eps = 1e-12, dev.output=F
) {
  
  bliss=clonal_synergy(
    R1, R2, RC,
    type = "bliss", sizes = sizes, t = t, scaling = scaling,
    calculate_permutation_results = T, n_perm = n_perm,
    synergy_type = synergy_type, eps = eps,
    calculate_relative_to_Baseline = calculate_relative_to_Baseline # or "relative"
  )
  
  HSA=clonal_synergy(
    R1, R2, RC,
    type = "HSA", sizes = sizes, t = t, scaling = scaling,
    calculate_permutation_results = T, n_perm = n_perm,
    synergy_type = synergy_type, eps = eps,
    calculate_relative_to_Baseline = calculate_relative_to_Baseline # or "relative"
  )
  if (dev.output){
    return(c(bliss=bliss, HSA=HSA))
  } else {
    return(c(bliss_z=bliss$z,HSA_z=HSA$z,zS=.stouffer(c(bliss$z,HSA$z)), HSA =HSA$clonal_synergy, bliss=bliss$clonal_synergy,cS=(bliss$clonal_synergy+HSA$clonal_synergy)/2  )  )
  }
  
}

clonal_z_bliss=function(
    R1, R2, RC, sizes, t , n_perm = 100, eps = 1e-12, dev.output=F
) {
  
  bliss=clonal_synergy(
    R1, R2, RC,
    type = "bliss", sizes = sizes, t = t, scaling = FALSE,
    calculate_permutation_results = T, n_perm = n_perm,
    synergy_type = "clonal_synergy_ratio", eps = eps,
    calculate_relative_to_Baseline = "absolute" # or "relative"
  )
  
  
  if (dev.output){
    return(c(bliss=bliss))
  } else {
    return(c(bliss_z=bliss$z, bliss=bliss$clonal_synergy )  )
  }
  
}
