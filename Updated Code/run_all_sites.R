# Call the MCMC function
source("mcmc.R") 

# Define Core MCMC and Data Parameters
site_names <- unique(genotypedata_latefailures$Site)

# Loop per site analysis
for (site in site_names) {
  jobname <- site
  
  # Site-Specific Data Preparation
  genotypedata_RR <- genotypedata_latefailures[genotypedata_latefailures$Site == site, -c(2)]
  additional_neutral <- additional_genotypedata[additional_genotypedata$Site == site, -c(2)]
  maxMOI = max(as.numeric(sapply(1:length(colnames(genotypedata_RR)), function (x)
    strsplit(colnames(genotypedata_RR)[x],"_")[[1]][2])),na.rm=TRUE)
  ids = unique(unlist(strsplit(genotypedata_RR$Sample.ID[grepl("Day 0",genotypedata_RR$Sample.ID)]," Day 0")))
  locinames = unique(sapply(colnames(genotypedata_RR)[-1],function(x) strsplit(x,"_")[[1]][1]))
  nloci = length(locinames)
  nids = length(ids)
  maxalleles=30
  k = rep(maxalleles, nloci)
  alleles_definitions_RR  = define_alleles(rbind(genotypedata_RR,additional_neutral),locirepeats,k)
  
  # Automated MCMC Loop
  full_loglik_history <- list()
  full_chain_results <- list()  
  total_iterations <- 0
  converged <- FALSE
  
  while (!converged && total_iterations < max_iterations) {
    total_iterations <- total_iterations + chunk_size
    
    chunk_results <- future_lapply(1:n_chains, function(id) {
      run_one_chain(
        chain_id = id, nruns = chunk_size, burnin = 0, record_interval = 10,
        nids = nids, nloci = nloci, maxMOI = maxMOI, locinames = locinames,
        genotypedata_RR = genotypedata_RR, 
        additional_neutral = additional_neutral, 
        alleles_definitions_RR = alleles_definitions_RR
      )
    }, future.seed = TRUE)
    
    # Append results to the full history
    if (length(full_chain_results) == 0) {
      full_chain_results <- chunk_results
      full_loglik_history <- lapply(chunk_results, `[[`, "loglikelihood")
    } else {
      for (i in 1:n_chains) {
        full_chain_results[[i]]$parameters <- cbind(full_chain_results[[i]]$parameters, chunk_results[[i]]$parameters)
        full_chain_results[[i]]$classification <- cbind(full_chain_results[[i]]$classification, chunk_results[[i]]$classification)
        full_chain_results[[i]]$alleles0 <- abind(full_chain_results[[i]]$alleles0, chunk_results[[i]]$alleles0, along = 3)
        full_chain_results[[i]]$allelesf <- abind(full_chain_results[[i]]$allelesf, chunk_results[[i]]$allelesf, along = 3)
        full_loglik_history[[i]] <- c(full_loglik_history[[i]], chunk_results[[i]]$loglikelihood)
      }
    }
    
    # Robust Convergence Check
    if (length(full_loglik_history) == 0 || length(full_loglik_history[[1]]) == 0) {
      cat("Log-likelihood history is not populated yet. Running another chunk...\n")
      next
    }
    
    mcmc_list_loglik <- coda::mcmc.list(lapply(full_loglik_history, coda::mcmc))
    n_samples <- nrow(mcmc_list_loglik[[1]])
    burn_in_end <- floor(burn_in_frac * n_samples)
    
    if (is.null(n_samples) || (n_samples - burn_in_end < 50)) { next }
    
    post_burn_mcmc <- window(mcmc_list_loglik, start = burn_in_end + 1)
    r_hat <- try(coda::gelman.diag(post_burn_mcmc)$psrf[1, 1], silent = TRUE)
    ess <- try(coda::effectiveSize(post_burn_mcmc)[1], silent = TRUE)
    
    if (inherits(r_hat, "try-error") || inherits(ess, "try-error")) { next }
    r_hat_ok <- !is.na(r_hat) && r_hat < R_hat_threshold
    ess_ok <- !is.na(ess) && ess > ESS_threshold
    
    cat(sprintf("--- Convergence Check (Log-Likelihood Only) ---\n"))
    cat(sprintf("  R-hat: %.4f (Threshold: < %.2f) -> %s\n", r_hat, R_hat_threshold, ifelse(r_hat_ok, "OK", "FAIL")))
    cat(sprintf("  ESS:   %.1f (Threshold: > %d) -> %s\n", ess, ESS_threshold, ifelse(ess_ok, "OK", "FAIL")))
    
    if (r_hat_ok && ess_ok) {
      converged <- TRUE
      cat(sprintf("\nSUCCESS: Model for site '%s' CONVERGED after %d total iterations.\n", site, total_iterations))
    } else if (total_iterations >= max_iterations) {
      cat(sprintf("\nWARNING: Model for site '%s' FAILED TO CONVERGE after reaching max iterations.\n", site, max_iterations))
    }
  }
  
  # Post-processing and Saving Results for the CURRENT SITE
  num_total_samples_per_chain <- ifelse(length(full_chain_results) > 0, ncol(full_chain_results[[1]]$parameters), 0)
  if (num_total_samples_per_chain == 0) {
    next
  }
  
  burn_in_samples_per_chain <- floor(burn_in_frac * num_total_samples_per_chain)
  if (burn_in_samples_per_chain >= num_total_samples_per_chain) {
    cat(sprintf("WARNING: Burn-in resulted in 0 samples for site '%s'. Skipping analysis.\n", site))
    next
  }
  
  keep_indices <- (burn_in_samples_per_chain + 1):num_total_samples_per_chain
  
  final_classification <- do.call(cbind, lapply(full_chain_results, function(x) x$classification[, keep_indices, drop=FALSE]))
  final_parameters <- do.call(cbind, lapply(full_chain_results, function(x) x$parameters[, keep_indices, drop=FALSE]))
  final_loglikelihood <- lapply(full_loglik_history, function(x) x[keep_indices])
  final_alleles0 <- do.call(abind::abind, c(lapply(full_chain_results, function(x) x$alleles0[,, keep_indices, drop=FALSE]), along = 3))
  final_allelesf <- do.call(abind::abind, c(lapply(full_chain_results, function(x) x$allelesf[,, keep_indices, drop=FALSE]), along = 3))
  
  generate_likelihood_diagnostics(all_chains_loglikelihood = final_loglikelihood, site_name = site)
  
  modealleles <- matrix("", 2 * nids, maxMOI * nloci)
  for (i in 1:nids) {
    for (j in 1:nloci) {
      for (x in 1:maxMOI) {
        idx <- (j - 1) * maxMOI + x
        tbl0 <- table(final_alleles0[i, idx, ])
        if (length(tbl0) > 0) modealleles[2 * (i - 1) + 1, idx] <- names(tbl0)[which.max(tbl0)]
        tblf <- table(final_allelesf[i, idx, ])
        if (length(tblf) > 0) modealleles[2 * (i - 1) + 2, idx] <- names(tblf)[which.max(tblf)]
      }
    }
  }
  
  rowMeans2 <- function(x) { if (is.null(dim(x))) mean(x) else rowMeans(x) }
  
  prob_rec_combined <- rowMeans2(final_classification)
  temp_combined <- rep(prob_rec_combined, each = 2)
  
  outputmatrix <- cbind(temp_combined, modealleles)
  colnames(outputmatrix) <- c("Prob Rec", sapply(1:nloci, function(x) paste(locinames[x], "_", 1:maxMOI, sep = "")))
  write.csv(outputmatrix, paste(jobname, "_posterior", ".csv", sep = ""))
  write.csv(final_parameters, paste(jobname, "_state_parameters.csv", sep = ""))
  summary_statisticsmatrix <- cbind(
    format(rowMeans(final_parameters), digits = 2),
    apply(format(t(sapply(1:nrow(final_parameters), function(x) quantile(final_parameters[x, ], c(0.25, 0.75)))), digits = 2), 1, function(x) paste(x, collapse = "–"))
  )
  mean_diversity_samples <- colMeans(final_parameters[(3 + nloci):(3 + 2 * nloci - 1), ])
  summary_statisticsmatrix <- rbind(
    summary_statisticsmatrix,
    c(format(mean(mean_diversity_samples), digits = 2), paste(format(quantile(mean_diversity_samples, c(0.25, 0.75)), digits = 2), collapse = "–"))
  )
  summary_statisticsmatrix <- as.matrix(sapply(1:nrow(summary_statisticsmatrix), function(x) paste(summary_statisticsmatrix[x, 1], " (", summary_statisticsmatrix[x, 2], ")", sep = "")))
  rownames(summary_statisticsmatrix) <- c("q", "d", locinames, locinames, "Mean diversity")
  write.csv(summary_statisticsmatrix, paste(jobname, "_summarystatistics.csv", sep = ""))
  
  all_sites_classification[[site]] <- final_classification
  all_sites_parameters[[site]] <- final_parameters
  all_sites_alleles0[[site]] <- final_alleles0
  all_sites_allelesf[[site]] <- final_allelesf
  all_sites_ids[[site]] <- ids
  
} 
plan(sequential)