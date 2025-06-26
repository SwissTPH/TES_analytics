# This file defines a function to perform a traditional match-counting analysis.

perform_match_counting <- function(genotypedata_latefailures, bin_sizes) {
  ids <- unique(unlist(strsplit(genotypedata_latefailures$Sample.ID[grepl("Day 0", genotypedata_latefailures$Sample.ID)], " Day 0")))
  locinames <- unique(sapply(colnames(genotypedata_latefailures)[-c(1, 2)], function(x) strsplit(x, "_")[[1]][1]))
  nloci <- length(locinames)
  nids <- length(ids)
  MOI0 <- rep(0, nids)
  MOIf <- rep(0, nids)
  for (i in 1:nids) {
    for (j in 1:nloci) {
      locicolumns <- grepl(paste(locinames[j], "_", sep = ""), colnames(genotypedata_latefailures))
      nalleles0 <- sum(!is.na(genotypedata_latefailures[grepl(paste(ids[i], "Day 0"), genotypedata_latefailures$Sample.ID), locicolumns]))
      nallelesf <- sum(!is.na(genotypedata_latefailures[grepl(paste(ids[i], "Day Failure"), genotypedata_latefailures$Sample.ID), locicolumns]))
      MOI0[i] <- max(MOI0[i], nalleles0)
      MOIf[i] <- max(MOIf[i], nallelesf)
    }
  }
  maxMOI <- max(c(MOI0, MOIf), na.rm = TRUE)
  alleles0 <- matrix(NA, nids, maxMOI * nloci)
  allelesf <- matrix(NA, nids, maxMOI * nloci)
  
  for (j in 1:nloci) {
    locus <- locinames[j]
    locicolumns <- grepl(paste0(locus, "_"), colnames(genotypedata_latefailures))
    oldalleles <- as.matrix(genotypedata_latefailures[, locicolumns])
    oldalleles <- apply(oldalleles, 2, as.numeric)
    day0_rows <- grepl("Day 0", genotypedata_latefailures$Sample.ID)
    dayf_rows <- grepl("Day Failure", genotypedata_latefailures$Sample.ID)
    alleles0[, (maxMOI * (j - 1) + 1):(maxMOI * (j - 1) + ncol(oldalleles))] <- oldalleles[day0_rows, , drop = FALSE]
    allelesf[, (maxMOI * (j - 1) + 1):(maxMOI * (j - 1) + ncol(oldalleles))] <- oldalleles[dayf_rows, , drop = FALSE]
  }

  number_matches <- rep(NA, nids)
  match_output <- matrix("", nids * 2, nloci)
  colnames(match_output) <- locinames
  number_loci <- rep(NA, nids)
  
  for (i in 1:nids) {
    nmatches_temp <- 0
    nloci_temp <- 0
    for (j in 1:nloci) {
      day0_alleles <- alleles0[i, (maxMOI * (j - 1) + 1):(maxMOI * (j - 1) + MOI0[i])]
      dayf_alleles <- allelesf[i, (maxMOI * (j - 1) + 1):(maxMOI * (j - 1) + MOIf[i])]
      
      if (sum(!is.na(day0_alleles)) > 0 && sum(!is.na(dayf_alleles)) > 0) {
        nloci_temp <- nloci_temp + 1
        allpossiblerecrud <- expand.grid(day0_alleles, dayf_alleles)
        allpossiblerecrud <- allpossiblerecrud[complete.cases(allpossiblerecrud), ]
        mindistance <- min(abs(allpossiblerecrud$Var1 - allpossiblerecrud$Var2))
        
        if (mindistance <= bin_sizes[j]) {
          nmatches_temp <- nmatches_temp + 1
          match_output[2 * (i - 1) + 1, j] <- "R" # Recrudescence
        } else {
          match_output[2 * (i - 1) + 1, j] <- "NI" # New Infection
        }
      } else {
        match_output[2 * (i - 1) + 1, j] <- "IND" # Indeterminate
      }
    }
    number_matches[i] <- nmatches_temp
    number_loci[i] <- nloci_temp
  }
  return(
    list(
      match_output_table = as.data.frame(match_output),
      number_matches = number_matches,
      number_loci = number_loci,
      nids = nids,
      nloci = nloci,
      locinames = locinames
    )
  )
}