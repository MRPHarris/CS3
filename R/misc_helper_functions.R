# Miscellaneous helper functions. Used mostly to extract model data.
# Many are adaptations of code from the staRdom package (Pucher et al., 2019).
# Some others are borrowed from my own package for convenience and to avoid dependency on that package, which I do not plan to list on CRAN: https://github.com/MRPHarris/eemUtils

#' Extract spectra at the peak wavelength position.
#'
#' @description Extract the spectra at the peak for a given component. Originally from https://github.com/MRPHarris/eemUtils.
#'
#' @param pfmodel a PARAFAC model object from staRdom::eem_parafaC()
#' @param component numeric. Target component.
#'
#' @noRd
#'
extrpf_peak_spectra_int <- function(pfmodel, component = 1){
  pfres <- pfmodel # Code adapted from staRdom package; naming will reflect this.
  names = TRUE
  c <- eempf_comp_mat(pfres) # get matrices of all the comps.
  if (is.null(names(c))) {
    names(c) <- paste0("model", seq(1:length(c)))
  }
  # Performing for only one model.
  c1 <- c # holdover code; lapply wrap removed here.
  nc1 <- length(c1)
  nc2 <- 0
  # for each component, pull out the data
  tab <- lapply(c1, function(c2) {
    nc2 <<- nc2 + 1
    c2 <- c2 %>% mutate(comps = nc1, comp = paste0(nc2))
  }) %>% bind_rows() %>%
    mutate_at(vars(ex, em, value, comp), as.numeric) # numeric vars
  # collate, organise, extract max spectra.
  comp_spectra <- tab %>%
    group_by(comp) %>%
    mutate(max_pos = which.max(value), max_em = em[max_pos], max_ex = ex[max_pos]) %>%
    mutate(exn = ifelse(em == max_em, ex, NA), emn = ifelse(ex == max_ex, em, NA)) %>%
    dplyr::filter(!is.na(emn) | !is.na(exn)) %>%
    ungroup() %>%
    mutate_at(vars(ex, em, value, comp), as.numeric) %>%
    dplyr::filter(comp == component)
  comp_spectra
}

#' Calculate and extract PARAFAC model residuals
#'
#' @description An expansion upon staRdom::eempf_residuals. Minor internal code changes and extra options added.
#'
#' @param pfmodel a PARAFAC model object from staRdom::eem_parafaC()
#' @param eem_list a list of eem objects compliant with the staRdom/eemR framework
#' @param select a character vector containing the names of the desired samples. See ?staRdom::eem_parafac
#' @param cores number of cores used for multi-threading. See ?staRdom::eem_parafac
#' @param denormalise is the PARAFAC model normalised? If so, use the eemlist to denormalise the loadings. Shouldn't have an influence on similarity scores.
#' @param extend_eemlist TRUE/FALSE to ensure equal EEM sizes before calculations. Uses staRdom::eem_extend2largest if differences are found.
#' @param verbose return various messages during operation.
#'
#' @noRd
#'
extrpf_residuals_int <- function(pfmodel, eem_list, select = NULL, cores = parallel::detectCores(logical = FALSE)-1,
                                denormalise = FALSE, extend_eemlist = TRUE, verbose = FALSE){
  # pfmodel <- norm2A(pfmodel)
  if (!is.null(select)) {
    eem_list <- eem_extract(eem_list, sample = select, keep = TRUE,
                            verbose = FALSE)
  }
  if(isTRUE(extend_eemlist)){
    message("Extending eemlist, if applicable")
    eem_list <- staRdom::eem_extend2largest(eem_list, interpolation = TRUE)
  }
  # revised normalisation
  if(isTRUE(denormalise)){
    pfmodel$A <- extrpf_loadings_denorm(pfmodel, eemlist = eem_list) %>%
      column_to_rownames(var = 'sample') %>%
      mutate_all(as.numeric)
  }
  if (!all(eem_names(eem_list) %in% rownames(pfmodel$A)) |
      length(eem_list) == 0) {
    pfmodel <- A_missing(eem_list, pfmodel, cores = cores)
  }
  what <- which(rownames(pfmodel$A) %in% (eem_list %>% eem_names()))
  pfmodel$A <- as.data.frame(pfmodel$A)[what, ]
  # building the residuals data. Multiple lapply layers.
  res_data <- lapply(pfmodel$A %>% rownames(), function(sample) {
    # First, main layer: iterate through each sample.
    comps <- lapply(pfmodel$A %>% colnames(), function(component) {
      # Another lapply! Prepare component matrices for this sample.
      s1 <- pfmodel$B[,component] %*% t(pfmodel$C[,component])
      s2 <- s1 * pfmodel$A[sample,component]
      s2
    })
    names(comps) <- pfmodel$C %>% colnames() # match excitation band names
    # Combine components to produce fit.
    fit <- comps %>% Reduce("+", .)
    eem <- eem_list[[which(eem_list %>% eem_names == sample)]]
    samp <- eem$x[eem$em %in% rownames(pfmodel$B), eem$ex %in%
                    rownames(pfmodel$C)]
    if(isTRUE(verbose)){
      message(sample)
    }
    res <- samp - fit
    comps <- lapply(pfmodel$A %>% colnames(), function(component) {
      comps[[component]] %>%
        data.frame() %>%
        mutate(type = component,em = rownames(pfmodel$B)) %>%
        gather(ex, value, -em, -type) %>% mutate(ex = substr(ex, 2, 4))
    }) %>% bind_rows()
    colnames(samp) <- rownames(pfmodel$C)
    samp <- samp %>% data.frame() %>% mutate(type = "sample",
                                             em = rownames(pfmodel$B)) %>% gather(ex, value, -em,
                                                                                  -type) %>% mutate(ex = substr(ex, 2, 4))
    rownames(res) <- rownames(pfmodel$B)
    colnames(res) <- rownames(pfmodel$C)
    res <- res %>% data.frame() %>% mutate(type = "residual",
                                           em = rownames(pfmodel$B)) %>% gather(ex, value, -em,
                                                                                -type) %>% mutate(ex = substr(ex, 2, 4))
    res2 <- bind_rows(list(comps, samp, res)) %>% mutate(sample = sample)
  }) %>% bind_rows()
}


#' Extract spectra from an EEM at a given Ex/Em wavelength.
#'
#' @description 'Slice' an EEM at a given Ex/Em coordinate pair, and extract the spectra at that location. Originally from https://github.com/MRPHarris/eemUtils.
#'
#' @param eem an eem object compliant with the staRdom/eemR framework
#' @param ex an excitation wavelength value, in nm
#' @param em an emission wavelength value, in nm
#'
#' @noRd
#'
slice_eem_int <- function(eem, ex, em){
  # input checks
  if (is.data.frame(eem)) {
  } else if (class(eem) == "eem") {
    eem_df <- as.data.frame(eem, gather = FALSE)
  } else {
    stop("Please pass the function an object of class 'eem'")
  }
  # Extract emission values slice at given ex wavelength
  if(ex %in% colnames(eem_df)){
    em_slice <- data.frame(matrix(NA, nrow = nrow(eem_df), ncol = 2))
    em_slice[, 1] <- as.numeric(rownames(eem_df))
    em_slice[, 2] <- as.numeric(as.matrix(eem_df[, which(colnames(eem_df) == ex)]))
    colnames(em_slice) <- c("emission", "intensity")
    em_slice <- pivot_longer(em_slice, cols = emission, values_to = "wavelength")
  } else {
    ex <- binary_search_nearest(data = colnames(eem_df), value = ex)
    em_slice <- data.frame(matrix(NA, nrow = nrow(eem_df), ncol = 2))
    em_slice[, 1] <- as.numeric(rownames(eem_df))
    em_slice[, 2] <- as.numeric(as.matrix(eem_df[, which(colnames(eem_df) == ex)]))
    colnames(em_slice) <- c("emission", "intensity")
    em_slice <- pivot_longer(em_slice, cols = emission, values_to = "wavelength")
  }
  # Extract excitation values slice at given em wavelength
  if(em %in% rownames(eem_df)){
    # Continue as normally
    ex_slice <- data.frame(matrix(NA, nrow = ncol(eem_df), ncol = 2))
    ex_slice[, 1] <- as.numeric(colnames(eem_df))
    ex_slice[, 2] <- as.numeric(as.matrix(t(eem_df[which(rownames(eem_df) == em), ])))
    colnames(ex_slice) <- c("excitation", "intensity")
    ex_slice <- pivot_longer(ex_slice, cols = excitation, values_to = "wavelength")
  } else {
    # Oh no! Find nearest em value with binary search via data.table.
    em <- binary_search_nearest(data = rownames(eem_df), value = em)
    ex_slice <- data.frame(matrix(NA, nrow = ncol(eem_df), ncol = 2))
    ex_slice[, 1] <- as.numeric(colnames(eem_df))
    ex_slice[, 2] <- as.numeric(as.matrix(t(eem_df[which(rownames(eem_df) == em), ])))
    colnames(ex_slice) <- c("excitation", "intensity")
    ex_slice <- pivot_longer(ex_slice, cols = excitation, values_to = "wavelength")
  }
  # Combine
  slices <- rbind(em_slice, ex_slice)
  slices %>% mutate_at(vars(intensity, wavelength), as.numeric)
  slices
}
