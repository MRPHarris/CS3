# Miscellaneous functions. Used mostly to extract model data.
# Many are adaptations of code from the staRdom package (Pucher et al., 2019).
# Some others are borrowed from my own package for convenience and to avoid dependency on that package, which I do not plan to list on CRAN: https://github.com/MRPHarris/eemUtils

#' Get SSC along with alpha and beta penalty terms
#'
#' @description The shift- and shape sensitive congruence (SSC) was developed by Wunsch et al., 2019 as an improvement
#'      upon the TCC metric. It incorporates two penalty terms, alpha and beta, to account for differences in the wavelength
#'      peak position and area. This function adds these terms to the data frame returned by staRdom::ssc(). Code from https://github.com/MRPHarris/eemUtils.
#'
#' @param mat1 a matrix
#' @param mat2 a matrix
#' @param tcc TRUE/FALSE to return only TCC value instead of SSC, alpha and beta.
#'
#' @noRd
#'
ssc_more_int <- function (mat1, mat2, tcc = FALSE) {
  if (any(is.null(mat1), is.na(mat1), is.null(mat2), is.na(mat2))) {
    a <- NA
  } else {
    a <- lapply(1:ncol(mat1), function(nc) {
      col1 <- mat1[, nc]
      apply(mat2, 2, function(col2) {
        tcc_cal <- sum(col1 * col2)/sqrt(sum(col1^2) *
                                           sum(col2^2))
        if (!tcc) {
          wl <- as.numeric(names(col1))
          if (any(is.na(wl)) | pracma::isempty(wl)) {
            stop("SSCs cannot be calculated. Please add wavelengths as rownames of the matrices!")
          }
          alpha <- abs((wl[which.max(col1)] - wl[which.max(col2)])/diff(range(wl)))
          beta <- abs((sum(col1/max(col1)) - sum(col2/max(col2)))/diff(range(wl)))
          ssc <- tcc_cal - alpha - beta
          ssc <- c(ssc, alpha, beta)
          # rownames(a) <- (c("ssc","alpha","beta"))
        } else {
          tcc_cal
        }
      })
    }) %>% setNames(colnames(mat1)) %>% do.call(rbind, .)
  }
  attr(a, "method") <- ifelse(tcc, "TCC", "SSC")
  if(!isTRUE(tcc)){
    rownames (a) <- c("ssc","alpha","beta")
  }
  a
}

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
    pfmodel$A <- extrpf_loadings_denorm_int(pfmodel, eemlist = eem_list) %>%
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

#' Multiply the loadings values of PARAFAC components by normalisation factors
#'
#' @description A PARAFAC model fed with normalised eems outputs by default loadings that are less useful for direct
#'       sample-to-sample comparison. This can be remedied by multiplying the loadings by each eem fmax value (i.e. the
#'       normalisation factors used to normalise the eems originally). Originally from https://github.com/MRPHarris/eemUtils.
#'
#' @param pfmodel a PARAFAC model object, typically an individual ouptut from staRdom::eem_parafaC()
#' @param eemlist a list of eem objects compliant with the staRdom/eemR framework
#' @param type short or long format. Short by default. Long format data works better for grouping in ggplot
#'
#' @noRd
#'
extrpf_loadings_denorm_int <- function(pfmodel, eemlist, type = "short"){
  maxvals <- eemlist_fmax_values(eemlist)
  fm <- extrpf_loadings(pfmodel)
  loadings_frame <- fm[,2:ncol(fm)]
  newframe <- apply(loadings_frame, 2, function(col) {
    col * maxvals
  }) %>% data.frame()
  newframe <- newframe %>% `colnames<-`(c(paste0("Comp.",
                                                 seq(1, ncol(newframe), 1)))) %>% `rownames<-`(unlist(lapply(eemlist,
                                                                                                             "[[", "sample"))) %>% rownames_to_column(var = "sample")
  if (type == "short") {
    newframe
  }
  else if (type == "long") {
    newframe <- newframe %>% pivot_longer(cols = c(2:ncol(newframe)))
  }
  else {
    stop("Unknown 'type'. Please input either 'short' or 'long'")
  }
}

