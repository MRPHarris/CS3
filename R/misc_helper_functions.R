# Miscellaneous functions. Used mostly to extract model data.
# Many are adaptations of code from the staRdom package (Pucher et al., 2019).
# Some others are borrowed from my own package for convenience and to avoid dependency on that package, which I do not plan to list on CRAN: https://github.com/MRPHarris/eemUtils

#' Create an EEM object representing the Long-Term Method Detection Limit LT-MDL
#'
#' @description The LT-MDL quantifies a given fluorometer's ability to detect signals
#'       at specific wavelengths by combining the response of a set of ultrapure water blanks
#'       run using the same parameters as the method of interest. Two options are available for
#'       calculating the LT-MDL: the average of a set of ultrapure water blanks plus 3* the standard
#'       deviation of the blanks (following Hansen et al., 2018), or simply 3* the standard deviation of
#'       the blanks (following Thomsen et al., 2003).
#'
#' @param blank_eemlist an eemlist compliant with the eem/eemR/staRdom framework, containing a set of ultrapure water blank EEMs created using identical measurement parameters.
#' @param remove_gamma_spikes uses forecast::tsclean to treat each emission scan in the LT-MDL as a time series and remove outliers. This gets rid of single-point scatter artefacts caused by gamma spikes.
#' @param excise_scatter four logical inputs, e.g. c(TRUE,TRUE,TRUE,TRUE). Passed to staRdom::eem_rem_scat. See ?staRdom::eem_rem_scat for more information.
#' @param scatter_widths four numeric inputs, specifying the width, in nm, of scatter excissions. See ?staRdom::eem_rem_scat.
#' @param method character vector, one of either 'Hansen' or 'Thomsen' corresponding to the two supproted methods for LT-MDL calculation.
#'
#' @importFrom staRdom eem_rem_scat
#' @importFrom parallel detectCores
#' @importFrom staRdom eem_interp
#' @importFrom magrittr %>%
#'
#' @export
#'
create_MDL_eem <- function(blank_eemlist,
                           remove_gamma_spikes = TRUE,
                           excise_scatter = c(TRUE,TRUE,TRUE,TRUE),
                           scatter_widths = c(12,12,10,10),
                           method = "Hansen"){
  ## Only Rayleigh excision performed on raw blank data.
  if((excise_scatter[3] == TRUE) || (excise_scatter[4] == TRUE)){
    message("Removing Rayleigh scatter of specified orders from blank eems.")
    scatter <- c(FALSE, FALSE, excise_scatter[3], excise_scatter[4])
    scatter_width <- c(0,0,scatter_widths[3],scatter_widths[4])
    ## Apply scatter removal
    mq_eems_masked <- blank_eemlist %>%
      eem_rem_scat(remove_scatter = scatter, remove_scatter_width = scatter_width) %>%
      eem_rayleigh_zero(order = 1, width = scatter_widths[3])
    ## Interpolate
    cores <- detectCores(logical = FALSE)
    mq_eems_interp <- eem_interp(mq_eems_masked, type = 1, extend = FALSE, cores = cores)
  } else {
    mq_eems_interp <- blank_eemlist
  }
  ### Create average of mq eems
  mq_average <- unlist(eemlist_average_int(mq_eems_interp), recursive = FALSE) %>% 'class<-'('eem')
  mq_average$sample <- "MQav"
  ### Generate 3*the standard deviation.
  mq_sd3 <- eemlist_sd_int(eemlist = mq_eems_interp, mult = 3)
  mq_sd3$sample = c("MQsd3")
  ## Generate the LT-MDL EEM object after Hansen
  if(method == "Hansen"){
    eemlist_avpsd <- list(mq_sd3,mq_average)
    class(eemlist_avpsd) <- c('eemlist')
    eem_avpsd <- eemlist_sum(eemlist_avpsd)
    eem_avpsd$sample <- c("LT-MDL for method m1")
  } else if(method == "Thomsen"){
    eem_avpsd <- mq_sd3
  } else {
    stop("Please supply 'method' as either 'Hansen' or 'Thomsen'")
  }
  eemlist_avpsd <- list(eem_avpsd)
  class(eemlist_avpsd) = c('eemlist')
  ## Additional scatter removal steps from LT_MDL.
  if((excise_scatter[3] == TRUE) || (excise_scatter[4] == TRUE)){
    message("Ensuring Rayleigh scatter of specified orders removed from MDL EEM.")
    scatter <- c(FALSE, FALSE, excise_scatter[3], excise_scatter[4])
    scatter_width <- c(0,0,scatter_widths[3],scatter_widths[4])
    ## Apply scatter removal
    eemlist_avpsd_masked <- eemlist_avpsd %>%
      eem_rem_scat(remove_scatter = scatter, remove_scatter_width = scatter_width) %>%
      eem_rayleigh_zero(order = 1, width = scatter_widths[3])
    ## Interpolate
    cores <- detectCores(logical = FALSE)
    eemlist_avpsd_interp <- eem_interp(eemlist_avpsd_masked, type = 1, extend = FALSE, cores = cores)
  } else {
    eemlist_avpsd_interp <- eemlist_avpsd
  }
  ## Raman removal from combined EEM
  if((excise_scatter[1] == TRUE) || (excise_scatter[2] == TRUE)){
    message("Removing raman scatter from MDL EEM.")
    scatter <- c(excise_scatter[1], excise_scatter[2], FALSE, FALSE)
    scatter_width <- c(scatter_widths[1],scatter_widths[2],0,0)
    # ## Apply scatter removal
    eemlist_avpsd_masked2 <- eemlist_avpsd_interp %>%
      eem_rem_scat(remove_scatter = scatter, remove_scatter_width = scatter_width)
    ## Interpolate
    cores <- detectCores(logical = FALSE)
    eemlist_avpsd_interp2 <- eem_interp(eemlist_avpsd_masked2, type = 1, extend = FALSE, cores = cores)
  } else {
    eemlist_avpsd_interp2 <- eemlist_avpsd_interp
  }
  ## Collate
  LT_MDL_m1 <- unlist(eemlist_avpsd_interp2, recursive = FALSE)
  class(LT_MDL_m1) <- c('eem')
  ## Gamma ray spike removal
  if(isTRUE(remove_gamma_spikes)){
    MDL_lst <- list(LT_MDL_m1)
    class(MDL_lst) <- 'eemlist'
    MDL_lst_dn <- eemlist_sp_denoise_int(MDL_lst)
    LT_MDL_m1 <- unlist(MDL_lst_dn, recursive = FALSE)
    class(LT_MDL_m1) <- c('eem')
    LT_MDL_m1
  } else {
    LT_MDL_m1
  }
}




#' Get SSC along with alpha and beta penalty terms
#'
#' @description The shift- and shape sensitive congruence (SSC) was developed by Wunsch et al., 2019 as an improvement
#'      upon the TCC metric. It incorporates two penalty terms, alpha and beta, to account for differences in the wavelength
#'      peak position and area. This function adds these terms to the data frame returned by staRdom::ssc(). Code from https://github.com/MRPHarris/eemUtils.
#'
#' @param mat1 a matrix
#' @param mat2 a matrix
#' @param tcc TRUE/FALSE to return only TCC value instead of SSC, alpha and beta.
#' @param alpha_start NULL or a value to set the starting point for peak detection. If a numeric value is supplied, in nanometres, this will be set as the start of the search for the peak position in each spectra. If the value is NULL or NA, the operation is ignored.
#' @param abs_beta TRUE/FALSE to use the absolute value of both mats when calculating the beta term. This prevents negative values in either mat from reducing the penalty score.
#'
#' @noRd
#'
ssc_more_int <- function (mat1, mat2, tcc = FALSE, alpha_start = NULL, abs_beta = TRUE) {
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
          if(!is.null(alpha_start) && !is.na(alpha_start)){
            if(!is.numeric(alpha_start)){
              stop("The trough wavelength value used for the modified alpha calculation must be numeric.")
            }
            trough_ind <- as.numeric(which(wl == alpha_start))
            wl_trim <- wl[trough_ind:length(wl)]
            alpha <- abs((wl_trim[which.max(col1[trough_ind:length(col1)])] - wl_trim[which.max(col2[trough_ind:length(col2)])])/diff(range(wl)))
          } else {
            # Normal alpha calculation
            alpha <- abs((wl[which.max(col1)] - wl[which.max(col2)])/diff(range(wl)))
          }
          if(isTRUE(abs_beta)){
            col2abs <- abs(col2)
            col1abs <- abs(col1)
            beta <- abs((sum(col1abs/max(col1abs)) - sum(col2abs/max(col2abs)))/diff(range(wl)))
          } else {
            beta <- abs((sum(col1/max(col1)) - sum(col2/max(col2)))/diff(range(wl)))
          }
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
#' @importFrom dplyr bind_rows
#' @importFrom dplyr mutate_at
#' @importFrom dplyr mutate
#' @importFrom dplyr vars
#' @importFrom dplyr group_by
#' @importFrom dplyr ungroup
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
#' @importFrom eemR eem_extract
#' @importFrom tibble column_to_rownames
#' @importFrom dplyr mutate_all
#' @importFrom eemR eem_names
#' @importFrom staRdom A_missing
#' @importFrom dplyr mutate
#' @importFrom tidyr gather
#' @importFrom dplyr bind_rows
#'
#' @noRd
#'
extrpf_residuals_int <- function(pfmodel, eem_list, select = NULL, cores = parallel::detectCores(logical = FALSE)-1,
                                denormalise = FALSE, extend_eemlist = TRUE, verbose = FALSE, force_names = TRUE){
  ncomps = ncol(pfmodel$A)
  colnames_pre <- colnames(pfmodel$A)
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
    if(isTRUE(force_names)){
      rownames(pfmodel$A) <- eem_names(eem_list)
    } else {
      pfmodel <- A_missing(eem_list, pfmodel, cores = cores)
    }
  }
  what <- which(rownames(pfmodel$A) %in% (eem_list %>% eem_names()))
  pfmodel$A <- as.matrix((pfmodel$A)[what, ]) %>% 'colnames<-'(colnames_pre)
  # building the residuals data. Multiple lapply layers.
  res_data <- lapply(pfmodel$A %>% rownames(), function(sample) {
    # Lapply over samples
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
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate_at
#' @importFrom dplyr vars
#' @importFrom stats setNames
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
    ex <- binary_search_nearest_int(data = colnames(eem_df), value = ex)
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
    em <- binary_search_nearest_int(data = rownames(eem_df), value = em)
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
#' @importFrom tidyr pivot_longer
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

#' Ensure that wavelengths where Ex > Em are set to 0 when removing 1st-order Rayleigh scattering.
#'
#' @description A modification to eemR's scatter removal function that ensures wavelengths proximal to
#'      the 1st order Rayleigh scatter line where Excitation > Emission are set to 0, ensuring there is
#'      no impossible fluorescence upon subsequent interpolation over scattering areas.
#'
#' @param eem an eem object compliant with the eem/eemR/staRdom framework
#' @param order integer, either 1 for first order or 2 for second order
#' @param width width in nm for the scatter excission
#' @param corr0 logical; TRUE/FALSE; if TRUE, ensures wavelengths where Ex > Em are set to an intensity of 0.
#'
#' @noRd
eem_rayleigh_zero <- function (eem, order = 1, width = 10, corr0 = TRUE){
  if (class(eem) == "eemlist") {
    res <- lapply(eem, eem_rayleigh_zero,
                  order = order, width = width)
    class(res) <- class(eem)
    return(res)
  }
  x <- eem$x
  em <- eem$em
  ex <- eem$ex
  ind1 <- mapply(function(x) em <= x, order * ex)
  ind2 <- mapply(function(x) em <= x, order * ex + width)
  ind3 <- ifelse(ind1 + ind2 == 1, NA, 1)
  x <- x * ind3
  if(isTRUE(corr0)){
    ind4 <- mapply(function(x) em <= x, ex)
    ind5 <- ifelse(ind4 == 1, 0, 1)
    x <- x * ind5
  }
  res <- eem
  res$x <- x
  attributes(res) <- attributes(eem)
  attr(res, "is_scatter_corrected") <- TRUE
  class(res) <- class(eem)
  return(res)
}

#' A single-point denoiser for long-term method detection and quantification limit EEM objects.
#'
#' @description EEMs are often punctuated by isolated, single points of scatter. These are potentially
#'        caused by a wide range of factors, including gamma-ray spikes. They can severely hamper normalisation
#'        attempts, and throw off method detection limit calculations. This function provides a simple denoising
#'        solution using the tsclean function from the forecast package. Data supplied to this function should not
#'        contain missing values. Borrowed from the package eemUtils.
#'
#' @param eemlist A list of EEMs in a format compliant with eemR/staRdom.
#' @param interp TRUE/FALSE to replace outliers with linear interpolated values instead of NA.
#'
#' @importFrom forecast tsclean
#'
#' @noRd
#'
eemlist_sp_denoise_int <- function(eemlist, interp = TRUE){
  denoise_eemlist <- lapply(eemlist, function(e){
    e_df <- as.data.frame(e, gather = FALSE)
    clean_df <- lapply(e_df, function(ex_wl){
      wl_clean <- tsclean(ex_wl, replace.missing = interp)
    }) %>%
      data.frame() %>%
      'rownames<-'(rownames(e_df)) %>%
      'colnames<-'(colnames(e_df))
    clean_df <- eemdf_to_eem(eemdf = clean_df,
                             file = e$file,
                             sample = e$sample,
                             location = e$location,
                             gathered = FALSE)
    clean_df
  })
  return(denoise_eemlist)
}

#' EEM object constructor
#'
#' @description borrowed from eemUtils, adapting staRdom's eem_csv importer to turn
#'        dataframes into eems.
#'
#' @param eemdf the dataframe to be coerced to an EEM object.
#' @param file filename of the EEM, if applicable.
#' @param sample the samplename of the EEM, if applicable.
#' @param location the location of the EEM file, if applicable.
#' @param gathered TRUE/FALSE is the eemdf in a short (not gathered; FALSE) or a long (gathered; TRUE) format?
#'
#' @importFrom magrittr %>%
#' @importFrom tidyr spread
#' @importFrom dplyr select
#'
#' @noRd
#'
eemdf_to_eem_int <- function(eemdf,
                         file = NULL,
                         sample = NULL,
                         location = NULL,
                         gathered = FALSE){
  # code adapted from staRdom's .eem_csv importer.
  x <- eemdf
  if(!isTRUE(gathered)){
    # The eem is in a short, non-gathered format.
    ex <- colnames(x)[] %>% as.numeric()
    em <- rownames(x) %>% as.numeric()
    x <- x[,] %>% as.matrix() %>% unname()
    x <- x[!is.na(em),!is.na(ex)]
    ex <- ex[!is.na(ex)]
    em <- em[!is.na(em)]
    l <- list(
      file = file,
      sample = sample,
      x = x,
      ex = ex,
      em = em,
      location = location
    )
    class(l) <- "eem"
    return(l)
  } else {
    # The eem is in a long or 'gathered' format.
    gath_df <- eemdf
    gath_df$value <- as.numeric(as.character(gath_df$value)) # factor management
    if("sample" %in% colnames(gath_df)){
      gath_df <- subset(gath_df, select = -c(sample))
    }
    if(colnames(gath_df)[3] == "z"){
      gath_df_short <- spread(data = gath_df, key = "ex", value = "z")
    } else {
      gath_df_short <- spread(data = gath_df, key = "ex", value = "value")
    }
    rnames <- as.matrix(gath_df_short[,1])
    rownames(gath_df_short) <- as.numeric(rnames)
    gath_df_short <- select(gath_df_short, -c(1))
    eem <- eemdf_to_eem(eemdf = gath_df_short,
                        file = file,
                        sample = sample,
                        location = location,
                        gathered = FALSE)
    return(eem)
  }
}

#' Generate standard deviations of an eemlist
#'
#' @description Produces an EEM object containing ex/em-pair standard deviations for all eems within the supplied
#'       eemlist. Borrowed from eemUtils.
#'
#' @param eemlist A list of EEMs in a format compliant with eemR/staRdom.
#' @param mult Numeric; a simple numeric multiplier applied to the resulting ex-em pair sd values.
#'
#' @noRd
#'
eemlist_sd_int <- function(eemlist, mult = 1){
  cols <- eemlist[[1]]$ex
  rows <- eemlist[[1]]$em
  # Coerce eemlist to sd-operable matrices
  eemlist_matrices <- lapply(eemlist, function(x){
    eem_df <- as.data.frame(x, gather = FALSE)
    eem_df <- as.matrix(eem_df)
  })
  # Produce standard deviations
  sd <- apply(simplify2array(eemlist_matrices), 1:2, sd)
  # Multiplication handling
  if(mult != 1){
    if(!is.numeric(mult)){
      stop("'mult' must be numeric.")
    } else {
      sd <- sd * mult
    }
  }
  sd <- data.frame(sd)
  colnames(sd) <- cols
  rownames(sd) <- rows
  sd_eem <- eemdf_to_eem_int(sd, sample = "eemlist_sd")
}

#' Average a set of EEMs.
#'
#' @description Average the EEMs within an eemlist. Borrowed from eemUtils.
#'
#' @param eemlist A list of EEMs, compliant with the eemR/staRdom framework.
#'
#' @export
#'
eemlist_average_int <- function(eemlist){
  if(length(eemlist) == 1){
    message("1 EEM passed to average_eems() for averaging. Returning unchanged.")
    new_eemlist <- eemlist
    return(new_eemlist)
  } else if(length(eemlist) > 1){
    ungathered_list <- vector(mode = "list",length = length(eemlist))
    for(i in seq_along(eemlist)){
      eem_it <- eemlist[[i]]
      file_it <- eem_it[['file']]
      sample_it <- eem_it[['sample']]
      location_it <- eem_it[['location']]
      eem_ungathered <- as.data.frame(eemlist[[i]], gather = FALSE)       # extract EEM, don't gather
      ungathered_list[[i]] <- eem_ungathered
    }
    # Average it
    averaged <- Reduce("+",ungathered_list)/length(ungathered_list)
    # Now back to eem
    averaged_eem <- eemdf_to_eem_int(averaged,
                                 file = "",
                                 sample = "averaged_eem",
                                 location = "")
    new_eemlist <- vector(mode = "list",length = 1)
    class(new_eemlist) <- "eemlist"
    message(paste0(length(eemlist)," EEMs passed to average_eems() for averaging."))
    new_eemlist[[1]] <- averaged_eem
    return(new_eemlist)
  }
}
