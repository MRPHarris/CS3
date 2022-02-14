# Functions used in CS3 for correcting spectra.

#' Interpolate spectra to 1nm wavelength bands.
#'
#' @description A simple interpolating function that ensures spectra are comprised of uniform 1nm bandwidths. Adapted from the Murphy (2011) MATLAB code.
#'
#' @param matrix a data.matrix with wavelengths as rownames
#' @param mat_out TRUE/FALSE to output a matrix of identical format to the input (i.e. a data.matrix with the wavelengths as the rownames)
#'
#' @importFrom stats approx
#' @importFrom tibble column_to_rownames
#'
#' @noRd
#'
interpolate_matrix_1nm <- function(mat, mat_out = FALSE){
  if(is.null(rownames(mat))){
    stop("Please add wavelengths as the rownames of the input matrix")
  }
  values <- as.numeric(mat)
  wavelengths <- as.numeric(rownames(mat))
  min_wl_ind <- as.numeric(which(wavelengths == min(wavelengths)))
  max_wl_ind <- as.numeric(which(wavelengths == max(wavelengths)))
  wl_for_interp <- seq(min(wavelengths),max(wavelengths),1)
  new_values <- approx(x = wavelengths, # linear interpolation using approx()
                       y = values,
                       xout = wl_for_interp,
                       method = "linear")
  new_values_df <- data.frame(matrix(t(unlist(new_values)), nrow = length(new_values$x), ncol = 2))
  colnames(new_values_df) <- c("wavelength","value")
  if(!isTRUE(mat_out)){
    new_values_df
  } else {
    # Return matrix for other purposes.
    new_values_mat <- new_values_df %>%
      column_to_rownames('wavelength') %>%
      data.matrix
    new_values_mat
  }
}

#' Remove the contribution of other components to sample fluorescence in the target spectra.
#'
#' @description Removes the contribution of other modelled components from a given set of ex/em spectra. Ensures that the target component is being compared with its constituent data (i.e. sample + residuals)
#'
#' @param grob residuals object containing ex or em data from the target component peak position.
#' @param sample_char a character vector matching the name of the EEM to be targeted
#' @param comp numeric, target component
#' @param type one of either 'ex' or 'em'
#' @param output_mat TRUE/FALSE to output either a simple matrix object used by other functions, or a longer form output.
#' @param neg_to_0 TRUE/FALSE to set all values of negative fluoresence to 0 after subtraction.
#' @param radd_nontarget TRUE/FALSE to add residuals and target component spectra to the output. Only applies if output_mat = FALSE.
#'
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr mutate
#' @importFrom dplyr across
#' @importFrom dplyr select
#' @importFrom dplyr bind_rows
#'
#' @noRd
#'
comp_correct_spectra <- function(grob = NULL, sample_char, comp, type = "ex",
                                 output_mat = TRUE, neg_to_0 = TRUE, add_nontarget = TRUE){
  target_comp <- comp
  target_comp_string <- paste0("Comp.", target_comp)
  if(!is.character(sample_char)){
    stop("Please supply a character vector for 'sample_char' matching that of the EEM.")
  }
  if(is.null(grob)){
    stop("Please supply a grob to correct spectra.")
  }
  # Pull the sample out of the grob.
  grob_sample <- grob[which(grob$sample == sample_char),]
  # Get all the types: components, residuals, sample
  types <- unique(grob_sample$type)
  # Which comps are *not* the target?
  non_target_comps <- types[which(types != target_comp_string & types != 'residual' & types != 'sample')]
  # Extract spectra for non-target components.
  non_target_spectra <- grob_sample[(apply(as.data.frame(grob_sample$type), 2, function(r) which(r %in% non_target_comps))),]
  # Sum them
  sum_spectra <- non_target_spectra %>%
    pivot_wider(names_from = type, values_from = value) %>%
    replace(is.na(.), 0) %>%
    mutate(combined.spectra = rowSums(across(.cols = 4:(4+length(non_target_comps)-1)))) %>%
    select(-(4:(4+length(non_target_comps)-1)))
  # Get the sample spectra
  sample_spectra <- grob_sample[(apply(as.data.frame(grob_sample$type), 2, function(r) which(r %in% 'sample'))),]
  # Perform the subtraction.
  corrected_spectra <- sample_spectra
  corrected_spectra$value <- sample_spectra$value - sum_spectra$combined.spectra
  if(isTRUE(output_mat)){
    # return a matrix that can replace the relevant matrix.
    mat <- data.matrix(corrected_spectra$value)
    if(isTRUE(neg_to_0)){
      mat[mat<0] <- 0
    }
    if(type == "ex"){
      rownames(mat) <- corrected_spectra$ex
    } else if(type =="em"){
      rownames(mat) <- corrected_spectra$em
    } else {
      stop("Please supply a type of 'ex' or 'em'.")
    }
    mat
  } else if(!isTRUE(output_mat)){
    if(isTRUE(add_nontarget)){
      residuals <- grob_sample[(apply(as.data.frame(grob_sample$type), 2, function(r) which(r %in% 'residual'))),]
      targetcomp <- grob_sample[(apply(as.data.frame(grob_sample$type), 2, function(r) which(r %in% target_comp_string))),]
      corrected_spectra_plus <- bind_rows(corrected_spectra, residuals, targetcomp)
      corrected_spectra_plus
    } else {
      corrected_spectra
    }
  }
}

#' Perform spectral smoothing with a Savitzky-Golay filter.
#'
#' @description A wrapper for signal::sgolayfilt that incorporates negative value handling.
#'
#' @param mat an input data.matrix
#' @param neg_to_0 TRUE/FALSE to set negative values to 0. Defaults to true
#' @param p filter order. See ?signal::sgolayfilt
#' @param n window filter size (must be odd). See ?signal::sgolayfilt
#' @param ... arguments passed on to signal::sgolayfilt
#'
#' @importFrom signal sgolayfilt
#'
#' @noRd
#'
sg_smooth <- function(mat, neg_to_0 = TRUE , p = 2, n = 21,...){
  # Apply Savitzky-Golay smoother.
  mat_smth <- sgolayfilt(mat, p = p, n = n,...) %>%
    data.matrix() %>%
    'rownames<-'(c(rownames(mat)))
  if(isTRUE(neg_to_0)){
    mat_smth[1][mat_smth[1] < 0] <- 0
  }
  mat_smth
}

#' Trim incomplete lower wavelength peaks from spectral data
#'
#' @description Apply a gradient detection method to identify and trim incomplete peaks towards higher wavenumbers/lower wavelengths. Utilises code from Murphy, K. R. (2011). notes A Note on Determining the Extent of the Water Raman Peak in Fluorescence Spectroscopy. Applied Spectroscopy, 65(2), 233–236. https://doi.org/10.1366/10-06136.
#'
#' @param mat1 a single spectra in the form of a data.matrix
#' @param mat2 a single spectra in the form of a data.matrix
#' @param tolerance numeric decimal value representing percentage tolerance. 0.05 (5%) by default.
#' @param verbose TRUE/FALSE to return various messages during operation.
#'
#' @noRd
#'
extract_complete_peak <- function(mat1, mat2, tolerance = 0.05, verbose = FALSE){
  # list matrices
  matrices <- list(mat1,mat2) %>%
    'names<-'(c("mat1","mat2"))
  matrices_interpolated <- lapply(matrices, interpolate_matrix_1nm) %>%
    'names<-'(c("mat1","mat2"))
  # get gradient segments
  gradient_values <- lapply(matrices_interpolated, get_gradient_vals, tolerance = tolerance) %>%
    'names<-'(c("mat1","mat2"))
  # get number of segments per segment object
  # In a typical, dual-peak scenario the EEM will feature two negative segments and one positive.
  segments <- lapply(gradient_values, nsegments)
  if(length(which(as.numeric(unlist(segments)) > 2)) > 0 & isTRUE(verbose)){
    message("More than 2 negative gradient sections. Either there's more than two peaks present, or the data is noisy (and you should increase the tolerance value).")
  }
  ## If the above message doesn't trip, proceed.
  # primary peak locations: max intensity value in the given spectra.
  primary_peak_locations <- lapply(matrices, function(mat){
    rownames(mat)[which.max(mat)]
  }) %>%
    `names<-`(c("mat1", "mat2"))
  # Secondary peak locations. Not needed explicitly, but useful.
  counter <- 0
  secondary_peak_locations <- lapply(matrices, function(mat){
    counter <<- counter + 1
    val <- binary_search_nearest_int(data = rownames(mat), value = max(gradient_values[[counter]][['pos']], na.rm = TRUE))
    val
  })
  ## trough locations. Looped; lapply wouldn't crack it.
  trough_wavelengths <- vector("list", length = 2)
  for(v in seq_along(trough_wavelengths)){
    mat <- matrices[[v]]
    neg_grads_it <- gradient_values[[v]][['neg']][!is.na(gradient_values[[v]][['neg']])]
    trough <- neg_grads_it[which(diff(neg_grads_it) != 1)]
    trough_location <- binary_search_nearest_int(data = rownames(mat), value = trough)
    trough_wavelengths[[v]] <- trough_location
  }
  # Check that the trough exists. I.e., there is a number for both entries.
  troughs <- unlist(trough_wavelengths)
  if(any(is.na(troughs))){
    missing <- which(is.na(troughs))
    if(length(missing) == 1){
      # There's one trough. In this case, set the NA to 0 and carry on.
      troughs[missing] <- 0 # set to 0; the max find below will ignore it.
      ## find the greater of the two wavelengths.
      starting_wavelength <- max(troughs)
      secondary_spectra <- lapply(matrices, function(mat){
        mat_new <- data.matrix(mat[which(rownames(mat) == starting_wavelength):nrow(mat),])
      })
      secondary_spectra
    } else if(length(missing) == 2){
      # No troughs.
      secondary_spectra <- list(mat1,mat2)
      secondary_spectra
    }
  } else {
    # Two troughs, get the greater one.
    ## find the greater of the two wavelengths.
    starting_wavelength <- max(troughs)
    # New spectra
    secondary_spectra <- lapply(matrices, function(mat){
      mat_new <- data.matrix(mat[which(rownames(mat) == starting_wavelength):nrow(mat),])
    })
    # Return
    secondary_spectra
  }
}



