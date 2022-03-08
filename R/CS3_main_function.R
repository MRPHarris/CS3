# Here lies the main function in the CS3 package.
# The function compares PARAFAC spectra with underlying sample spectra at the target peak coordinates using spectral similarity metrics (TCC and/or SSC). Various corrections may be applied to the data.
# The function is written with the goal of being understandable - it is certainly not optimised for speed/efficiency.

#' Compare PARAFAC and EEM spectra with spectral similarity metrics
#'
#' @description Extract the spectra from a PARAFAC model (B and C mode loadings) and compare them with underlying
#'       sample data. Various options exist for corrections, including adjusting excitation spectra so that SSC might
#'       be applied rather than just TCC.
#'
#' @param pfmodel a PARAFAC model object. An output from staRdom::eem_parafac()
#' @param eemlist a group of EEMs compliant with the staRdom/EEM/eemR package framework
#' @param comp numeric, singular - the target component for comparisons
#' @param tcc TRUE/FALSE to extract only TCC rather than SSC
#' @param terms TRUE/FALSE to extract the alpha and beta penalty term values alongside SSC
#' @param modified_metrics TRUE/FALSE to add the modified (ex + em combined) TCC and SSC metrics after Parr et al., 2014.
#' @param spectral_correct either NULL, 'all', 'ex', or 'em'. Subtract the loadings of other components from the raw sample spectra of the type specified
#' @param interp_1nm either NULL, 'all', 'ex', or 'em'.  Interpolate spectra to 1nm bandwidths. Recommended. Applied after spectral correction
#' @param smooth_sg either NULL, 'all', 'ex', or 'em'. Applies a Savitzky-Golay filter to the data prior to metric calculation. Only recommended if interpolation is performed. Uses signal::sgolay(). Default values are a 2nd order polynomial, n = 21 for emission spectra and n = 11 for excitation spectra
#' @param complete_peak either NULL, 'all', 'ex', or 'em'. Use a gradient detection method to identify incomplete peaks in spectra, and prevent them from interfering with the SSC peak position calculation. Set options with the parameter 'SSC_trim_method'.
#' @param SSC_trim_method if !is.null(complete_peak), choose one of either 'mod_alpha' to alter the peak position penalty term, or 'trim_spectra' to remove the incomplete peak entirely
#' @param verbose TRUE/FALSE to return various messages during the function's opperation. Useful for error checking or to keep track of how things are proceeding. Only used for spectral correction.
#' @param denormalise_residuals TRUE/FALSE to denormalise residuals using the max fluorescence value of the supplied eemlist. Default to FALSE.
#' @param exclude_negative_residuals TRUE/FALSE to set negative values to 0 during component overlap correction. Setting this to TRUE may inflate similarity scores in component:sample comparisons where the model has over-fitted the component.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select
#' @importFrom dplyr mutate_at
#' @importFrom dplyr vars
#' @importFrom utils data
#'
#' @export
#'
per_eem_ssc <- function(pfmodel,
                        eemlist,
                        comp,
                        tcc = FALSE,
                        modified_metrics = TRUE,
                        spectral_correct = "all",
                        interp_1nm = "all",
                        smooth_sg = "all",
                        complete_peak = "ex",
                        SSC_trim_method = 'mod_alpha',
                        verbose = FALSE,
                        denormalise_residuals = FALSE,
                        exclude_negative_residuals = FALSE){
  # get PARAFAC spectra
  pf_peak_spectra <- extrpf_peak_spectra_int(pfmodel, component = comp)
  # get PARAFAC Emission 'B' mode - emission
  mat_pf_em_main <- as.matrix(pfmodel$B[,comp])
  # get PARAFAC Excitation 'C' mode - excitation
  mat_pf_ex_main <- as.matrix(pfmodel$C[,comp])
  # Interpolate PARAFAC spectra, if specified.
  if(!isTRUE(interp_1nm)){
    if(isTRUE(verbose)){
      message("Interpolating PARAFAC spectra, type: ", interp_1nm)
    }
    if(interp_1nm == "ex"){
      mat_pf_ex_main <- interpolate_matrix_1nm(mat = mat_pf_ex_main, mat_out = TRUE)
    } else if(interp_1nm == "em"){
      mat_pf_em_main <- interpolate_matrix_1nm(mat = mat_pf_em_main, mat_out = TRUE)
    } else if(interp_1nm == "all"){
      mat_pf_ex_main <- interpolate_matrix_1nm(mat = mat_pf_ex_main, mat_out = TRUE)
      mat_pf_em_main <- interpolate_matrix_1nm(mat = mat_pf_em_main, mat_out = TRUE)
    } else {
      stop("Please specify interp_1nm as NULL, 'ex', 'em' or 'all'")
    }
  }
  # Confirm only one component was extracted
  if(length(unique(pf_peak_spectra$max_ex)) > 1 | length(unique(pf_peak_spectra$max_em)) > 1 ){
    message("Oops! More than one spectra extracted.")
  }
  # Get target wavelength pair
  target_em <- pf_peak_spectra$max_em[1]
  target_ex <- pf_peak_spectra$max_ex[1]
  # If correcting spectra, extract the residuals and ex/em grobs.
  if(!is.null(spectral_correct)){
    if(isTRUE(verbose)){
      message("Extracting residuals data for component overlap correction...")
    }
    residuals <- extrpf_residuals_int(pfmodel = pfmodel, eem_list = eemlist, denormalise = denormalise_residuals, verbose = verbose)
    # Get the grobs for ex and em
    grob_ex <- residuals[which(residuals$em == target_em),]
    grob_em <-  residuals[which(residuals$ex == target_ex),]
    if(isTRUE(verbose)){
      message("Residuals extracted.")
    }
  }
  # Fill SSC table
  types <- c("tcc","ssc","alpha","beta")
  cats <- c("excitation","emission")
  types_pst <- paste(rep(cats, each = length(types)), types, sep = "_")
  SSC_table <- data.frame(matrix(NA, nrow = length(eemlist), ncol = length(types_pst)))
  colnames(SSC_table) <- types_pst
  if(isTRUE(modified_metrics)){
    SSC_table$mSSC <- NA
    SSC_table$mTCC <- NA
  }
  SSC_table$sample <- unlist(lapply(eemlist,"[[","sample"))
  # Calculate metrics for each EEM.
  for(e in seq_along(eemlist)){
    # For this EEM, pull out the emission and excitation slices.
    target_eem <- eemlist[[e]]
    name <- target_eem$sample
    mat_pf_em <- mat_pf_em_main # reassign the pf spectra to ensure no modification between loops
    mat_pf_ex <- mat_pf_ex_main
    eem_slice <- slice_eem_int(eem = target_eem, ex = target_ex, em = target_em)
    mat_em_it <- eem_slice %>%
      dplyr::filter(name == "emission") %>%
      tibble::column_to_rownames('wavelength') %>%
      select(-2) %>%
      mutate_at(vars(intensity), as.numeric) %>%
      data.matrix()
    mat_ex_it <- eem_slice %>%
      dplyr::filter(name == "excitation") %>%
      tibble::column_to_rownames('wavelength') %>%
      select(-2) %>%
      mutate_at(vars(intensity), as.numeric) %>%
      data.matrix()
    # EEM Correction step 1: Removal of fluorescence contribution from non-target components; component spectral overlap correction.
    if(!is.null(spectral_correct)){
      if(spectral_correct == "ex"){
        mat_ex_it <- comp_correct_spectra(grob = grob_ex, sample_char = name, comp = comp, type = "ex", neg_to_0 = exclude_negative_residuals)
      } else if(spectral_correct == "em"){
        mat_em_it <- comp_correct_spectra(grob = grob_em, sample_char = name, comp = comp, type = "em", neg_to_0 = exclude_negative_residuals)
      } else if(spectral_correct == "all"){
        mat_ex_it <- comp_correct_spectra(grob = grob_ex, sample_char = name, comp = comp, type = "ex", neg_to_0 = exclude_negative_residuals)
        mat_em_it <- comp_correct_spectra(grob = grob_em, sample_char = name, comp = comp, type = "em", neg_to_0 = exclude_negative_residuals)
      } else {
        stop("Please specify spectral_correct as NULL, 'ex', 'em' or 'all'")
      }
    }
    # Correction step 2: interpolation to 1nm bandwidth
    if(!isTRUE(interp_1nm)){
      if(interp_1nm == "ex"){
        mat_ex_it <- interpolate_matrix_1nm(mat = mat_ex_it, mat_out = TRUE)
      } else if(interp_1nm == "em"){
        mat_em_it <- interpolate_matrix_1nm(mat = mat_em_it, mat_out = TRUE)
      } else if(interp_1nm == "all"){
        mat_ex_it <- interpolate_matrix_1nm(mat = mat_ex_it, mat_out = TRUE)
        mat_em_it <- interpolate_matrix_1nm(mat = mat_em_it, mat_out = TRUE)
      } else {
        stop("Please specify interp_1nm as NULL, 'ex', 'em' or 'all'")
      }
    }
    # Correction step 3: smoothing with Savitzky-Golay filter
    if(!isTRUE(smooth_sg)){
      data("sg_terms_deftab") # default smoothing terms. Eventually tweakable (I hope)
      if(smooth_sg == "ex"){
        mat_ex_it <- sg_smooth(mat = mat_ex_it,  n = sg_terms_deftab['ex','n'], p = sg_terms_deftab['ex','p'], m = sg_terms_deftab['ex','m'], ts = sg_terms_deftab['ex','ts'])
      } else if(smooth_sg == "em"){
        mat_em_it <- sg_smooth(mat = mat_em_it,  n = sg_terms_deftab['em','n'], p = sg_terms_deftab['em','p'], m = sg_terms_deftab['em','m'], ts = sg_terms_deftab['em','ts'])
      } else if(smooth_sg == "all"){
        mat_ex_it <- sg_smooth(mat = mat_ex_it,  n = sg_terms_deftab['ex','n'], p = sg_terms_deftab['ex','p'], m = sg_terms_deftab['ex','m'], ts = sg_terms_deftab['ex','ts'])
        mat_em_it <- sg_smooth(mat = mat_em_it,  n = sg_terms_deftab['em','n'], p = sg_terms_deftab['em','p'], m = sg_terms_deftab['em','m'], ts = sg_terms_deftab['em','ts'])
      } else {
        stop("Please specify smooth_sg as NULL, 'ex', 'em' or 'all'")
      }
    }
    # Secondary peak handling. Gradient peak detection to remove incomplete peaks.
    if(!is.null(complete_peak)){
      # Note: if complete_peak is specified (i.e. not NULL), then the peak handling occurs in one of two ways, based upon the SSC_trim method var.
      if(SSC_trim_method == 'mod_alpha'){
        # The default (and recommended) method is to adjust the peak position calculation, rather than discarding spectra associated with an incomplete peak.
        if(complete_peak == "ex"){
          trough_val_ex <- extract_complete_peak(mat1 = mat_pf_ex, mat2 = mat_ex_it, return_trough = TRUE)
        } else if(complete_peak == "em"){
          trough_val_em <- extract_complete_peak(mat1 = mat_pf_em, mat2 = mat_em_it, return_trough = TRUE)
        } else if(complete_peak == "all"){
          trough_val_ex <- extract_complete_peak(mat1 = mat_pf_ex, mat2 = mat_ex_it, return_trough = TRUE)
          trough_val_em <- extract_complete_peak(mat1 = mat_pf_em, mat2 = mat_em_it, return_trough = TRUE)
        } else {
          stop("Please specify complete_peak as NULL, 'ex', 'em' or 'all'")
        }
      } else if(SSC_trim_method == 'trim_spectra'){
        # Alternatively, the 'original' method can be used, where incomplete peaks are simply removed altogether.
        if(complete_peak == "ex"){
          ex_mats <- extract_complete_peak(mat1 = mat_pf_ex, mat2 = mat_ex_it)
          mat_pf_ex <- ex_mats$mat1
          mat_ex_it <- ex_mats$mat2
        } else if(complete_peak == "em"){
          em_mats <- extract_complete_peak(mat1 = mat_pf_em, mat2 = mat_em_it)
          mat_pf_em <- em_mats$mat1
          mat_em_it <- em_mats$mat2
        } else if(complete_peak == "all"){
          em_mats <- extract_complete_peak(mat1 = mat_pf_em, mat2 = mat_em_it)
          mat_pf_em <- em_mats$mat1
          mat_em_it <- em_mats$mat2
          ex_mats <- extract_complete_peak(mat1 = mat_pf_ex, mat2 = mat_ex_it)
          mat_pf_ex <- ex_mats$mat1
          mat_ex_it <- ex_mats$mat2
        } else {
          stop("Please specify complete_peak as NULL, 'ex', 'em' or 'all'")
        }
      } else {
        stop("Please specify SSC_trim_method as 'mod_alpha' or 'trim_spectra'")
      }
    }
    # Populating SSC table.
    # SSC table handling
    if(SSC_trim_method == 'trim_peak' || is.null(SSC_trim_method)){
      ## Use mats (trimmed or un-trimmed) to find the SSC value.
      ssc_more_em <-  ssc_more_int(mat_pf_em, mat_em_it, tcc = FALSE)
      tcc_em <- staRdom::ssc(mat_pf_em, mat_em_it, tcc = TRUE)
      ssc_more_ex <-  ssc_more_int(mat_pf_ex, mat_ex_it, tcc = FALSE)
      tcc_ex <- staRdom::ssc(mat_pf_ex, mat_ex_it, tcc = TRUE)
      SSC_table[e,c("excitation_tcc", "excitation_ssc", "excitation_alpha", "excitation_beta")] <- c(tcc_ex,ssc_more_ex[1],ssc_more_ex[2],ssc_more_ex[3])
      SSC_table[e,c("emission_tcc", "emission_ssc", "emission_alpha", "emission_beta")] <- c(tcc_em,ssc_more_em[1],ssc_more_em[2],ssc_more_em[3])
      if(isTRUE(modified_metrics)){
        SSC_table[e,'mSSC'] <- sqrt(as.numeric(ssc_more_ex['ssc',1]) * as.numeric(ssc_more_em['ssc',1]))
        SSC_table[e,'mTCC'] <- sqrt(as.numeric(tcc_ex) * as.numeric(tcc_em))
      }
    } else if(SSC_trim_method == 'mod_alpha'){
      # Modified alpha term calculation using the supplied trough value.
      # Depending on which spectra may be trimmed.
      if(complete_peak == "ex"){
        ssc_more_ex <- ssc_more_int(mat_pf_ex, mat_ex_it, tcc = FALSE, alpha_start = trough_val_ex)
        ssc_more_em <-  ssc_more_int(mat_pf_em, mat_em_it, tcc = FALSE)
      } else if(complete_peak == "em"){
        ssc_more_ex <- ssc_more_int(mat_pf_ex, mat_ex_it, tcc = FALSE)
        ssc_more_em <- ssc_more_int(mat_pf_em, mat_em_it, tcc = FALSE, alpha_start = trough_val_em)
      } else if(complete_peak == "all"){
        ssc_more_ex <- ssc_more_int(mat_pf_ex, mat_ex_it, tcc = FALSE, alpha_start = trough_val_ex)
        ssc_more_em <- ssc_more_int(mat_pf_em, mat_em_it, tcc = FALSE, alpha_start = trough_val_em)
      } else {
        stop("Please specify complete_peak as NULL, 'ex', 'em' or 'all'")
      }
      tcc_ex <- staRdom::ssc(mat_pf_ex, mat_ex_it, tcc = TRUE)
      tcc_em <- staRdom::ssc(mat_pf_em, mat_em_it, tcc = TRUE)
      SSC_table[e,c("excitation_tcc", "excitation_ssc", "excitation_alpha", "excitation_beta")] <- c(tcc_ex,ssc_more_ex[1],ssc_more_ex[2],ssc_more_ex[3])
      SSC_table[e,c("emission_tcc", "emission_ssc", "emission_alpha", "emission_beta")] <- c(tcc_em,ssc_more_em[1],ssc_more_em[2],ssc_more_em[3])
      if(isTRUE(modified_metrics)){
        # missing value handling to prevent the occurrence of positive value calculation were two negatives are supplied.
        if((ssc_more_ex['ssc',1] < 0) || (ssc_more_em['ssc',1] < 0) || (is.na(ssc_more_ex['ssc',1])) || (is.na(ssc_more_em['ssc',1]))){
          SSC_table[e,'mSSC'] <- NA
        } else {
          SSC_table[e,'mSSC'] <- sqrt(as.numeric(ssc_more_ex['ssc',1]) * as.numeric(ssc_more_em['ssc',1]))
        }
        if((as.numeric(tcc_ex) < 0) || (as.numeric(tcc_em) < 0) || (is.na(tcc_ex)) || (is.na(tcc_em))){
          SSC_table[e,'mTCC'] <- NA
        } else {
          SSC_table[e,'mTCC'] <- sqrt(as.numeric(tcc_ex) * as.numeric(tcc_em))
        }
      }
    }
    if(isTRUE(verbose)){
      message(name)
    }
    if(e %% 50 == 0){
      message("EEM ", e, " complete")
    }
  }
  message("Calculations complete")
  SSC_table
}
