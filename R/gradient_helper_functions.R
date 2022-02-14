# Functions used in the correction function extract_complete_peaks()
# These are all variously modified ports of code originally used in the RamanIntegrationRange function from
# the drEEM MATLAB package.
# Function originally presented in: Utilises code from Murphy, K. R. (2011). notes A Note on Determining the Extent of the Water Raman Peak in Fluorescence Spectroscopy. Applied Spectroscopy, 65(2), 233–236. https://doi.org/10.1366/10-06136.
# MATLAB toolbox: Murphy, K. R., Stedmon, C. A., Graeber, D. and Bro, R.: Fluorescence spectroscopy and multi-way techniques. PARAFAC, Anal. Methods, 5, 6557-6566, 2013.


#' Find the wavelength values of positive and negative gradients in a matrix.
#'
#' @description Use a discrete numerical gradient to find all gradient values in a given matrix.
#'
#' @param mat a data.matrix
#' @param tolerance a decimal representing percentage tolerance. 0.01 (1%) by default.
#'
#' @importFrom stats smooth
#'
#' @noRd
#'
get_gradient_vals <- function(mat, tolerance = 0.01){
  RAMmat <- mat
  RAMmat <- as.matrix(t(RAMmat))
  col_ok = apply(RAMmat,2,function(x)!any(is.na(x))) # identify columns containing NAs, assign logical T/F
  RAMmat <- as.matrix(RAMmat[,col_ok]) # remove NAs
  # Next: passing to MATLAB code.
  # Perform checks and read functions.
  if(ncol(RAMmat) < nrow(RAMmat)){ # check that the data is in rows. If not, transpose.
    RAMmat = t(RAMmat)
  } # Ensure R is in rows, not columns. Transposes if in columns.
  # Use ispos to confirm wavelengths are in the first row of RAMmat.
  if(ispos(Xmatrix = pracma::gradient(RAMmat[1,])) == 0){
    warning("check that wavelengths are in the first row of RAMmat, and that scans for different samples are in subsequent rows!")
  }
  # Next: iterate along n samples, determining gradients.
  NoSamples = as.numeric(length(RAMmat[,1])) - 1 # n samples. Leave for now.
  gradient_list <- vector("list", length = NoSamples)
  # for(z in seq_along(gradient_list)){
  # zn <- formatC(z, width = 2, format = "d", flag = "0") # formatted sample numbers (for <100 sample groups)
  # tolcheck <- tolerance.
  # tolerance; what's the minimum acceptable gradient threshold? 1% to start.
  tolcheck <- tolerance
  # R <- RAMmat[c(1,z+1),] # pull out sample data.
  R <- RAMmat
  x = R[1,]
  y = R[2,]
  g <- pracma::gradient(R[2,]) # determine gradient, return y output.
  sg = forecast::ma(g, 5, TRUE) #5 point moving average
  sg2 = smooth(g)
  sg <- as.matrix(sg)
  #Identify sequences of positive and negative gradients
  # wtf is sequence?
  sequence = 8
  p = matrix(0,nrow=length(sg),ncol=1)
  n = p
  sq = sequence - 1
  for (i in seq_along(1:length(sg))){
    sgsub_pos = sg[i]
    p[i] = ispos(sgsub_pos)
  } # logical 1/0 for positives
  for (i in seq_along(1:length(sg))){
    sgsub_neg = sg[i]
    n[i] = isneg(sgsub_neg)
  } # logical 1/0 for negatives
  sgp = sg*p
  sgn = sg*n
  # restrict gradient bounds to avoid rayleigh scatter in determining the maximum gradient
  min_wl_ind <- as.numeric(which(R[1,] == min(R[1,])))
  max_wl_ind <- as.numeric(which(R[1,] == max(R[1,])))
  sgmax = max(sg[x > min_wl_ind], na.rm = TRUE)
  # incorporate tolerance.
  xi = x[sgp - tolcheck*sgmax > 0]
  sgpi = sgp[sgp - tolcheck*sgmax > 0]
  xf = x[sgn + tolcheck*sgmax < 0]
  sgnf = sgn[sgn + tolcheck*sgmax < 0]
  # xi are the indices of positive gradients. xf of negative gradient locations.
  # format, add to list
  grad_list_it <- list(xi,xf)
  names(grad_list_it) <- c("pos","neg")
  grad_list_it
}

#' Determine the number of discrete gradient segments.
#'
#' @description given a set of gradients, return the number of uninterrupted positive and negative gradient sequences.
#'
#' @param grads gradient values, returned from get_gradient_vals()
#'
#' @noRd
#'
nsegments <- function(grads){
  segments <- lapply(1:length(grads), function(g){
    type <- grads[[g]]
    diffs <- diff(type)[which(!is.na(diff(type)))]
    nseg <- as.numeric(length(diffs[diffs > 1])+1)
  })
  names(segments) <- names(grads)
  segments
}

#' Use data.table to find the nearest value in a vector.
#'
#' @description Using binary search provided by the data table package, find the nearest value. Originally from https://github.com/MRPHarris/eemUtils.
#'
#' @param data a vector of numeric values to search within
#' @param value to drive binary search. i.e. output will be closest value to this input.
#'
#' @importFrom data.table data.table
#' @importFrom data.table setattr
#' @importFrom data.table setkey
#'
#' @noRd
#'
binary_search_nearest_int <- function(data, value){
  if(!is.numeric(data)){
    data <- as.numeric(data)
  }
  if(!is.numeric(value)){
    value <- is.numeric(value)
  }
  data_tab <- data.table(data, val = data)
  setattr(data_tab, "sorted","data")
  setkey(data_tab, data)
  nearest <- as.numeric(data_tab[J(value), roll = "nearest"][1,2])
  nearest
}

#' Determine if a matrix contains only positive values
#'
#' @description Positive value determination in a matrix
#'
#' @param Xmatrix the input matrix
#'
#' @noRd
#'
ispos <- function(Xmatrix){
  has.neg <- apply(as.data.frame(Xmatrix), 1, function(col) any(col <= 0))
  if(any(has.neg == TRUE, na.rm = TRUE)){
    TrueOrFalse = 0
  } else {
    TrueOrFalse = 1
  }
  TrueOrFalse
}

#' Determine if a matrix contains only negative values
#'
#' @description Negative value determination in a matrix
#'
#' @param Xmatrix the input matrix
#'
#' @noRd
#'
isneg <- function(Xmatrix){
  has.pos <- apply(as.data.frame(Xmatrix), 1, function(col) any(col > 0))
  if(any(has.pos == TRUE, na.rm = TRUE)){
    TrueOrFalse = 0
  } else {
    TrueOrFalse = 1
  }
  TrueOrFalse
}


#
