#' @importFrom Rcpp evalCpp
#' @useDynLib adjVUS, .registration = TRUE

#' @export
vus <- function(x, y, z){
  if(any(is.na(x)) | any(is.na(y)) | any(is.na(z))) return(NA)
  else return(vusC_full(x, y, z))
}
