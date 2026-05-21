#' Functions associated with 'dsgn' objects
#'
#' @param x A data.frame or coercible object with shifted time series columns labeled 'VAR_tminus#'
#'
#' @returns A 'dsgn' object or TRUE/FALSE
#' @export
#'
#' @examples
#' design <- matrix(nrow = 2, ncol = 2)
#' is.dsgn(design)
#'
#' # as.dsgn will throw an error if columns do not represent time series (i.e. if no column name contains 'tminus')
#' try(design <- as.dsgn(design))
#'
#' colnames(design) <- c("var", "var_tminus1")
#' design <- as.dsgn(design)
#' is.dsgn(design)
dsgn <- function(x) {
  class(x) <- c("dsgn", class(x))
  x
}
#' @export
as.dsgn <- function(x){
  if(!is.data.frame(x)){
    x <- as.data.frame(x)
  }
  if(grepl("tminus", colnames(x))){
    class(x) <- c("dsgn", class(x))}
  else{stop("Provided object does not contain shifted time series columns")}
  x
}
#' @export
is.dsgn <- function(x){
  "dsgn" %in% class(x)
}
