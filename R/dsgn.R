#' Create and check for `dsgn` objects
#'
#' @description
#' `dsgn` objects are recognized as design matrices with shifted, spatially distributed time series by other functions in the `stcd` package. Create `dsgn` objects using `as.dsgn()` and check if something is a `dsgn` object with `is.dsgn()`.
#'
#'
#' @param x A data.frame or coercible object with shifted time series columns labeled `VAR_tminus#`
#'
#' @returns
#' * `as.dsgn()` returns a `dsgn` object
#'
#' * `is.dsgn()` returns TRUE if its argument is a design matrix object (that is, has `"dsgn"` amongst its classes) and FALSE otherwise
#'
#' @export
#'
#' @examples
#' design <- matrix(nrow = 2, ncol = 2)
#' is.dsgn(design)
#'
#' # as.dsgn will throw an error if columns do not represent time series (i.e. if no column name contains `tminus`)
#' try(design <- as.dsgn(design))
#'
#' colnames(design) <- c("var", "var_tminus1")
#' design <- as.dsgn(design)
#' is.dsgn(design)
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
#' @export
dsgn <- function(x) {
  class(x) <- c("dsgn", class(x))
  x
}
