#' Coerce to and check for `dsgn` objects
#'
#' @description
#' `dsgn` objects are recognized as design matrices with shifted, spatially
#' distributed time series by other functions in the `stcd` package. Coerce objects
#' to the `dsgn` class using `as.dsgn()` and check if something is a `dsgn` object
#' with `is.dsgn()`.
#'
#' NOTE: To create new `dsgn` objects from data, see [shift_to_design()]
#'
#'
#' @param x A data.frame or coercible object with shifted time series columns labeled `VAR_tminus#`
#'
#' @returns
#' * `as.dsgn()` returns a `dsgn` object
#'
#' * `is.dsgn()` returns TRUE if its argument is a design matrix object (that is,
#' has `"dsgn"` amongst its classes) and FALSE otherwise
#'
#' @seealso
#' * [as.discrete_dsgn()] and [as.coord_dsgn()] for the coercion to the subclasses of
#' discrete-indexed spatial design matrices and coordinate-indexed spatial design matrices,
#' respectively.
#'
#' * [shift_to_design()] for creation of new design matrices from data
#' @name dsgn
#' @rdname dsgn
#' @export
#'
#' @examples
#' design <- matrix(nrow = 2, ncol = 2)
#' is.dsgn(design)
#'
#' # as.dsgn will throw an error if columns do not represent time series
#' # (i.e. if no column name contains `tminus`)
#' try(design <- as.dsgn(design))
#'
#' colnames(design) <- c("var", "var_tminus1")
#' design <- as.dsgn(design)
#' is.dsgn(design)
as.dsgn <- function(x){
  if(!is.data.frame(x)){
    x <- as.data.frame(x)
  }
  if(any(grepl("tminus", colnames(x)))){
    class(x) <- c("dsgn", class(x))}
  else{stop("Provided object does not contain shifted time series columns")}
  x
}
#' @rdname dsgn
#' @export
is.dsgn <- function(x){
  "dsgn" %in% class(x)
}

#barebones helper for other functions
dsgn <- function(x) {
  class(x) <- c("dsgn", class(x))
  x
}
