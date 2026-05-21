#' Coerce to and check for `discrete_dsgn` objects
#'
#' @description
#' `discrete_dsgn` objects are recognized as design matrices with shifted time series in discrete space by other functions in the `stcd` package. Coerce objects to the `discrete_dsgn` class using `as.dsgn()` and check if something is a `dsgn` object with `is.dsgn()`.
#'
#' NOTE: To create a new `discrete_dsgn` object from data with a site ID column, use [`shift_to_design()`] with a non-NULL value for `site_col`.
#'
#'
#' @param x A data.frame or coercible object with shifted time series columns labeled `VAR_tminus#`
#'
#' @returns
#' * `as.discrete_dsgn()` returns a `discrete_dsgn` object
#'
#' * `is.discrete_dsgn()` returns TRUE if its argument is a design matrix object
#' (that is, has `"discrete_dsgn"` amongst its classes) and FALSE otherwise
#'
#' @seealso
#' * [as.dsgn()] for the umbrella class
#'
#' * [as.coord_dsgn()] for coercion to coordinate-indexed spatial design matrices
#'
#' * [shift_to_design()] for creation of new design matrices from data
#'
#' @name discrete_dsgn
#' @rdname discrete_dsgn
#' @export
#'
#' @examples
#' design <- matrix(nrow = 2, ncol = 2)
#' is.discrete_dsgn(design)
#'
#' # as.discrete_dsgn will throw an error if columns do not represent time series
#' # (i.e. if no column name contains `tminus`)
#' try(design <- as.discrete_dsgn(design))
#'
#' # it will not, however, throw any errors if you accidentally coerce something
#' # with a non-discrete spatial structure.
#' # to create a new discrete_dsgn object from data, see shift_to_design()
#' colnames(design) <- c("var", "var_tminus1")
#' design <- as.discrete_dsgn(design)
#' is.discrete_dsgn(design)
#'
#' # discrete_dsgn objects are also dsgn objects
#' is.dsgn(design)
as.discrete_dsgn <- function(x){
  if(!is.data.frame(x)){
    x <- as.data.frame(x)
  }
  if(any(grepl("tminus", colnames(x)))){
    class(x) <- c("discrete_dsgn", "dsgn", class(x))}
  else{stop("Provided object does not contain shifted time series columns")}
  x
}
#' @rdname discrete_dsgn
#' @export
is.discrete_dsgn <- function(x){
  "discrete_dsgn" %in% class(x)
}

#barebones helper for other functions
discrete_dsgn <- function(x) {
  class(x) <- c("discrete_dsgn", "dsgn", class(x))
  x
}
