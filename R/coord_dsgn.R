#' Coerce to and check for `coord_dsgn` objects
#'
#' @description
#' `coord_dsgn` objects are recognized as design matrices with shifted time series in coordinate/continuous space by other functions in the `stcd` package. Coerce objects to the `coord_dsgn` class using `as.dsgn()` and check if something is a `dsgn` object with `is.dsgn()`.
#'
#' NOTE: To create a new `coord_dsgn` object from data with coordinates defined by two columns, use [`shift_to_design()`] with a vector naming the two axes of coordinates for `coord`.
#'
#'
#' @param x A data.frame or coercible object with shifted time series columns labeled `VAR_tminus#`
#'
#' @returns
#' * `as.coord_dsgn()` returns a `coord_dsgn` object
#'
#' * `is.coord_dsgn()` returns TRUE if its argument is a design matrix object (that is, has `"coord_dsgn"` amongst its classes) and FALSE otherwise
#'
#' @seealso
#' * [as.dsgn()] for the umbrella class
#'
#' * [as.discrete_dsgn()] for coercion to discrete-indexed spatial design matrices
#'
#' * [shift_to_design()] for creation of new design matrices from data
#'
#' @name coord_dsgn
#' @rdname coord_dsgn
#' @export
#'
#' @examples
#' design <- matrix(nrow = 2, ncol = 2)
#' is.coord_dsgn(design)
#'
#' # as.coord_dsgn will throw an error if columns do not represent time series
#' # (i.e. if no column name contains `tminus`)
#' try(design <- as.coord_dsgn(design))
#'
#' # it will not, however, throw any errors if you accidentally coerce something
#' # with a non-coordinate spatial structure.
#' # to create a new coord_dsgn object from data, see shift_to_design()
#' colnames(design) <- c("var", "var_tminus1")
#' design <- as.coord_dsgn(design)
#' is.coord_dsgn(design)
#'
#' # coord_dsgn objects are also dsgn objects
#' is.dsgn(design)
as.coord_dsgn <- function(x){
  if(!is.data.frame(x)){
    x <- as.data.frame(x)
  }
  if(any(grepl("tminus", colnames(x)))){
    class(x) <- c("coord_dsgn", "dsgn", class(x))}
  else{stop("Provided object does not contain shifted time series columns")}
  x
}
#' @rdname coord_dsgn
#' @export
is.coord_dsgn <- function(x){
  "coord_dsgn" %in% class(x)
}

#barebones helper for other functions
coord_dsgn <- function(x) {
  class(x) <- c("coord_dsgn", "dsgn", class(x))
  x
}
