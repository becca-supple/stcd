#' CPDAG class -- WIP
#'
#' @param x A data frame representing nodes and edges in a time series complete
#' partially directed acyclic graph
#'
#' @returns A `CPDAG` object
#' @export
#'
#' @name CPDAG
#' @rdname CPDAG
#' @examples
as.CPDAG <- function(x){
  class(x) <- c("CPDAG", "t_graph", class(x))
  x
}
#' @rdname CPDAG
#' @export
is.CPDAG <- function(x){
  "CPDAG" %in% class(x)
}

#barebones helper for other functions
CPDAG <- function(x){
  class(x) <- c("CPDAG", "t_graph", class(x))
  x
}
