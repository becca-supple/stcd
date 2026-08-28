#' PAG class -- WIP
#'
#' @param x A data frame representing nodes and edges in a time series partial ancestral graph
#'
#' @returns A `PAG` object
#' @export
#'
#' @name PAG
#' @rdname PAG
#' @examples
as.PAG <- function(x){
  class(x) <- c("PAG", "t_graph", class(x))
  x
}
#' @rdname PAG
#' @export
is.PAG <- function(x){
  "PAG" %in% class(x)
}

#barebones helper for other functions
PAG <- function(x){
  class(x) <- c("PAG", "t_graph", class(x))
  x
}
