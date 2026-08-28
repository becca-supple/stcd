#' t_graph class -- WIP
#'
#' @param x A data frame representing nodes and edges in a time series graph
#'
#' @returns A `t_graph` object
#' @export
#'
#' @name t_graph
#' @rdname t_graph
#' @examples
as.t_graph <- function(x){
  class(x) <- c("t_graph", class(x))
  x
}
#' @rdname t_graph
#' @export
is.t_graph <- function(x){
  "t_graph" %in% class(x)
}

#barebones helper for other functions
t_graph <- function(x){
  class(x) <- c("t_graph", class(x))
  x
}
