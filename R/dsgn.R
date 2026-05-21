#' Functions associated with 'dsgn' objects
#'
#' @param x A data.frame or coercible object
#'
#' @returns A 'dsgn' object or TRUE/FALSE
#' @export
#'
#' @examples
#' design <- matrix(nrow = 2, ncol = 2)
#' design <- dsgn(design)
#' class(design)
dsgn <- function(x) {
  class(x) <- c("dsgn", class(x))
  x
}

