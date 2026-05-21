as.dsgn <- function(x) {
  class(x) <- c("dsgn", class(x))
  x
}
