dsgn <- function(x) {
  class(x) <- c("dsgn", class(x))
  x
}

as.dsgn <- function(x) {
  class(x) <- c("dsgn", class(x))
  x
}

is.dsgn <- function(x){
  "dsgn" %in% class(x)
}
