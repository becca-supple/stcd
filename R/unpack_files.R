#' Unpack simulations from bulk files into the environment
#'
#' @param type A vector of characters indicating what to unpack - any
#' combination of c("seed", "design", "ind", "or") or "all"
#' @param path String indicating path to folder where bulk files are kept
#' @param n_sims Number of simulated scenarios
#' @param n_reps Number of repetitions per scenario
#' @param discrete TRUE for discrete space, FALSE for continuous
#' @param sufficiency TRUE for sufficient simulations, FALSE for nonsufficient
#'
#' @returns Adds to the environment lists of all the indicated files of the
#' indicated types
#' @export
unpack_files <- function(type = "all",
                         path,
                         n_sims,
                         n_reps,
                         discrete = FALSE,
                         sufficiency = TRUE){

  if(type == "all"){
    type <- c("seed", "design", "ind", "or")
  }

  if("seed" %in% type){
    seedlist <- list()
  }

  if("design" %in% type){
    designlist <- list()
  }

  if("ind" %in% type){
    indlist <- list()
  }

  if("or" %in% type){
    orlist <- list()
  }

  for(i in seq(1, n_sims)){

    for(r in seq(1, n_reps)){

      name <- paste0("data_", i, "_", r)

      if(file.exists(paste0(path, name, ".RData")) == FALSE){
        next
      }

      load(paste0(path, name, ".RData"))
      seed <- get(name)[[1]]
      if(discrete){
        design <- discrete_dsgn(get(name)[[2]])
      }else{
        design <- coord_dsgn(get(name)[[2]])
      }
      ind <- get(name)[[3]]
      if(sufficiency){
        or <- CPDAG(get(name)[[4]])
      }else{
        or <- PAG(get(name)[[4]])
      }

      if("seed" %in% type){

        seedlist[[paste0(i, "_", r)]] <- seed

      }

      if("design" %in% type){

        designlist[[paste0(i, "_", r)]] <- design

      }

      if("ind" %in% type){

        indlist[[paste0(i, "_", r)]] <- ind

      }

      if("or" %in% type){

        orlist[[paste0(i, "_", r)]] <- or

      }

      rm(list = name)

    }
  }

  to_return <- list()

  for(t in type){

    to_return[[length(to_return) + 1]] <- get(paste0(t, "list"))

  }

  names(to_return) <- paste0(type, "list")

  list2env(to_return, .GlobalEnv)

}
