#' Find conditional independencies of spatially-structured time series
#'
#' @param data Either a data.frame of time series and spatial information, or a
#' `dsgn` object as created by [shift_to_design]
#' @param max_lag Maximum time lag at which to test relationships
#' @param sites String naming column of site names, if in discrete space.
#' Otherwise, `NULL`.
#' @param coords Vector of strings naming columns with coordinate information,
#' if in coordinate space. Otherwise, `NULL`.
#' @param time String naming column with time information, if data.frame provided
#' @param family Ordered vector of strings indicating the distributions/family
#' of each time series. If one string is provided, it is assumed that is the
#' family of all variables.
#' @param modeltype One of `"separated"`, `"tensor"`, or `"both"`. See
#' [fit_custom_gam].
#' @param neighborhood List defining neighborhood structure (only in discrete
#' space). See [build_neighborhood].
#' @param alpha Significance level
#' @param collider_rule One of `"conservative"`, `"majority"`, or `"none"`. See
#' details.
#' @param independence_test Either `"GCM"` or `"WGCM"`. See [WGCM_fix].
#' @param k0 Number of quantiles for [WGCM_fix].
#' @param nsim Number of simulations from which to create null distribution for
#' conditional independence tests.
#' @param sufficiency `TRUE` if all relevant variables observed. `FALSE` by
#' default.
#' @param silent If `FALSE`, reports more verbosely on function status. `TRUE` by
#' default.
#' @param skeleton_names If `TRUE`, returns a skeleton with nodes labeled by
#' variable name rather than index. `FALSE` by default. Do not set this to `TRUE`
#' if you plan to run the output through [orient_edges].
#' @param skeleton_only If `TRUE`, interrupts the function before any collider
#' detection and returns the skeleton. `FALSE` by default.
#'
#' @description Function to find the conditional independencies and identify
#' colliders in a graph of time series with spatial structure. Heavily inspired
#' by the logic of PCMCI+ (Runge, 2020).
#'
#' @section Collider rules:
#' To improve computational and statistical efficiency, PC-based algorithms
#' search for the minimum separating set for each pair of variables, rather than
#' for all separating sets. That is, the exclusion of a variable Xk from a
#' separating set Aij found in a PC-based algorithm does not in all cases
#' indicate that Xi ⊥ Xj | Xk, but just that a smaller set was found first.
#'
#' Users can indicate how strictly they want to determine if Xk is a collider in
#' the structure Xi--Xk--Xj:
#'
#' * `collider_rule = "conservative"`: Xk is not in any separating sets. Default
#' for `conditional_independencies`.
#'
#' * `collider_rule = "majority"`: Xk is not in the majority of separating sets
#'
#' * `collider_rule = "none"`: Xk is not in the smallest/first found separating
#' set. Original rule and standard for PC-based algorithms.
#'
#'
#' @returns A list:
#' \describe{
#'    \item{`skeleton`}{data.frame describing graph structure. `Xj` and `Xi`
#'    indicate all pairs of nodes. `Xj_mark` and `Xi_mark` indicate the edge marks
#'    at each end of the edge connecting the pair, or are `NA` if no edge exists.}
#'    \item{`independencies`}{data.frame listing the independence tests run and
#'    their results. Test statistics associated with tests of `Xi` is independent
#'    of `Xj` given `S` are stored in the `Imin` column.}
#'    \item{`triples`}{data.frame with triples `Xi_tau`--`Xk_t`--`Xj_t`. If
#'    direction is yet to be resolved, `ambiguous` = 1. The `flag` column is
#'    used internally to handle potential conflicts.}
#'    \item{`names`}{vector of column names, if `skeleton_names = FALSE`}
#'    \item{`sep_sets`}{data.frame listing, for all pairs of nodes `Xj`, `Xi`,
#'    the set conditional on which they are independent. The `S_sep` column
#'    contains either the vector of indices comprising the separating set, `"null"`
#'    if unconditionally independent, or `NA` if no such set was found.}
#' }
#'
#' @export
#'
#' @examples
#'
#' head(example_coord_data)
#' #reducing dimensions for example speed
#' smaller_data <- example_coord_data[,-c(1:5)]
#' head(smaller_data)
#'
#'
#' \dontrun{
#' example_inds <- conditional_independencies(data = smaller_data,
#'                   max_lag = 2,
#'                   coords = c("x", "y"),
#'                   family = c("gaussian", "gaussian", "betar"),
#'                   modeltype = "separated",
#'                   alpha = 0.025,
#'                   collider_rule = "none", #for example speed, not recommended
#'                   independence_test = "GCM",
#'                   sufficiency = FALSE #since we removed the Agr time series
#'                   )
#'
#' str(example_inds)
#' }
#'
conditional_independencies <- function(data,
                                       max_lag = 2,
                                       sites = NULL,
                                       coords = NULL,
                                       time = "Year",
                                       family = "gaussian",
                                       modeltype = "separated",
                                       neighborhood = NULL,
                                       alpha = 0.05,
                                       collider_rule = "conservative",
                                       independence_test = "GCM",
                                       k0 = 7,
                                       nsim = 499,
                                       sufficiency = FALSE,
                                       silent = TRUE,
                                       skeleton_names = FALSE,
                                       skeleton_only = FALSE){

  ####
  # Force variables needed for later functions
  ####

  modeltype <- modeltype
  k0 <- k0
  nsim <- nsim
  max_lag <- max_lag
  if(independence_test %in% c("GCM", "WGCM") == FALSE){
    stop("\nError specifying independence test. Only GCM or WGCM independence tests supported.")
  }

  ####
  # Create variable set S & conditional independencies dataframe
  ####

  #Check for columns with NAs
  if(any(is.na(data))){
    stop("One or more NAs have been detected in the data.")
  }

  if("dsgn" %in% class(data) == FALSE){

    if(is.null(sites)){

      if(is.null(coords)){
        stop("Provided data must be a 'dsgn' object or be coercible to one via the sites or coord arguments.")
      }else{

        X <- seq(1, sum(colnames(data) %in% c(coords, time) == FALSE))
        n_X <- length(X)

        design <- shift_to_design(data,
                                  coord_cols = coords,
                                  time_col = time,
                                  max_lag = max_lag)

        space <- list(coords = coords)

      }

    }else{ #if column containing sites is provided

      #Set of variables indexed
      X <- seq(1, sum(colnames(data) %in% c(sites, time) == FALSE))
      n_X <- length(X)

      design <- shift_to_design(data,
                                site_col = sites,
                                time_col = time,
                                max_lag = max_lag)

      space <- list(sites = sites, nb = neighborhood)

    }
  }else{

    design <- data

    #Set of variables indexed
    if("discrete_dsgn" %in% class(design)){
      n_X <- (ncol(design) - 1) / (1 + (max_lag * 2))
      space <- list(sites = sites, nb = neighborhood)
    }else{
      if("coord_dsgn" %in% class(design)){
        n_X <- (ncol(design) - 2) / (1 + (max_lag * 2))
        space <- list(coords = coords)
      }
    }

    X <- seq(1, n_X)

  }

  #Each variable assigned a distribution
  if(length(family) == 1){
    distributions <- rep(family, (2*max_lag + 1)*n_X)
  }else{
    if(length(family) == n_X){
      distributions <- unlist(lapply(family, rep, (2*max_lag + 1)))
    }else{
      stop("\nError. Family list length does not match number of variables.")
    }
  }

  #Residual storing environment
  res_store <- new.env()

  #All present versions of variables
  X_t <- 1 + (1 + 2*max_lag)*(X - 1)

  #List of lagged sets for each Xj_t
  Bj_ts <- list()


  #List of valid past indices
  pasts <- numeric(max_lag * n_X)
  ind <- 1

  for(tau in seq(1, max_lag)){
    for(i in X){
      pasts[ind] <- 1 + tau + (1 + 2*max_lag)*(i - 1)
      ind <- ind + 1
    }
  }

  #Separating sets dataframe
  sep_sets <- data.frame(Xj = c(rep(X_t, each = length(pasts)),
                                combn(X_t, 2)[1,]),
                         Xi = c(rep(pasts, n_X),
                                combn(X_t, 2)[2,]),
                         S_sep = NA)

  ###
  # Lagged Adjacencies
  ###

  for(j in seq(1:length(X_t))){

    #Define column in design relevant to Xj_t and its distribution

    Xj_t <- X_t[j]
    Xj_t_dist <- distributions[Xj_t]

    #Initialize lagged conditioning set
    Bj_t <- numeric(max_lag * n_X)
    ind <- 1

    for(tau in seq(1, max_lag)){
      for(i in X){
        Bj_t[ind] <- 1 + tau + (1 + 2*max_lag)*(i - 1)
        ind <- ind + 1
      }
    }

    #Dataframe to store independence test results
    store_ind <- data.frame(matrix(nrow = length(Bj_t),
                                   ncol = 4))

    colnames(store_ind) <- c("Xi_past", "Xj_t", "Imin", "remove")

    #Define test pairs, initialize test stat at Inf
    store_ind$Xj_t <- Xj_t
    store_ind$Xi_past <- Bj_t
    store_ind$Imin <- Inf
    store_ind$remove <- 0

    ####
    # Find lagged set
    ####

    #Initialize conditioning set size at 0
    p <- 0

    #While there are variables left to test
    while(length(Bj_t) - 1 >= p){

      #For each variable in the set
      for(Xi_past in Bj_t){

        #are there non-Xi_past vars in Bj_t?
        if(length(setdiff(Bj_t, Xi_past)) >= p){
          Bj_t_i <- setdiff(Bj_t, Xi_past)
        }else{
          next
        }

        #relevant index of independence storing dataframe
        i <- which(store_ind$Xi_past == Xi_past)

        #relevant distribution
        Xi_past_dist <- distributions[Xi_past]

        #Find the next largest conditioning set
        if(p == 0){
          #if empty set set model_on = "null"
          S <- "null"
        }else{
          #or define as first p values of Bj_t
          S <- Bj_t_i[1:p]
        }

        #Test if Xi_past _||_ Xj_t | S

        if(Xi_past %in% S || Xj_t %in% S){
          stop("Self-model in lagged testing")
        }

        if(!silent){cat("\nFitting ", colnames(design)[Xi_past], " ~ ", S)}

        #Fit and extract residuals for a null model of Xi_past
        Xi_past_S <- find_or_fit(design = design,
                                 response = Xi_past,
                                 model_on = S,
                                 distribution = Xi_past_dist,
                                 Env = res_store,
                                 modeltype = modeltype,
                                 spatial_str = space)

        if(!silent){cat("\nFitting ", colnames(design)[Xj_t], " ~ ", S)}

        #Fit and extract residuals for Xj_t | S
        Xj_t_S <- find_or_fit(design = design,
                              response = Xj_t,
                              model_on = S,
                              distribution = Xj_t_dist,
                              Env = res_store,
                              modeltype = modeltype,
                              spatial_str = space)

        if(independence_test == "GCM"){
          test_result <-  gcm.test(resid.XonZ = Xi_past_S,
                                   resid.YonZ =  Xj_t_S)
        }else{
          if(is.character(S)){
            zed <- as.numeric(design[[sites]])
          }else{
            zed <- design[,S]
          }

          test_result <- WGCM_fix(resid.XonZ = Xi_past_S,
                                  resid.YonZ =  Xj_t_S,
                                  Z = zed,
                                  k0 = k0,
                                  nsim = nsim)
        }


        #Replace test statistic if lower than previously calculated
        if(test_result$test.statistic < store_ind$Imin[i]){
          store_ind$Imin[i]  <- test_result$test.statistic
        }

        #Store separating set and mark for removal if p-value > alpha
        if(test_result$p.value > alpha){

          sep_sets <- sep_sets |>
            mutate(S_sep = case_when(Xj == Xj_t & Xi == Xi_past ~ paste(S, collapse = ";"),
                                     TRUE ~ S_sep))

          store_ind$remove[i]  <- 1
        }
      }

      #Remove those variables that were marked from Bj_t
      Bj_t <- setdiff(Bj_t, store_ind$Xi_past[store_ind$remove == 1])

      #store_ind <- store_ind[store_ind$Xi_past %in% Bj_t, ]

      #Sort the remaining set by Imin from largest to smallest
      Bj_t <- store_ind  |>
        filter(Xi_past %in% Bj_t) |>
        arrange(desc(Imin)) |>
        pull(Xi_past)

      #Increment the size of the conditioning set
      p <- p + 1

    }

    #Store the finalized lagged set
    Bj_ts[[j]] <- Bj_t

  }

  ###
  # Contemporaneous adjacencies
  ###

  #Initialize the graph
  G <- data.frame(Xj = c(rep(X_t, each = length(pasts)),
                         combn(X_t, 2)[1,]),
                  Xi = c(rep(pasts, n_X),
                         combn(X_t, 2)[2,]),
                  Xj_mark = NA,
                  Xi_mark = NA)

  #Add links for all contemporaneous variables
  G <- G |>
    mutate(Xj_mark = ifelse(Xi %in% X_t, "o", NA),
           Xi_mark = ifelse(Xi %in% X_t, "o", NA))

  #Add links from lagged adjacency step
  for(j in X){

    #if there were any
    if(length(Bj_ts[[j]]) > 0){

      #for each adjacency
      for(Xi in Bj_ts[[j]]){

        if(sufficiency == TRUE){
          #add a fully directed edge from it to the present var in G
          G$Xj_mark[G$Xj == X_t[j] & G$Xi == Xi] <- "<"
          G$Xi_mark[G$Xj == X_t[j] & G$Xi == Xi] <- "-"
        }else{
          #if not assuming sufficiency, add only an arrowhead into present var
          G$Xj_mark[G$Xj == X_t[j] & G$Xi == Xi] <- "<"
          G$Xi_mark[G$Xj == X_t[j] & G$Xi == Xi] <- "o"
        }

      }
    }
  }

  #Initialize all contemporaneous adjacency sets as all other present variables
  Aj_ts <- list()

  for(j in X){
    Aj_ts[[j]] <- X_t[X[-j]]
  }

  #Set up dataframe to store independence test results for all linked pairs
  store_ind <- G |>
    filter(is.na(Xj_mark) == FALSE) |>
    mutate(Imin = Inf,
           S = NA) |>
    select(-c(Xj_mark, Xi_mark))

  ###
  # Contemporaneous adjacency set testing
  ###

  p <- 0

  #While there is any linked pair for which the contemporaneous adjacency of Xj
  #is at least 1 greater than the separating set size
  while(any(sapply(Aj_ts, function(x) {
    if (length(x) == 0) return(0)  # empty set treated as zero count
    length(x[x %in% G[!is.na(G$Xj_mark), 1]])
  }) - 1 >= p)) {

    js <- which(sapply(Aj_ts, function(x) {
      if(length(x) == 0) return(0)
      length(x[x %in% G[!is.na(G$Xj_mark), 1]])
    }) - 1 >= p)

    pairs <- which(is.na(G$Xj_mark) == FALSE & G$Xj %in% X_t[js])

    #For each pair in that condition
    for(pair in pairs){

      Xj_t <- G$Xj[pair]
      j <- which(X_t == G$Xj[pair])
      Xi_tau <- G$Xi[pair]
      i <- ((Xi_tau - 1) %/% (1 + 2 * max_lag)) + 1

      #are there non-Xi_tau vars in Aj_t?
      Aj_t_i <- setdiff(Aj_ts[[j]], Xi_tau)

      #index in independence storing df
      ind <- which(store_ind$Xj == Xj_t & store_ind$Xi == Xi_tau)

      #Find all separating sets of size p in order
      if(p > 0){
        S_sets <- combn(as.list(Aj_t_i), p)
      }else{
        S_sets <- "null"
      }

      tested <- 0


      #While there are still untested sets
      while(tested < length(tested)){
        for(S in S_sets){

          #Set Z as the set containing S, lagged adj of Xj - Xi_tau, and lagged adj of Xi_tau
          Bj_t <- setdiff(Bj_ts[[j]], Xi_tau)

          #if tau = 0
          if(Xi_tau %in% X_t){
            Bi_t <- Bj_ts[[i]]

            #or if greater than 0
          }else{
            lag <- Xi_tau - X_t[i]
            Bi_t <- Bj_ts[[i]]  + lag

          }

          if(S == "null"){
            Z <- unique(c(Bj_t, Bi_t))
            if(length(Z) == 0){
              Z <- "null"
            }
          }else{
            Z <- unique(c(S, Bj_t, Bi_t))
          }

          #Test Xi_tau _||_ Xj_t | Z

          if(Xi_tau %in% Z  || Xj_t %in% Z){
            stop("Self-model in comtemp. testing")
          }

          if(!silent){cat("\nFitting ", colnames(design)[Xi_tau], " ~ ", Z)}

          #Fit and extract residuals for a null model of Xi_past
          Xi_tau_Z <- find_or_fit(design = design,
                                  response = Xi_tau,
                                  model_on = Z,
                                  distribution = distributions[Xi_tau],
                                  Env = res_store,
                                  modeltype = modeltype,
                                  spatial_str = space)

          if(!silent){cat("\nFitting ", colnames(design)[Xj_t], " ~ ", Z)}

          #Fit and extract residuals for Xj_t | Z
          Xj_t_Z <- find_or_fit(design = design,
                                response = Xj_t,
                                model_on = Z,
                                distribution = distributions[Xj_t],
                                Env = res_store,
                                modeltype = modeltype,
                                spatial_str = space)

          if(independence_test == "GCM"){
            test_result <-  gcm.test(resid.XonZ = Xi_tau_Z,
                                     resid.YonZ =  Xj_t_Z)
          }else{
            if(is.character(Z)){
              zed <- as.numeric(design[[sites]])
            }else{
              zed <- design[,Z]
            }
            test_result <- WGCM_fix(resid.XonZ = Xi_tau_Z,
                                    resid.YonZ =  Xj_t_Z,
                                    Z = zed,
                                    k0 = k0,
                                    nsim = nsim)
          }


          #Replace test statistic if lower than previously calculated
          if(test_result$test.statistic < store_ind$Imin[ind]){
            store_ind$Imin[ind]  <- test_result$test.statistic
          }

          #If p-value > alpha, delete link in G
          if(test_result$p.value > alpha){
            G[pair, c(3,4)] <- NA

            #and store S as the separating set
            sep_sets <- sep_sets |>
              mutate(S_sep = case_when(Xj == Xj_t & Xi == Xi_tau ~ paste(Z, collapse = ";"),
                                       Xj == Xi_tau & Xj == Xj_t ~ paste(Z, collapse = ";"),
                                       TRUE ~ S_sep))

            store_ind$S[ind] <- paste(Z, collapse = ";")
          }

          #Mark this set as tested
          tested <- tested + 1

          #Break if set was already found
          if(is.na(store_ind$S[ind]) == FALSE){
            tested <- length(tested)
            break
          }

        }

      }


    }

    #Increment p
    p <- p + 1

    #Reset contemporaneous sets based on the remaining links in G
    for(j in X){
      result <- store_ind |>
        filter(is.na(S)) |>
        filter(Xj == X_t[j] | Xi == X_t[j]) |>
        filter(Xi %in% X_t & Xj %in% X_t) |>
        arrange(desc(Imin)) |>
        mutate(other = ifelse(Xj == X_t[j], Xi, Xj)) |>
        pull(other)

      #Remove self (response) explicitly
      result <- setdiff(result, X_t[j])

      if(length(result) == 0){
        Aj_ts[[j]] <- numeric(0)  # assign empty numeric vector if no result
      }else{
        Aj_ts[[j]] <- result
      }
    }

  }

  ###
  # Format and return if you just want skeleton
  ###

  if(skeleton_only == TRUE){
    if(skeleton_names == TRUE){
      #replace indices with names for interpretability
      G <- G |>
        mutate(Xj = colnames(design)[Xj],
               Xi = colnames(design)[Xi])
    }

    list2return <- list(skeleton = G, independencies = store_ind, sep_sets = sep_sets)
    return(list2return)
  }

  G_old <- G #store skeleton for later

  ###
  # Otherwise, ID unshielded triples and test for colliders
  ###

  #Set up df to store indices of unsh triples
  unsh_triples <- data.frame(Xi_tau = NA, Xk_t = NA, Xj_t = NA, ambiguous = NA,
                             flag = 1)

  testing <- TRUE

  #While there are unshielded triples to test
  while(testing == TRUE){

    ###
    # ID unshielded triples
    ###

    contemp_adj <- G |>
      filter(is.na(Xj_mark) == FALSE) |>
      filter(Xi %in% X_t & Xj %in% X_t)

    if(nrow(contemp_adj) < 1){
      unsh_triples <- unsh_triples[-1,]
      testing <- FALSE
      break
    }

    for(row in seq(1, nrow(contemp_adj))){

      #trying with this orientation of Xj/Xi first, then will do opposite
      Xk_t <- contemp_adj$Xj[row]
      Xj_t <- contemp_adj$Xi[row]

      #Is there any Xi_tau that has an arrow or edge into Xk_t but not adj to Xj_t?

      Xi_tau <- G |>
        filter(Xj == Xk_t | Xi == Xk_t) |>
        filter(Xj != Xj_t & Xi != Xj_t) |>
        mutate(Xk_mark = case_when(Xj == Xk_t ~ Xj_mark,
                                   Xi == Xk_t ~ Xi_mark)) |>
        filter(Xk_mark == "<") |> #arrow or edge into Xk_t?
        mutate(other = ifelse(Xj == Xk_t, Xi, Xj)) |> #all other values
        pull(other)

      if(length(Xi_tau) == 0){
        next
      }

      Xi_tau <- G |> #no arrow or edge into Xj_t?
        filter(Xj == Xj_t | Xi == Xj_t) |>
        filter(is.na(Xj_mark)) |> #only those with no edge
        mutate(other = ifelse(Xj == Xj_t, Xi, Xj)) |> #all other values
        filter(other %in% Xi_tau) |>  #is Xi_tau in this list?
        pull(other)

      if(length(Xi_tau) == 0){
        next
      }

      #there are instance(s) of an unshielded triple, add to df
      for(Xi in Xi_tau){

        if(is.na(unsh_triples[1,1])){ #replace 1st row if need be
          ind <- 1
        }else{
          ind <- nrow(unsh_triples) + 1 #otherwise add to end of df
        }

        unsh_triples[ind,] <- c( #store in df
          #Xi_tau
          Xi,
          #Xk_t
          Xk_t,
          #Xj_t
          Xj_t,
          #ambiguous
          0,
          #flag
          0
        )
      }

      #Reattempt with alternative Xk/Xj assignment
      Xk_t <- contemp_adj$Xi[row]
      Xj_t <- contemp_adj$Xj[row]

      #Is there any Xi_tau that has an arrow or edge into Xk_t but not adj to Xj_t?

      Xi_tau <- G |>
        filter(Xj == Xk_t | Xi == Xk_t) |>
        filter(Xj != Xj_t & Xi != Xj_t) |>
        mutate(Xk_mark = case_when(Xj == Xk_t ~ Xj_mark,
                                   Xi == Xk_t ~ Xi_mark)) |>
        filter(Xk_mark == "<") |> #arrow or edge into Xk_t?
        mutate(other = ifelse(Xj == Xk_t, Xi, Xj)) |> #all other values
        pull(other)

      if(length(Xi_tau) == 0){
        next
      }

      Xi_tau <- G |> #no arrow or edge into Xj_t?
        filter(Xj == Xj_t | Xi == Xj_t) |>
        filter(is.na(Xj_mark)) |> #only those with no edge
        mutate(other = ifelse(Xj == Xj_t, Xi, Xj)) |> #all other values
        filter(other %in% Xi_tau) |>  #is Xi_tau in this list?
        pull(other)

      if(length(Xi_tau) == 0){
        next
      }

      #there are instance(s) of an unshielded triple, add to df
      for(Xi in Xi_tau){

        if(is.na(unsh_triples[1,1])){ #replace 1st row if need be
          ind <- 1
        }else{
          ind <- nrow(unsh_triples) + 1 #otherwise add to end of df
        }

        unsh_triples[ind,] <- c( #store in df
          #Xi_tau
          Xi,
          #Xk_t
          Xk_t,
          #Xj_t
          Xj_t,
          #ambiguous
          0,
          #flag
          0
        )
      }
    }

    #if no unsh_triples were found at all, break
    if(is.na(unsh_triples[1,1])){
      testing <- FALSE
      break
    }

    ###
    # Test unshielded triples for colliders
    ###

    for(row in seq(1, nrow(unsh_triples))){

      Xi_tau <- unsh_triples$Xi_tau[row]
      i <- ((Xi_tau - 1) %/% (1 + 2 * max_lag)) + 1

      Xk_t <- unsh_triples$Xk_t[row]

      Xj_t <- unsh_triples$Xj_t[row]
      j <- which(X_t == Xj_t)

      if(collider_rule == "none"){

        sep_set <- sep_sets |>
          filter((Xi == Xi_tau & Xj == Xj_t) | (Xi == Xj_t & Xj == Xi_tau)) |>
          pull(S_sep)

        if(length(sep_set) > 0){
          #If Xk_t is not in the separating set
          if(Xk_t %in% suppressWarnings(as.numeric(str_split(sep_set, ";", simplify = TRUE))) == FALSE){

            #Orient Xj_t -> Xk_t if sufficient
            if(sufficiency){
              G <- G |>
                mutate(
                  Xj_mark = case_when(
                    Xj == Xj_t & Xi == Xk_t ~ "-",
                    Xj == Xk_t & Xi == Xj_t ~ "<",
                    TRUE ~ Xj_mark  # retain existing value
                  ),
                  Xi_mark = case_when(
                    Xj == Xj_t & Xi == Xk_t ~ "<",
                    Xj == Xk_t & Xi == Xj_t ~ "-",
                    TRUE ~ Xi_mark  # retain existing value
                  )
                )
            }else{ #Or just add arrowhead into Xk_t if not sufficient
              G <- G |>
                mutate(
                  Xj_mark = case_when(
                    Xj == Xk_t & Xi == Xj_t ~ "<",
                    TRUE ~ Xj_mark  # retain existing value
                  ),
                  Xi_mark = case_when(
                    Xj == Xj_t & Xi == Xk_t ~ "<",
                    TRUE ~ Xi_mark  # retain existing value
                  )
                )

            }

            #And flag to remove from unsh_triples df
            unsh_triples$flag[row] <- 1
          }
        }

        if(row == nrow(unsh_triples)){
          break
        }
        next
      }

      Aj_t <- Aj_ts[[j]]
      if(Xi_tau %in% X_t){ #include contemp set of i if relevant
        Ai_t <- Aj_ts[[i]]
        C_set <- setdiff(unique(c(Aj_t, Ai_t)), Xi_tau)
      }else{
        C_set <- setdiff(Aj_t, Xi_tau)
      }

      C_set <- as.vector(C_set)

      #All subsets
      n <- length(C_set)
      S_s <- list()

      if (n > 0) {
        for (w in 1:(2^n - 1)) {  # skip 0 to avoid empty set
          mask <- as.logical(intToBits(w))[1:n]
          S_s[[length(S_s) + 1]] <- C_set[mask]
        }
      }

      for(S in S_s){

        #Set Z as the set containing S, lagged adj of Xj - Xi_tau, and lagged adj of Xi_tau
        Bj_t <- setdiff(Bj_ts[[j]], Xi_tau)

        #if tau = 0
        if(Xi_tau %in% X_t){
          Bi_t <- Bj_ts[[i]]
          #or if greater than 0
        }else{
          lag <- Xi_tau - X_t[i]
          Bi_t <- Bj_ts[[i]]  + lag
        }

        Z <- unique(c(S, Bj_t, Bi_t))

        if(Xi_tau %in% Z || Xj_t %in% Z){
          stop("Self model in collider testing")
        }

        if(!silent){cat("\nFitting ", colnames(design)[Xi_tau], " ~ ", Z)}

        #Test Xi_tau _||_ Xj_tau | Z:
        #Fit and extract residuals for a model of Xi_past
        Xi_tau_Z <- find_or_fit(design = design,
                                response = Xi_tau,
                                model_on = Z,
                                distribution = distributions[Xi_tau],
                                Env = res_store,
                                modeltype = modeltype,
                                spatial_str = space)

        if(!silent){cat("\nFitting ", colnames(design)[Xi_tau], " ~ ", Z)}

        #Fit and extract residuals for Xj_t | Z
        Xj_t_Z <- find_or_fit(design = design,
                              response = Xj_t,
                              model_on = Z,
                              distribution = distributions[Xj_t],
                              Env = res_store,
                              modeltype = modeltype,
                              spatial_str = space)

        if(independence_test == "GCM"){
          test_result <-  gcm.test(resid.XonZ = Xi_tau_Z,
                                   resid.YonZ =  Xj_t_Z)
        }else{
          if(is.character(Z)){
            zed <- as.numeric(design[[sites]])
          }else{
            zed <- design[,Z]
          }
          test_result <- WGCM_fix(resid.XonZ = Xi_tau_Z,
                                  resid.YonZ =  Xj_t_Z,
                                  Z = zed,
                                  k0 = k0,
                                  nsim = nsim)
        }

        if(test_result$p.value > alpha){ #if independent

          #store as sep. set
          store_ind[nrow(store_ind) + 1,] <- c(
            #Xj
            Xj_t,
            #Xi
            Xi_tau,
            #Imin
            test_result$test.statistic,
            #S
            paste(S, collapse = ";")
          )

        }

      }

      #list of all separating sets
      all_sets <- store_ind |>
        filter(Xj == Xj_t & Xi == Xi_tau) |>
        pull(S)

      #If no separating sets found, mark as ambiguous:
      if(length(all_sets) == 0){
        unsh_triples$ambiguous[row] <- 1
      }else{ #If some found:

        n_k <- mean(grepl(paste0("(?<!\\d)", Xk_t, "(?!\\d)"),
                          unique(all_sets), perl=TRUE))

        if(collider_rule == "conservative"){

          if(n_k == 0){
            #Orient Xj_t -> Xk_t if sufficient
            if(sufficiency){
              G <- G |>
                mutate(
                  Xj_mark = case_when(
                    Xj == Xj_t & Xi == Xk_t ~ "-",
                    Xj == Xk_t & Xi == Xj_t ~ "<",
                    TRUE ~ Xj_mark  # retain existing value
                  ),
                  Xi_mark = case_when(
                    Xj == Xj_t & Xi == Xk_t ~ "<",
                    Xj == Xk_t & Xi == Xj_t ~ "-",
                    TRUE ~ Xi_mark  # retain existing value
                  )
                )
            }else{ #Or just add arrowhead into Xk_t if not sufficient
              G <- G |>
                mutate(
                  Xj_mark = case_when(
                    Xj == Xk_t & Xi == Xj_t ~ "<",
                    TRUE ~ Xj_mark  # retain existing value
                  ),
                  Xi_mark = case_when(
                    Xj == Xj_t & Xi == Xk_t ~ "<",
                    TRUE ~ Xi_mark  # retain existing value
                  )
                )

            }

            #And flag to remove from unsh_triples df
            unsh_triples$flag[row] <- 1

          }else{
            if(n_k == 1){
              #leave as is
              next
            }else{
              #if any other value, mark as ambiguous
              unsh_triples$ambiguous[row] <- 1
            }
          }
        }else{ #if collider_rule is majority or unset
          if(n_k < 0.5){

            #Orient Xj_t -> Xk_t if sufficient
            if(sufficiency){
              G <- G |>
                mutate(
                  Xj_mark = case_when(
                    Xj == Xj_t & Xi == Xk_t ~ "-",
                    Xj == Xk_t & Xi == Xj_t ~ "<",
                    TRUE ~ Xj_mark  # retain existing value
                  ),
                  Xi_mark = case_when(
                    Xj == Xj_t & Xi == Xk_t ~ "<",
                    Xj == Xk_t & Xi == Xj_t ~ "-",
                    TRUE ~ Xi_mark  # retain existing value
                  )
                )
            }else{ #Or just add arrowhead into Xk_t if not sufficient
              G <- G |>
                mutate(
                  Xj_mark = case_when(
                    Xj == Xk_t & Xi == Xj_t ~ "<",
                    TRUE ~ Xj_mark  # retain existing value
                  ),
                  Xi_mark = case_when(
                    Xj == Xj_t & Xi == Xk_t ~ "<",
                    TRUE ~ Xi_mark  # retain existing value
                  )
                )

            }

            #And flag to remove from unsh_triples df
            unsh_triples$flag[row] <- 1

          }else{
            if(n_k > 0.5){
              #leave as is
              next
            }else{
              #if = 0.5 or NA, mark as ambiguous
              unsh_triples$ambiguous[row] <- 1
            }
          }
        }
      }
    }

    #Once all unsh_triples tested, remove those flagged
    unsh_triples <- unsh_triples |>
      filter(flag == 0)

    testing <- FALSE
  }

  if(skeleton_names == TRUE){
    #replace indices with names for interpretability
    G <- G |>
      mutate(Xj = colnames(design)[Xj],
             Xi = colnames(design)[Xi])

    list2return <- list(skeleton = G, independencies = store_ind, triples = unsh_triples, sep_sets = sep_sets)
  }else{
    list2return <- list(skeleton = G, independencies = store_ind, triples = unsh_triples, names = colnames(design), sep_sets = sep_sets)
  }

  return(list2return)

}
