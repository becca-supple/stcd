#' Orient edges
#'
#' @param independencies Outputs from conditional independence algorithm [conditional_independencies]
#' @param max_lag Maximum considered time lag
#' @param names If = TRUE, replace indices with variable names in output
#' @param sufficiency FALSE by default
#'
#' @returns A PAG object if sufficiency = FALSE, a CPDAG if TRUE
#' @export
#'
#' @examples
#'
#' #' head(example_coord_data)
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
#'
#' example_or <- orient_edges(independencies = example_inds,
#'                 max_lag = 2,
#'                 names = TRUE,
#'                 sufficiency = FALSE
#'                 )
#'
#' str(example_or)
#'
#' }
orient_edges <- function(independencies,
                         max_lag = 2,
                         names = TRUE,
                         sufficiency = FALSE
){

  G <- independencies$skeleton
  trips <- independencies$triples
  X <- seq(1, (length(independencies$names) - 1)/ (2 * max_lag + 1))
  X_t <- 1 + (1 + 2*max_lag)*(X - 1)


  ###
  # ID unambiguous triples
  ###

  #Start with unambiguous triples ID'd in last algorithm

  unamb_trips <- trips |>
    filter(ambiguous == 0) |>
    mutate(cycle = 0)

  #Add triples that form a potential cycle

  contemp_adj <- G |>
    filter(is.na(Xj_mark) == FALSE) |>
    filter(Xi %in% X_t & Xj %in% X_t)

  if(nrow(contemp_adj) > 0){ #if there are contemporaneous adjacencies

    for(row in seq(1, nrow(contemp_adj))){

      #trying with this orientation of Xj/Xi first, then will do opposite
      Xk_t <- contemp_adj$Xj[row]
      Xj_t <- contemp_adj$Xi[row]

      #Is there any Xi_tau that has an arrow or edge into Xk_t and adj to Xj_t?

      Xi_tau <- G |> #arrow or edge into Xk_t?
        filter(Xj == Xk_t | Xi == Xk_t) |>
        filter(is.na(Xj_mark) == FALSE) |> #any edge
        mutate(other = ifelse(Xj == Xk_t, Xi, Xj)) |> #all other values
        pull(other)

      Xi_tau <- G |> #arrow or edge into Xj_t?
        filter(Xj == Xj_t | Xi == Xj_t) |>
        filter(is.na(Xj_mark) == FALSE) |> #any edge
        mutate(other = ifelse(Xj == Xj_t, Xi, Xj)) |> #all other values
        filter(other %in% Xi_tau) |>  #is Xi_tau in this list?
        pull(other)

      if(length(Xi_tau) == 0){
        next
      }

      #if there are relevant Xi_tau, add to the triples df
      for(Xi in Xi_tau){

        unamb_trips[nrow(unamb_trips) + 1,] <- c( #store in df
          #Xi_tau
          Xi,
          #Xk_t
          Xk_t,
          #Xj_t
          Xj_t,
          #ambiguous
          0,
          #flag
          0,
          #cycle
          1
        )
      }

      #Swap Xk_t and Xj_t
      Xk_t <- contemp_adj$Xi[row]
      Xj_t <- contemp_adj$Xj[row]

      #Is there any Xi_tau that has an arrow or edge into Xk_t and adj to Xj_t?

      Xi_tau <- G |> #arrow or edge into Xk_t?
        filter(Xj == Xk_t | Xi == Xk_t) |>
        filter(is.na(Xj_mark) == FALSE) |> #any edge
        mutate(other = ifelse(Xj == Xk_t, Xi, Xj)) |> #all other values
        pull(other)

      Xi_tau <- G |> #arrow or edge into Xj_t?
        filter(Xj == Xj_t | Xi == Xj_t) |>
        filter(is.na(Xj_mark) == FALSE) |> #any edge
        mutate(other = ifelse(Xj == Xj_t, Xi, Xj)) |> #all other values
        filter(other %in% Xi_tau) |>  #is Xi_tau in this list?
        pull(other)

      if(length(Xi_tau) == 0){
        next
      }

      #if there are relevant Xi_tau, add to the triples df
      for(Xi in Xi_tau){

        unamb_trips[nrow(unamb_trips) + 1,] <- c( #store in df
          #Xi_tau
          Xi,
          #Xk_t
          Xk_t,
          #Xj_t
          Xj_t,
          #ambiguous
          0,
          #flag
          0,
          #cycle
          1
        )
      }
    }
  }

  #Make sure each row of unamb_trips is unique
  unamb_trips <- unamb_trips |>
    filter_out(duplicated(unamb_trips))

  while(nrow(unamb_trips) > 0){

    for(tr in seq(1, nrow(unamb_trips))){

      Xi_tau <- unamb_trips$Xi_tau[tr]
      Xk_t <- unamb_trips$Xk_t[tr]
      Xj_t <- unamb_trips$Xj_t[tr]

      ###
      # Rule 1: Non-colliders with an arrow into the middle must be chains
      ###

      if(unamb_trips$cycle[tr] != 1){

        Xi.k_mark <- G |>
          filter(Xi == Xi_tau & Xj == Xk_t |
                   Xi == Xk_t & Xj == Xi_tau) |>
          mutate(mark = case_when(Xi == Xk_t ~ Xi_mark,
                                  Xj == Xk_t ~ Xj_mark)) |>
          pull(mark)

        if(Xi.k_mark == "<"){

          Xj.k_mark <- G |>
            filter(Xi == Xj_t & Xj == Xk_t |
                     Xi == Xk_t & Xj == Xj_t) |>
            mutate(mark = case_when(Xi == Xk_t ~ Xi_mark,
                                    Xj == Xk_t ~ Xj_mark)) |>
            pull(mark)

          if(Xj.k_mark == "o"){

            G <- G |>
              mutate(Xj_mark = replace_when(Xj_mark,
                                            Xj == Xk_t & Xi == Xj_t ~ "-",
                                            Xj == Xj_t & Xi == Xk_t ~ "<"),
                     Xi_mark = replace_when(Xi_mark,
                                            Xj == Xk_t & Xi == Xj_t ~ "<",
                                            Xj == Xj_t & Xi == Xk_t ~ "-"))

          }

        }

        #Flag for removal

        unamb_trips$flag[tr] <- 1

      }

      ###
      # Rule 2: Avoid cycles
      ###

      if(unamb_trips$cycle[tr] == 1){ #For a potential cycle

        Xi.k_mark <- G |>
          filter(Xi == Xi_tau & Xj == Xk_t |
                   Xi == Xk_t & Xj == Xi_tau) |>
          mutate(mark = case_when(Xi == Xk_t ~ Xi_mark,
                                  Xj == Xk_t ~ Xj_mark)) |>
          pull(mark)

        if(Xi.k_mark == "<"){ #if arrowhead into Xk

          Xk.j_mark <- G |>
            filter(Xi == Xj_t & Xj == Xk_t |
                     Xi == Xk_t & Xj == Xj_t) |>
            mutate(mark = case_when(Xi == Xj_t ~ Xi_mark,
                                    Xj == Xj_t ~ Xj_mark)) |>
            pull(mark)

          if(Xk.j_mark == "<"){ #and into Xj

            Xk.i_mark <- G |>
              filter(Xi == Xi_tau & Xj == Xk_t |
                       Xi == Xk_t & Xj == Xi_tau) |>
              mutate(mark = case_when(Xi == Xi_tau ~ Xi_mark,
                                      Xj == Xi_tau ~ Xj_mark)) |>
              pull(mark)

            Xj.k_mark <- G |>
              filter(Xi == Xj_t & Xj == Xk_t |
                       Xi == Xk_t & Xj == Xj_t) |>
              mutate(mark = case_when(Xi == Xk_t ~ Xi_mark,
                                      Xj == Xk_t ~ Xj_mark)) |>
              pull(mark)

            if((sufficiency && Xk.i_mark == "-" && Xj.k_mark == "-") || #if fully directed arrows from Xi->Xk->Xj & sufficiency assumed
               (!sufficiency && Xk.i_mark == "-") || #or if not sufficient but we have Xi->Xk
               (!sufficiency && Xj.k_mark == "-")){ #or if not sufficient but we have Xk->Xj

              Xi.j_mark <- G |>
                filter(Xi == Xj_t & Xj == Xi_tau |
                         Xi == Xi_tau & Xj == Xj_t) |>
                mutate(mark = case_when(Xi == Xj_t ~ Xi_mark,
                                        Xj == Xj_t ~ Xj_mark)) |>
                pull(mark)

              Xj.i_mark <- G |>
                filter(Xi == Xj_t & Xj == Xi_tau |
                         Xi == Xi_tau & Xj == Xj_t) |>
                mutate(mark = case_when(Xi == Xi_tau ~ Xi_mark,
                                        Xj == Xi_tau ~ Xj_mark)) |>
                pull(mark)

              if((sufficiency && Xi.j_mark == "o" && Xj.i_mark == "o")){ #if sufficiency assumed and Xi o-o Xj

                G <- G |>  #orient Xi --> Xj
                  mutate(Xj_mark = replace_when(Xj_mark,
                                                Xj == Xi_tau & Xi == Xj_t ~ "-",
                                                Xj == Xj_t & Xi == Xi_tau ~ "<"),
                         Xi_mark = replace_when(Xi_mark,
                                                Xj == Xi_tau & Xi == Xj_t ~ "<",
                                                Xj == Xj_t & Xi == Xi_tau ~ "-"))

                #Flag for removal

                unamb_trips$flag[tr] <- 1

              }else{

                if(!sufficiency){ #if no assumed sufficiency

                  if(Xi.j_mark == "o"){ #and unresolved mark into Xj

                    G <- G |>  #orient Xi *-> Xj
                      mutate(Xj_mark = replace_when(Xj_mark,
                                                    Xj == Xj_t & Xi == Xi_tau ~ "<"),
                             Xi_mark = replace_when(Xi_mark,
                                                    Xj == Xi_tau & Xi == Xj_t ~ "<"))

                    #Flag for removal

                    unamb_trips$flag[tr] <- 1

                  }

                }

              }

            }

          }

        }

      }

      ###
      # Rule 3: Orient central edges in double cycles
      ###

      if(unamb_trips$cycle[tr] == 1){

        #Are there any adjacent cycles?
        adj_cycles <- unamb_trips |>
          filter((Xi_tau == Xi_tau & Xj_t == Xj_t & Xk_t != Xk_t) |
                   (Xi_tau == Xj_t & Xi_tau == Xj_t & Xk_t != Xk_t))

        if(nrow(adj_cycles) > 0){

          for(l in seq(1, nrow(adj_cycles))){

            Xl_t <- adj_cycles$Xk_t[l]

            #Are there arrowheads into Xj from Xk?
            Xk.j_mark <- G |>
              filter((Xi == Xj_t & Xj == Xk_t) |
                       (Xi == Xk_t & Xj == Xj_t)) |>
              mutate(mark = case_when(Xj == Xj_t ~ Xj_mark,
                                      Xi == Xj_t ~ Xi_mark)) |>
              pull(mark)

            if(Xk.j_mark == "<"){

              #And Xl?
              Xl.j_mark <- G |>
                filter((Xi == Xj_t & Xj == Xl_t) |
                         (Xi == Xl_t & Xj == Xj_t)) |>
                mutate(mark = case_when(Xj == Xj_t ~ Xj_mark,
                                        Xi == Xj_t ~ Xi_mark)) |>
                pull(mark)

              if(Xl.j_mark == "<"){

                #Are there unresolved marks at Xi from Xl?
                Xl.i_mark <- G |>
                  filter((Xi == Xi_tau & Xj == Xl_t) |
                           (Xi == Xl_t & Xj == Xi_tau)) |>
                  mutate(mark = case_when(Xj == Xi_tau ~ Xj_mark,
                                          Xi == Xi_tau ~ Xi_mark)) |>
                  pull(mark)

                if(Xl.i_mark == "o"){

                  #And from Xk?
                  Xk.i_mark <- G |>
                    filter((Xi == Xi_tau & Xj == Xk_t) |
                             (Xi == Xk_t & Xj == Xi_tau)) |>
                    mutate(mark = case_when(Xj == Xi_tau ~ Xj_mark,
                                            Xi == Xi_tau ~ Xi_mark)) |>
                    pull(mark)

                  if(Xk.i_mark == "o"){

                    #Is the mark Xi -o Xj unresolved?
                    Xi.j_mark <- G |>
                      filter((Xi == Xj_t & Xj == Xi_tau) |
                               (Xi == Xi_tau & Xj == Xj_t)) |>
                      mutate(mark = case_when(Xj == Xj_t ~ Xj_mark,
                                              Xi == Xj_t ~ Xi_mark)) |>
                      pull(mark)

                    if(Xi.j_mark == "o"){

                      if(sufficiency){ #If assuming sufficiency, orient Xi --> Xj

                        G <- G |>
                          mutate(Xj_mark = replace_when(Xj_mark,
                                                        Xj == Xj_t & Xi == Xi_tau ~ "<",
                                                        Xj == Xi_tau & Xi == Xj_t ~ "-"),
                                 Xi_mark = replace_when(Xi_mark,
                                                        Xj == Xi_tau & Xi == Xj_t ~ "<",
                                                        Xj == Xj_t & Xi == Xi_tau ~ "-"))

                      }

                      if(!sufficiency){#If not assuming sufficiency, change mark at Xj from Xi to "<"

                        G <- G |>
                          mutate(Xj_mark = replace_when(Xj_mark,
                                                        Xj == Xj_t & Xi == Xi_tau ~ "<"),
                                 Xi_mark = replace_when(Xi_mark,
                                                        Xj == Xi_tau & Xi == Xj_t ~ "<"))

                      }

                      #Flag for removal

                      unamb_trips$flag[tr] <- 1

                    }

                  }

                }

              }

            }


          }

        }

      }

    }

    if(!sufficiency){

      ###
      # Rule 4: Discriminating paths orient bidirected edges
      ###

      #If shielded
      if(unamb_trips$cycle[tr] == 1){

        #Is there a circle mark at Xk from Xj?
        Xj.k_mark <- G |>
          filter((Xi == Xj_t & Xj == Xk_t) |
                   (Xi == Xk_t & Xj == Xj_t)) |>
          mutate(mark = case_when(Xj == Xk_t ~ Xj_mark,
                                  Xi == Xk_t ~ Xi_mark)) |>
          pull(mark)

        if(Xj.k_mark == "o"){

          #Is there an edge Xi -> Xj?
          Xi.j_mark <- G |>
            filter((Xi == Xj_t & Xj == Xi_tau) |
                     (Xi == Xi_tau & Xj == Xj_t)) |>
            mutate(mark = case_when(Xj == Xj_t ~ Xj_mark,
                                    Xi == Xj_t ~ Xi_mark)) |>
            pull(mark)

          if(Xi.j_mark == "<"){

            Xj.i_mark <- G |>
              filter((Xi == Xj_t & Xj == Xi_tau) |
                       (Xi == Xi_tau & Xj == Xj_t)) |>
              mutate(mark = case_when(Xj == Xi_tau ~ Xj_mark,
                                      Xi == Xi_tau ~ Xi_mark)) |>
              pull(mark)

            if(Xj.i_mark == "-"){

              #Is there an arrowhead into Xi from Xk?
              Xk.i_mark <- G |>
                filter((Xi == Xi_tau & Xj == Xk_t) |
                         (Xi == Xk_t & Xj == Xi_tau)) |>
                mutate(mark = case_when(Xj == Xi_tau ~ Xj_mark,
                                        Xi == Xi_tau ~ Xi_mark)) |>
                pull(mark)

              if(Xk.i_mark == "<"){

                #Then search for discriminating path
                p <- find_disc_path(graph = G, i = Xi_tau, k = Xk_t, j = Xj_t)

                #If a path exists
                if(any(is.na(p)) == FALSE){

                  #Is Xk in the separating set for the first node in p and Xj?
                  all_sets <- independencies$independencies |>
                    filter(Xj == Xj_t & Xi == p[1] |
                             Xj == p[1] & Xi == Xj_t) |>
                    pull(S)

                  if(any(grepl(paste0("(?<!\\d)", Xk_t, "(?!\\d)"),
                               unique(all_sets), perl=TRUE))){

                    #If yes, orient Xk -> Xj
                    G <- G |>
                      mutate(Xj_mark = replace_when(Xj_mark,
                                                    Xj == Xj_t & Xi == Xk_t ~ "<",
                                                    Xj == Xk_t & Xi == Xj_t ~ "-"),
                             Xi_mark = replace_when(Xi_mark,
                                                    Xj == Xk_t & Xi == Xj_t ~ "<",
                                                    Xj == Xj_t & Xi == Xk_t ~ "-"))

                    #and flag for removal
                    unamb_trips$flag[tr] <- 1

                  }else{

                    #If no, orient Xi <-> Xk <-> Xj

                    G <- G |>
                      mutate(Xj_mark = replace_when(Xj_mark,
                                                    Xj == Xj_t & Xi == Xk_t ~ "<",
                                                    Xj == Xk_t & Xi == Xj_t ~ "<"),
                             Xi_mark = replace_when(Xi_mark,
                                                    Xj == Xk_t & Xi == Xj_t ~ "<",
                                                    Xj == Xj_t & Xi == Xk_t ~ "<"))

                    G <- G |>
                      mutate(Xj_mark = replace_when(Xj_mark,
                                                    Xj == Xk_t & Xi == Xi_tau ~ "<",
                                                    Xj == Xi_tau & Xi == Xk_t ~ "<"),
                             Xi_mark = replace_when(Xi_mark,
                                                    Xj == Xi_tau & Xi == Xk_t ~ "<",
                                                    Xj == Xk_t & Xi == Xi_tau ~ "<"))

                    #and flag for removal
                    unamb_trips$flag[tr] <- 1

                  }

                }

              }

            }

          }

        }

      }

    }

    #Remove any flagged entries
    unamb_trips_new <- unamb_trips |>
      filter(flag == 0)

    #If no change in the df, break the loop
    if(identical(unamb_trips_new, unamb_trips)){
      unamb_trips <- matrix(nrow = 0, ncol = 1)
    }else{ #if changes, update to newer version
      unamb_trips <- unamb_trips_new
    }
  }

  ###
  # Format and return
  ###

  #Replace indices with variable names
  if(names){
    G <- G |>
      mutate(Xj = independencies$names[Xj],
             Xi = independencies$names[Xi])
  }

  if(sufficiency){
    G <- CPDAG(G)
  }else{
    G <- PAG(G)
  }

  return(G)

}
