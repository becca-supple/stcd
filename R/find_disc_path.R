#' Find discriminating paths as part of edge orientation
#
#' @param graph Graph object as produced within [orient_edges]
#' @param i Indices of nodes to check for discriminating paths between
#' @param k Indices of nodes to check for discriminating paths between
#' @param j Indices of nodes to check for discriminating paths between
#'
#' @returns A vector of node indices forming the discriminating path, or `NA` if
#' no path is found
#' @export
#'
find_disc_path <- function(graph, i, k, j){

  #set up path p
  p <- NA

  #initiate search length
  search_length <- 0

  #Do any variables other than Xk have an edge with an arrowhead into Xi?
  Xl <- graph |>
    filter(Xj == i | Xi == i) |>
    mutate(mark = case_when(Xj == i ~ Xj_mark,
                            Xi == i ~ Xi_mark)) |>
    filter(mark == "<") |>
    mutate(other = case_when(Xj == i ~ Xi,
                             Xi == i ~ Xj)) |>
    filter(other != k) |>
    pull(other)

  if(length(Xl) > 0){

    while(any(is.na(p))){

      for(l in Xl){

        #Is the new node adjacent to Xj?
        Xl.j_mark <- graph |>
          filter((Xi == l & Xj == j) |
                   (Xi == j & Xj == l)) |>
          mutate(mark = case_when(Xj == j ~ Xj_mark,
                                  Xi == j ~ Xi_mark)) |>
          pull(mark)

        if(is.na(Xl.j_mark)){

          p <- c(l, i, k, j)
          break

        }else{

          #is the arrow from the new node directed into Xj
          if(Xl.j_mark == "<"){

            Xj.l_mark <- graph |>
              filter((Xi == l & Xj == j) |
                       (Xi == j & Xj == l)) |>
              mutate(mark = case_when(Xj == l ~ Xj_mark,
                                      Xi == l ~ Xi_mark)) |>
              pull(mark)

            if(Xj.l_mark == "-"){

              #is there an arrowhead into the new node from Xi?
              Xi.l_mark <- graph |>
                filter((Xi == l & Xj == i) |
                         (Xi == i & Xj == l)) |>
                mutate(mark = case_when(Xj == l ~ Xj_mark,
                                        Xi == l ~ Xi_mark)) |>
                pull(mark)

              if(Xi.l_mark == "<"){

                #Then this node could be an intermediate on a disc path
                #Need to search for nodes with arrows into this node:
                Xm <- graph |>
                  filter(Xj == l | Xi == l) |>
                  mutate(mark = case_when(Xj == l ~ Xj_mark,
                                          Xi == l ~ Xi_mark)) |>
                  filter(mark == "<") |>
                  mutate(other = case_when(Xj == l ~ Xi,
                                           Xi == l ~ Xj)) |>
                  filter(other != i) |>
                  pull(other)

                if(length(Xm) > 0){

                  for(m in Xm){

                    #Is the new node adjacent to Xj?
                    Xm.j_mark <- graph |>
                      filter((Xi == m & Xj == j) |
                               (Xi == j & Xj == m)) |>
                      mutate(mark = case_when(Xj == j ~ Xj_mark,
                                              Xi == j ~ Xi_mark)) |>
                      pull(mark)

                    if(is.na(Xm.j_mark)){

                      p <- c(m, l, i, k, j)
                      break

                    }else{

                      #is the arrow from the new node directed into Xj
                      if(Xm.j_mark == "<"){

                        Xj.m_mark <- graph |>
                          filter((Xi == m & Xj == j) |
                                   (Xi == j & Xj == m)) |>
                          mutate(mark = case_when(Xj == m ~ Xj_mark,
                                                  Xi == m ~ Xi_mark)) |>
                          pull(mark)

                        if(Xj.m_mark == "-"){

                          #is there an arrowhead into the new node from Xl?
                          Xl.m_mark <- graph |>
                            filter((Xi == m & Xj == l) |
                                     (Xi == l & Xj == m)) |>
                            mutate(mark = case_when(Xj == m ~ Xj_mark,
                                                    Xi == m ~ Xi_mark)) |>
                            pull(mark)

                          if(Xl.m_mark == "<"){

                            #Then this node could be an intermediate on a disc path
                            #Need to search for nodes with arrows into this node:
                            Xn <- graph |>
                              filter(Xj == m | Xi == m) |>
                              mutate(mark = case_when(Xj == m ~ Xj_mark,
                                                      Xi == m ~ Xi_mark)) |>
                              filter(mark == "<") |>
                              mutate(other = case_when(Xj == m ~ Xi,
                                                       Xi == m ~ Xj)) |>
                              filter(other != l) |>
                              pull(other)

                            if(length(Xn) > 0){

                              for(n in Xn){

                                #Is the new node adjacent to Xj?
                                Xn.j_mark <- graph |>
                                  filter((Xi == n & Xj == j) |
                                           (Xi == j & Xj == n)) |>
                                  mutate(mark = case_when(Xj == j ~ Xj_mark,
                                                          Xi == j ~ Xi_mark)) |>
                                  pull(mark)

                                if(is.na(Xn.j_mark)){

                                  p <- c(n, m, l, i, k, j)
                                  break

                                }else{

                                  #is the arrow from the new node directed into Xj
                                  if(Xn.j_mark == "<"){

                                    Xj.n_mark <- graph |>
                                      filter((Xi == n & Xj == j) |
                                               (Xi == j & Xj == n)) |>
                                      mutate(mark = case_when(Xj == n ~ Xj_mark,
                                                              Xi == n ~ Xi_mark)) |>
                                      pull(mark)

                                    if(Xj.n_mark == "-"){

                                      #is there an arrowhead into the new node from Xm?
                                      Xm.n_mark <- graph |>
                                        filter((Xi == n & Xj == m) |
                                                 (Xi == m & Xj == n)) |>
                                        mutate(mark = case_when(Xj == n ~ Xj_mark,
                                                                Xi == n ~ Xi_mark)) |>
                                        pull(mark)

                                      if(Xm.n_mark == "<"){

                                        #Then this node could be an intermediate on a disc path
                                        #Need to search for nodes with arrows into this node:
                                        Xo <- graph |>
                                          filter(Xj == n | Xi == n) |>
                                          mutate(mark = case_when(Xj == n ~ Xj_mark,
                                                                  Xi == n ~ Xi_mark)) |>
                                          filter(mark == "<") |>
                                          mutate(other = case_when(Xj == n ~ Xi,
                                                                   Xi == n ~ Xj)) |>
                                          filter(other != m) |>
                                          pull(other)

                                        if(length(Xo) > 0){

                                          for(o in Xo){

                                            #Is the new node adjacent to Xj?
                                            Xo.j_mark <- graph |>
                                              filter((Xi == o & Xj == j) |
                                                       (Xi == j & Xj == o)) |>
                                              mutate(mark = case_when(Xj == j ~ Xj_mark,
                                                                      Xi == j ~ Xi_mark)) |>
                                              pull(mark)

                                            if(is.na(Xo.j_mark)){

                                              p <- c(o, n, m, l, i, k, j)
                                              break

                                            }else{

                                              search_length <- search_length + 1

                                              warning("\nMaximum search length reached without eliminating potential of discriminating path containing p = <..., ",
                                                      o, ", ",  n, ", ", m, ", ", l, ", ", i, ", ",  k, ", ", j,
                                                      ">. Manual search suggested.")

                                              if(search_length == length(Xl) * length(Xm) * length(Xn) * length(Xo)){

                                                return(NA)

                                              }
                                            }
                                          }
                                        }
                                      }
                                    }
                                  }
                                }
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  return(p)

}
