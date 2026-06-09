#' Build neighborhood object from a dataframe of contiguity information
#'
#' @param site_names Vector of site names
#' @param contiguity_df data.frame with two columns for all site combinations
#' and a column indicating binary contiguity with 0/1
#' @param site1 String indicating name of the column with site 1
#' @param site2 String indicating name of the column with site 2
#' @param contiguity_col String indicating the name of the column with
#' contiguity information
#'
#' @returns A named list of vectors, one for each site, which lists all its
#' neighboring sites
#'
#' NOTE: For use in [fit_custom_gam], all sites must have at least one neighbor
#'
#' @export
#'
#' @examples
#'
#' #Set up contiguity df for the boroughs of NYC
#' boroughs <- c("Manhattan", "The Bronx", "Queens", "Brooklyn", "Staten Island")
#' nyc_cont <- as.data.frame(t(combn(boroughs, 2, simplify = TRUE)))
#' nyc_cont$neighbor <- c(rep(0, 7), 1, rep(0, 2))
#'
#' head(nyc_cont)
#'
#' #Build neighborhood
#' nyc_nb <- build_neighborhood(site_names = boroughs,
#'               contiguity_df = nyc_cont,
#'               site1 = "V1",
#'               site2 = "V2",
#'               contiguity_col = "neighbor")
#'
#' nyc_nb
build_neighborhood <- function(site_names,
                               contiguity_df, #data frame with all site combinations and a column that == 1 if contiguous
                               site1 = "state1ab", #name of column with site 1
                               site2 = "state2ab", #name of column with site 2
                               contiguity_col = "conttype" #name of column with contiguity info
){

  nb <- list()

  #Filter to just contiguities
  contiguity_df <- contiguity_df[contiguity_df[[contiguity_col]] == 1,]

  for(c in site_names){

    #find all contiguous other sites
    connections <- contiguity_df |>
      dplyr::mutate(connections = ifelse(contiguity_df[[site1]] == c, paste(contiguity_df[[site2]]),
                                  ifelse(contiguity_df[[site2]] == c, paste(contiguity_df[[site1]]), NA))) |>
      dplyr::select(connections)

    #filter to sites we're interested in
    connections <- unique(connections$connections[is.na(connections$connections) == F])
    connections <- connections[connections %in% site_names]

    nb[[c]] <- connections
  }

  return(nb)

}
