#' Example coordinate space design matrix
#'
#' A simulated `coord_dsgn` object with shifted time series of agricultural
#' activity and probability of a species' occupancy at randomly sampled sites
#' across Great Britain.
#'
#' @format ## `example_coord_data`
#' A data frame with 5,226 rows and 22 columns:
#' \describe{
#'   \item{Agr, Agr_tminus1, ..., Agr_tminus4}{% of land used for agriculture}
#'   \item{Fence, Fence_tminus1, ..., Fence_tminus4}{standardized amount of fencing}
#'   \item{Fert, Fert_tminus1, ..., Fert_tminus4}{standardized amount of fertilizer used}
#'   \item{Occ, Occ_tminus1, ..., Occ_tminus4}{probability of occupancy}
#'   \item{x, y}{site coordinates, in the British National Grid CRS}
#'   ...
#' }
"example_coord_data"
