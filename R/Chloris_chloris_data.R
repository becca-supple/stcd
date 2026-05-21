#' Example spatio-temporal data for the European greenfinch (Chloris chloris)
#'
#' Time series of greenfinch abundance and environmental, land use, and agricultural variables from 1992 to 2016 across Europe.
#'
#' @format ## `Chloris_chloris_data`
#' A data frame with 7,240 rows and 60 columns:
#' \describe{
#'   \item{CountryGroup}{Country name}
#'   \item{Year}{Year}
#'   \item{Abundance}{Estimated/imputed number of individuals during breeding season.}
#'   \item{Agr}{Percent surface area used for agriculture.}
#'   \item{Art}{Percent urban or artificial surface area.}
#'   \item{Wood}{Percent surface area covered by woodlands.}
#'   \item{Fert}{Total nitrogen-based fertilizer used in kg.}
#'   \item{Pest_tot}{Total pesticide use in tonnes.}
#'   \item{Temp}{Mean temperature in degrees Celcius}
#' }
#'
#' @source Modified from Rigal, S. (2023). R scripts and data for the following article: ”Farmland practices are driving bird populations decline across Europe.”.
#' <https://github.com/StanislasRigal/Drivers_European_bird_decline_public/blob/main/README.md.>
#'
"Chloris_chloris_data"
