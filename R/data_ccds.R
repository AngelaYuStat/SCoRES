#' Trajectories of pupil response to light after cannabis use
#'
#' Dataset contains functional observation of pupil size percent change after a light stimulus.
#' Participants in the cannabis use group smoked cannabis flower or concentrate 40 minutes prior to the pupillometry measurement.
#' Goal of this data is to understand differences in pupil response to light driven by acute cannabis users.
#' Measurements were collected on the right eye.
#'
#' @format A tibble with N rows and 4 variables:
#' \describe{
#'   \item{subject}{Character. Subject identifier.}
#'   \item{seconds}{Numeric. Time in seconds since light stimulus}
#'   \item{use}{Numeric. Binary indicator of cannabis use 40 minute prior to the light stimulus}
#'   \item{percent_change}{Numeric. Percent change in the outcome of interest.}
#' }
#'
#' @source Processed from \code{data-raw/ccds_load.R} using the \code{readr} and \code{dplyr} packages.
#' @usage data(ccds)
#' @keywords dataset
"ccds"
