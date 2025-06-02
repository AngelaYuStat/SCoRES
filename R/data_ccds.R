#' Raw longitudinal retinal data from clinical study
#'
#' This dataset contains repeated measures of percent change in some biomarker over time
#' for multiple subjects under different use categories. It includes both user and non-user groups,
#' time points, and metadata related to eye side and frame timing.
#'
#' @format A data frame with N rows and 8 variables:
#' \describe{
#'   \item{subject_id}{Character. Subject ID, e.g., "001-002"}
#'   \item{tp}{Character. Timepoint of measurement, e.g., "post"}
#'   \item{eye}{Character. Eye side measured, "Right" or "Left"}
#'   \item{frame}{Integer. Frame index of the measurement}
#'   \item{percent_change}{Numeric. Percent change}
#'   \item{seconds}{Numeric. Time in seconds}
#'   \item{user_cat}{Factor. Usage category: \code{daily}, \code{non-user}, or \code{occasional}}
#'   \item{use}{Character. Simplified binary usage label: "use" or no use}
#' }
#'
#' @details
#' This dataset is useful for functional regression and longitudinal modeling involving different user behaviors.
#' The `user_cat` column categorizes participants into `"daily"`, `"occasional"`, and `"non-user"` groups based on usage frequency.
#' The `use` column is a simplified binary indicator, where `"use"` corresponds to subjects actively using the treatment/device,
#' typically derived from the `user_cat` value.
#'
#' @usage data(ccds_raw)
#' @keywords dataset
"ccds_raw"

#' Cleaned CCDS dataset (Right Eye, Post Timepoint)
#'
#' This dataset is a cleaned subset of `ccds_raw`, including only measurements taken from the right eye
#' at the post-intervention timepoint (`tp == "post"`). It is preprocessed and ready for functional outcome modeling.
#'
#' @format A tibble with N rows and 4 variables:
#' \describe{
#'   \item{subject_id}{Character. Subject identifier.}
#'   \item{seconds}{Numeric. Time in seconds since initial measurement.}
#'   \item{use}{Numeric. Binary indicator of usage status; \code{1} indicates "use", \code{0} otherwise.}
#'   \item{percent_change}{Numeric. Percent change in the outcome of interest.}
#'   \item{subject}{Character. Equal to subject_id.}
#' }
#'
#' @details
#' The data was filtered to include only right-eye measurements at the post timepoint from the original \code{ccds_raw} dataset.
#' The \code{use} variable is derived from a character column, with \code{"use"} mapped to \code{1}, and all other values to \code{0}.
#'
#' @source Processed from \code{data-raw/ccds_preprocess.R} using the \code{here} and \code{dplyr} packages.
#' @usage data(ccds)
#' @keywords dataset
"ccds"

#' CCDS Data with FPCA Eigenfunctions
#'
#' This dataset augments the cleaned CCDS data (`ccds`) with estimated eigenfunctions obtained from
#' functional principal component analysis (FPCA) on the residuals of a fitted GAM mean model.
#'
#' @format A tibble with N rows and multiple columns:
#' \describe{
#'   \item{subject_id}{Character. Subject identifier.}
#'   \item{seconds}{Numeric. Time in seconds.}
#'   \item{use}{Numeric. Binary usage indicator; 1 = "use", 0 = otherwise.}
#'   \item{percent_change}{Numeric. Observed percent change in outcome.}
#'   \item{subject}{Character. Equal to subject_id.}
#'   \item{Phi1, Phi2, ..., Phi9}{Numeric. Estimated eigenfunctions from FPCA evaluated at each time point.}
#' }
#'
#' @details
#' - A GAM model was first fit on `percent_change` with smooth terms for `seconds` and interaction with `use`.
#' - Residuals were extracted and reshaped into a matrix with subjects as rows and time points as columns.
#' - FPCA was performed on this residual matrix using \code{\link[refund]{fpca.face}}.
#' - The eigenfunctions were merged back into the original CCDS data by `seconds`.
#'
#' This object is useful for functional regression and GAMM+FPCA modeling workflows.
#'
#' @source See script: \code{data-raw/ccds_fpca.R}
#' @usage data(ccds_fpca)
#' @keywords dataset
"ccds_fpca"

#' CMA Output for Functional Regression on CCDS Data
#'
#' This dataset contains the results of contrast-based model assessment (CMA)
#' applied to a functional regression model fit to the CCDS data using subject-level
#' functional principal components.
#'
#' @format A tibble with multiple rows and columns depending on the number of time points:
#' \describe{
#'   \item{time}{Numeric. The time index, corresponding to `seconds`.}
#'   \item{yhat}{Numeric. Estimated difference in fitted values between groups (e.g., use vs non-use).}
#'   \item{scb_low}{Numeric. Lower bound of the simultaneous confidence band for the contrast.}
#'   \item{scb_up}{Numeric. Upper bound of the simultaneous confidence band for the contrast.}
#' }
#'
#' @details
#' - A functional regression model was fit to `percent_change` using subject-level random effects
#'   weighted by the first four eigenfunctions from FPCA.
#' - The fitted model was passed to the `cma()` function.
#' - Simultaneous confidence bands were constructed for the group interested across time.
#'
#' This object is useful for visualizing and testing group-level effects in functional regression models.
#'
#' @source See script: \code{data-raw/ccds_cma.R}
#' @usage data(ccds_cma)
#' @keywords dataset
"ccds_cma"
