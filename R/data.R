#' Simulated Deep Vein Thrombosis (DVT) Dataset
#'
#' A simulated dataset of 9,805 patients from seven studies, used to illustrate
#' the development of a clinical prediction model for diagnosing deep vein
#' thrombosis (DVT). All predictors are binary (0/1).
#'
#' @format A data frame with 9,805 rows and 12 columns:
#' \describe{
#'   \item{male}{Male sex (1 = male)}
#'   \item{OCuse}{Oral contraceptive use (1 = yes)}
#'   \item{malignancy}{Active malignancy/cancer (1 = yes)}
#'   \item{surgery}{Recent surgery (1 = yes)}
#'   \item{ab_leg_trauma}{Recent leg trauma (1 = yes)}
#'   \item{vein_distension}{Leg vein distension/swelling (1 = yes)}
#'   \item{calf_diff}{Calf width difference >3 cm (1 = yes)}
#'   \item{ab_ddimer}{Abnormal D-dimer blood test result (1 = yes)}
#'   \item{pregnancy}{Recent pregnancy (1 = yes)}
#'   \item{immobile}{Recent immobility >24 hours (1 = yes)}
#'   \item{clotting}{Recent anticoagulant use (1 = yes)}
#'   \item{dvt}{DVT diagnosis (1 = DVT present; outcome variable)}
#' }
#' @details
#' Outcome prevalence is approximately 21%. The dataset spans seven studies of
#' varying size. It is intended as a realistic example of a modest-sized
#' development dataset for use with the generalised sample size approach
#' ([pmgen_reference()], [pmgen_simulate()]).
#'
#' @source Simulated dataset provided for illustrative purposes.
"dvt_data"
