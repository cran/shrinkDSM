#' Survival times of gastric cancer patients
#'
#' A data set of survival times of patients with locally advanced,
#' nonresectable gastric carcinoma. The patients were either
#' treated with chemotherapy plus radiation or chemotherapy alone.
#'
#' @format A data frame with 90 rows and 4 variables:
#' \describe{
#'   \item{id}{patient id}
#'   \item{radiation}{dummy variable indicating which treatment was employed,
#'   0 = chemotherapy, 1 = combined chemotherapy/radiation}
#'   \item{time}{time survived by patient in days}
#'   \item{status}{dummy variable indicating whether death of the patient was observed,
#'   0 = death not observed (i.e. censored), 1 = death observed.}
#' }
#' @source Moreau, T., O'Quigley, J., and Mesbah M. (1985) A global goodness-of-fit statistic for
#' the proportional hazards model Appl. Statist., 34, 212:218 (p 213)
#'
#' \url{https://www.mayo.edu/research/documents/gastrichtml/DOC-10027680}
"gastric"
