
# use in summarise
# returns tibble
binary_measures <- function(real, predicted, na.rm = TRUE) {

  require(tidyverse)

  assertthat::assert_that(
    rlang::is_logical(real),
    rlang::is_logical(predicted),
    na.rm || !(any(is.na(real)) | any(is.na(predicted))))

  if (na.rm) {
    na_at <- is.na(real) | is.na(predicted)
    real <- real[!na_at]
    predicted <- predicted[!na_at]
  }

  TP <- as.numeric(sum(real & predicted))
  TN <- as.numeric(sum((!real) & (!predicted)))
  FP <- as.numeric(sum((!real) & predicted))
  FN <- as.numeric(sum(real & (!predicted)))

  assertthat::assert_that(sum(TP, TN, FP, FN) == (length(real)))

  # see https://doi.org/10.1186/s12864-019-6413-7 for limits on MCC and f1_score when rows/columns in confusion matrix sum to zero
  # limits for precision, recall and specificity equal to 0.5, which is perfomance of random classifier
  res <-
    tibble(
      accuracy = (TP + TN) / (TP + FP + TN + FN),
      precision = TP / (TP + FP),
      recall = TP / (TP + FN),
      specificity = TN / (TN + FP),
      MCC = ((TP * TN) - (FP * FN)) / (sqrt((TP+FN)*(TP+FP)*(TN+FP)*(TN+FN)))) %>%
    mutate(
      MCC         = if_else(is.nan(MCC),           0, MCC),        # limit
      precision   = if_else(is.nan(precision),   0.5, precision),  # limit
      recall      = if_else(is.nan(recall),      0.5, recall),     # limit
      specificity = if_else(is.nan(specificity), 0.5, specificity),# limit
      f1_score = case_when(
        (TP == 0 & FP == 0 & FN == 0) ~ 1,
        (TP == 0 & FP >  0 & FN >  0) ~ 0,
        TRUE ~ (2 * precision * recall) / (precision + recall)))

  return(res)
}

RMSE <- function(a, b, real=NULL, predicted=NULL) {
  resid <- a-b
  if (!is.null(real) & !is.null(predicted)) {
    resid <- resid[real & predicted]
  }
  sqrt(mean((resid)**2))
}

SMAPE <- function(a, b, real=NULL, predicted=NULL) {
  resid <- a-b
  if (!is.null(real) & !is.null(predicted)) {
    resid <- resid[real & predicted]
  }
  mean(abs(a-b) / (abs(a) + abs(b)))
}

# absolute average fold error
AAFE <- function(a, b,  real=NULL, predicted=NULL) {
  if (!is.null(real) & !is.null(predicted)) {
    a <- a[real & predicted]
    b <- b[real & predicted]
  }
  10^(mean(abs(log(b/a))))
}

MAE <- function(a, b, real=NULL, predicted=NULL, na.rm = TRUE) {
  resid <- a-b
  if (!is.null(real) & !is.null(predicted)) {
    resid <- resid[real & predicted]
  }
  mean(abs(resid), na.rm = na.rm)
}

MedAE <- function(a, b,real=NULL, predicted=NULL,  na.rm = TRUE) {
  resid <- a-b
  if (!is.null(real) & !is.null(predicted)) {
    resid <- resid[real & predicted]
  }
  median(abs(resid), na.rm = na.rm)
}

