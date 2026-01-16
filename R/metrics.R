#' Evaluate Selection Metrics
#'
#' Computes various classification metrics for variable selection.
#'
#' @param estimate Vector of estimated coefficients.
#' @param truth Vector of true coefficients.
#' @param threshold Threshold for considering a coefficient as "selected" (non-zero).
#'
#' @return A list containing:
#' \itemize{
#'   \item TP: True Positives
#'   \item FP: False Positives
#'   \item TN: True Negatives
#'   \item FN: False Negatives
#'   \item TPR: True Positive Rate (Sensitivity/Recall)
#'   \item TNR: True Negative Rate (Specificity)
#'   \item FPR: False Positive Rate
#'   \item FDR: False Discovery Rate
#'   \item Precision: Positive Predictive Value
#'   \item Accuracy: Overall Accuracy
#'   \item F1: F1 Score
#'   \item MCC: Matthews Correlation Coefficient
#' }
#' @export
evaluate_selection_metrics <- function(estimate, truth, threshold = 1e-5) {
    if (length(estimate) != length(truth)) {
        cli::cli_abort("Length of estimate ({length(estimate)}) must match length of truth ({length(truth)}).")
    }

    is_selected <- abs(estimate) > threshold
    is_true <- abs(truth) > threshold

    # Confusion Matrix
    TP <- sum(is_selected & is_true)
    FP <- sum(is_selected & !is_true)
    TN <- sum(!is_selected & !is_true)
    FN <- sum(!is_selected & is_true)

    # Derived Metrics
    # Avoid division by zero
    safe_div <- function(n, d) ifelse(d == 0, 0, n / d)

    TPR <- safe_div(TP, TP + FN) # Sensitivity / Recall
    TNR <- safe_div(TN, TN + FP) # Specificity
    FPR <- safe_div(FP, FP + TN) # 1 - Specificity
    FDR <- safe_div(FP, TP + FP) # False Discovery Rate
    Precision <- safe_div(TP, TP + FP)
    Accuracy <- safe_div(TP + TN, TP + TN + FP + FN)

    F1 <- safe_div(2 * Precision * TPR, Precision + TPR)

    # MCC
    mcc_num <- (as.numeric(TP) * TN) - (as.numeric(FP) * FN)
    mcc_den <- sqrt(TP + FP) * sqrt(TP + FN) * sqrt(TN + FP) * sqrt(TN + FN)
    MCC <- safe_div(mcc_num, mcc_den)
    # Clamp MCC to range [-1, 1] if valid, though safe_div prevents NaN
    if (!is.na(MCC)) MCC <- max(-1, min(1, MCC))

    return(list(
        TP = TP, FP = FP, TN = TN, FN = FN,
        TPR = TPR, TNR = TNR, FPR = FPR, FDR = FDR,
        Precision = Precision, Accuracy = Accuracy,
        F1 = F1, MCC = MCC
    ))
}
