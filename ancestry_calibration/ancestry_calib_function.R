ancestry_calib <- function(data, pc_prefix = "pc", num_pcs = 10, score_col = "SCORE") {
    pcs <- paste(pc_prefix, 1:num_pcs, sep = "")

    lm1 <- lm(as.formula(paste0(score_col, " ~ ", paste0(pcs, collapse = " + "))), data = data)
    data$calib1 <- data[[score_col]] - predict(lm1)

    resid_var <- (resid(lm1) - mean(resid(lm1)))^2
    lm2 <- lm(as.formula(paste0("resid_var ~ ", paste0(pcs, collapse = " + "))), data = data)
    predicted_var <- predict(lm2)

    if(any(predicted_var < 0)) {
        print("WARNING: negative variance detected. Linear calibration might not be appropriate for this dataset.")
    }

    data$calib2 <- (data[[score_col]] - predict(lm1)) / sqrt(predicted_var) #if there are negative variances, this will produce NaNs for those individuals

    return(data)

}
