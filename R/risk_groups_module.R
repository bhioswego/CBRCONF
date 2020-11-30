#' Risk Groups Module
#'
#' This module constructs 3 prognostic risk groups, taken from the ile function. In other words, each risk group is derived from
#' a tertile of the prognostic score. A KM plot can be constructed at this point from the actual survival data.
#'
#' @param CBRCONF A CBRCONF object
#' @param KM.Plot If KM.Plot is true, a KM plot of 3 risk groups (low, medium, high) will be created and stored in the current work directory
#' @examples
#' @export

risk_groups_module <- function(CBRCONF, KM.Plot = TRUE) {
  library(tictoc)
  cat("-- Entering the Risk Groups Module -- \n")
  tic()
  feature_data <- CBRCONF$dataset
  prognostic_scores <- CBRCONF$prognostic_scores
  risk_groups <- data.frame(matrix(nrow = nrow(prognostic_scores), ncol = 1))
  rownames(risk_groups) <- rownames(prognostic_scores)
  colnames(risk_groups) <- "Risk Group"

  cat("-- Forming the risk groups -- \n")
  risk_groups$`Risk Group` <- ile(prognostic_scores$`prognostic score`, 3)

  if(KM.Plot == TRUE) {
    SurvTimes <- Surv(feature_data$time, feature_data$status)
    log1 = survdiff(SurvTimes ~ risk_groups$`Risk Group`)
    p = pchisq(log1$chisq, 1, lower.tail=FALSE)


    fit <- survfit(SurvTimes ~ risk_groups$`Risk Group`)
    group <- risk_groups$`Risk Group`
    n1 <- sum(group == 1)
    n2 <- sum(group == 2)
    n3 <- sum(group == 3)
    leg1 <- paste("Low Risk(" ,n1, ")", sep = "")
    leg2 <- paste("Medium Risk(" ,n2, ")", sep = "")
    leg3 <- paste("High Risk(", n3, ")", sep = "")
    cat("-- Writing the Kaplan Meier Plot to KM-plot.png -- ")
    png(filename = "KM-Plot.png", width = 5.5, height = 5.5, units = "in", res = 300, pointsize = 7)
    plot(fit, mark.time=TRUE, xlab = "Months", ylab = "Survival", main=paste("Survival Curve",sep=""),lty = 1:3, col = 1:3, cex = 0.5)
    grid()
    legend(x = "bottomright", legend = c(leg1, leg2, leg3), lty = 1:3,
           col = 1:3, cex = 0.65)
    text(10,0.1,paste("Log-Rank \n p = ", formatC(p, format="g", digits = 3), sep = ""),
         pos = 4, cex = 1)
    dev.off()
  }

  cat(" -- Completed -- \n")
  cat(" -- The next module to run is the Prototype Module -- \n")
  toc()
  CBRCONF$risk_groups <- risk_groups
  return(CBRCONF)

}

ile <- function (x, levels.out) {
  x <- as.double(x)
  ## if no levels.out specified, use three levels named 1,2,3
  if (missing(levels.out)) levels.out <- 1:3
  ## if levels.out is an integer of length 1, assume numbered levels
  if (identical(length(levels.out), 1L)) levels.out <- 1:as.integer(levels.out)
  n <- length(levels.out)
  xd <- ceiling(n * (rank(x, na.last = "keep") - 0.5)/sum(!is.na(x)))
  return(factor(levels.out[xd], levels.out))
}

