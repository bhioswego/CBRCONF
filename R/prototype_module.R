#' Prototype Module
#'
#' This module creates 3 prototypical representations of the risk groups. This is performed by using all samples within a risk
#' group and taking an average of each feature. Specifically, a vector is constructed for low, medium and high risk that has
#' averaged values for each feature.
#' @param CBRCONF A CBRCONF object
#' @examples
#' @export


prototype_module <- function(CBRCONF) {
  tic()
  cat("-- Entering the Prototype Module -- \n")
  risk_groups <- CBRCONF$risk_groups

  sig_feature_data <- CBRCONF$sig_feature_data

  low_risk_group <- subset(risk_groups, risk_groups$`Risk Group` == 1)
  medium_risk_group <- subset(risk_groups, risk_groups$`Risk Group` == 2)
  high_risk_group <- subset(risk_groups, risk_groups$`Risk Group` == 3)

  low_risk_features <- sig_feature_data[match(rownames(low_risk_group), rownames(sig_feature_data)), ]
  medium_risk_features <- sig_feature_data[match(rownames(medium_risk_group), rownames(sig_feature_data)), ]
  high_risk_features <- sig_feature_data[match(rownames(high_risk_group), rownames(sig_feature_data)), ]

  cat("-- Making the low risk prototype -- \n")
  low_risk_prototype <- as.data.frame(colMeans(low_risk_features))
  colnames(low_risk_prototype) <- "Low Risk Prototype"

  cat("-- Making the medium risk prototype -- \n")
  medium_risk_prototype <- as.data.frame(colMeans(medium_risk_features))
  colnames(medium_risk_prototype) <- "Medium Risk Prototype"

  cat("-- Making the high risk prototype -- \n")
  high_risk_prototype <- as.data.frame(colMeans(high_risk_features))
  colnames(high_risk_prototype) <- "High Risk Prototype"

  group_means <- cbind(low_risk_prototype, medium_risk_prototype, high_risk_prototype)


  cat(" -- Completed -- \n")
  cat(" -- The next module to run is the Prototype Distance Module -- \n")
  CBRCONF$low_risk_features <- low_risk_features
  CBRCONF$medium_risk_features <- medium_risk_features
  CBRCONF$high_risk_features <- high_risk_features


  CBRCONF$low_risk_prototype <- low_risk_prototype
  CBRCONF$medium_risk_prototype <- medium_risk_prototype
  CBRCONF$high_risk_prototype <- high_risk_prototype

  CBRCONF$group_means <- group_means

  return(CBRCONF)

}
