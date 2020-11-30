#' Prototype Distance Module
#'
#' This module builds the distance matrices that determine how similar the samples are. Euclidean distance is always used.
#' Confidence weights will start to be built in this module (at this stage, only the distance from the prototype is calculated)
#' @param CBRCONF A CBRCONF object
#' @examples
#' @export

prototype_distance_module <- function(CBRCONF) {
  tic()

  cat(" -- Entering the prototype distance module -- \n")

  low_risk_prototype <- CBRCONF$low_risk_prototype
  medium_risk_prototype <- CBRCONF$medium_risk_prototype
  high_risk_prototype <- CBRCONF$high_risk_prototype

  low_risk_prototype <- as.data.frame(t(low_risk_prototype))
  medium_risk_prototype <- as.data.frame(t(medium_risk_prototype))
  high_risk_prototype <- as.data.frame(t(high_risk_prototype))


  low_risk_features <- CBRCONF$low_risk_features
  medium_risk_features <- CBRCONF$medium_risk_features
  high_risk_features <- CBRCONF$high_risk_features

  low_risk <- rbind(low_risk_prototype, low_risk_features)
  medium_risk <- rbind(medium_risk_prototype, medium_risk_features)
  high_risk <- rbind(high_risk_prototype, high_risk_features)

  cat(" -- Calculating the distances for low, medium and high risk -- \n")
  low_risk_distance <- as.matrix(dist(low_risk))
  medium_risk_distance <- as.matrix(dist(medium_risk))
  high_risk_distance <- as.matrix(dist(high_risk))

  distance_from_prototype_low <- as.data.frame(low_risk_distance[,1])
  distance_from_prototype_medium <- as.data.frame(medium_risk_distance[,1])
  distance_from_prototype_high <- as.data.frame(high_risk_distance[,1])

  distance_from_prototype_low <- as.data.frame(normalize(distance_from_prototype_low))
  distance_from_prototype_medium <- as.data.frame(normalize(distance_from_prototype_medium))
  distance_from_prototype_high <- as.data.frame(normalize(distance_from_prototype_high))

  cat("-- Calculating the confidence weights -- \n")
  confidence_weights_low <- as.data.frame(1 - distance_from_prototype_low)
  confidence_weights_medium <- as.data.frame(1 - distance_from_prototype_medium)
  confidence_weights_high <- as.data.frame(1 - distance_from_prototype_high)

  colnames(confidence_weights_low) <- "Confidence Weight"
  colnames(confidence_weights_medium) <- "Confidence Weight"
  colnames(confidence_weights_high) <- "Confidence Weight"

  confidence_weights <- rbind(confidence_weights_low, confidence_weights_medium, confidence_weights_high)


  confidence_weights <- confidence_weights * 0.5


  cat(" -- Completed -- \n")
  cat(" -- The next module to run is the CBR module -- \n")

  CBRCONF$confidence_weights_low <- confidence_weights_low
  CBRCONF$confidence_weights_medium <- confidence_weights_medium
  CBRCONF$confidence_weights_high <- confidence_weights_high
  CBRCONF$confidence_weights <- confidence_weights

  return(CBRCONF)

}

normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}
