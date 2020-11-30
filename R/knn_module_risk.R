#' KNN Module
#'
#' This module will use K-nearest neighbor to predict survival.
#'
#'
#' @param CBRContainer A CBRContainer object
#' @param K The value of K for K Nearest Neighbors
#' @param n.folds The number of folds
#' @param KM.Plot If KM.Plot is true, a kaplan-meier plot will be constructed from 3 risk groups (low, medium and high)
#' @examples
#' @export

knn_module_risk <- function(CBRContainer, k = 5, n.folds = 10, blocksize = 100, cores = 7) {
  feature_data <- CBRContainer$dataset

  set.seed(25)
  cat("-- Setting up the training and testing folds -- \n")
  data.folds <- createFolds(y=feature_data$time, k = n.folds)

  features <- feature_data[,3:ncol(feature_data)]

  # Setting up lists for the concordance index values
  concordance_training_list <- rep(NA, n.folds)
  concordance_testing_list <- rep(NA, n.folds)

  prognostic_scores <- data.frame(matrix(nrow = nrow(feature_data), ncol = 1))
  rownames(prognostic_scores) <- rownames(feature_data)
  colnames(prognostic_scores) <- "prognostic score"

  # Setting up a dataframe to contain all of the individual feature scores. This is the feature's expression value multiplied by the mean beta coefficient
  feature_scores <- data.frame(matrix(nrow = nrow(features), ncol = ncol(features)))
  rownames(feature_scores) <- rownames(features)
  colnames(feature_scores) <- colnames(features)

  cat("-- Determining beta coefficients for features so we can find risk groups -- \n")
  res <- RegParallel(
    data = feature_data,
    formula = 'Surv(time, status) ~ [*]',
    FUN = function(formula, data)
      coxph(formula = formula,
            data = data,
            ties = 'efron',
            singular.ok = TRUE),
    FUNtype = 'coxph',
    variables = colnames(feature_data)[3:ncol(feature_data)],
    blocksize = blocksize,
    cores = cores,
    nestedParallel = FALSE,
    conflevel = 95)

  res <- res[!is.na(res$P),]
  res <- res[order(res$LogRank, decreasing = FALSE),]

  cat("-- Calculating feature scores for all samples -- \n")
  for(row in 1:nrow(features)) {
    for(column in 1:ncol(features)) {
      feature_score <- features[row, column] * res$Beta[column]
      feature_scores[row, column] <- feature_score
    }
  }

  # Calculate the prognostic score for each sample
  for(sample in 1:nrow(prognostic_scores)) {
    prognostic_scores[sample,] <- rowSums(feature_scores[sample,])
  }

  risk_groups <- data.frame(matrix(nrow = nrow(prognostic_scores), ncol = 1))
  rownames(risk_groups) <- rownames(prognostic_scores)
  colnames(risk_groups) <- "Risk Group"

  risk_groups$`Risk Group` <- ile(prognostic_scores$`prognostic score`, 3)


  for(i in 1:n.folds) {
    # Extracting testing and training fold
    testing_fold <- feature_data[data.folds[[i]],]
    training_fold <- feature_data[-data.folds[[i]],]



    # Setting up a dataframe to hold the mean survival times and survival statuses for training
    training_neighbor_times <- data.frame(matrix(nrow = nrow(training_fold), ncol = 2))
    rownames(training_neighbor_times) <- rownames(training_fold)
    colnames(training_neighbor_times)[1] <- "time"
    colnames(training_neighbor_times)[2] <- "status"

    # Similar dataframe as above, but for the testing data
    testing_neighbor_times <- data.frame(matrix(nrow = nrow(testing_fold), ncol = 2))
    rownames(testing_neighbor_times) <- rownames(testing_fold)
    colnames(testing_neighbor_times)[1] <- "time"
    colnames(testing_neighbor_times)[2] <- "status"




    # TRAINING
    # -------------------------------------------------

    # Getting only the features and removing time and status
    training_features <- training_fold[,3:ncol(training_fold)]


    # Building a distance matrix
    training_dist <- as.matrix(dist(training_features))

    cat("-- Running KNN on training fold ",i," with ",k," neighbors -- \n")
    for (j in 1:nrow(training_fold)) {
      indices <- knn_for_one(j, training_dist, k = k)
      neighbor_pheno <- training_fold[indices,c(1,2)]
      training_neighbor_times[j,] <- mean(neighbor_pheno$time)
      training_neighbor_times$status <- training_fold$status
    }


    cat(" -- Calculating the risk groups from survival times for fold ",i," -- \n")
    training_risk_groups <- data.frame(matrix(nrow = nrow(training_features), ncol = 1))
    rownames(training_risk_groups) <- rownames(training_features)
    colnames(training_risk_groups) <- "Risk Group"
    training_risk_groups$`Risk Group` <- risk_groups[match(rownames(training_fold), rownames(risk_groups)), ]

    cat(" -- Running a cox regression for fold ",i," -- \n")
    training_statistics <- coxph(Surv(training_neighbor_times$time, training_neighbor_times$status) ~ training_risk_groups$`Risk Group`, data = training_features)

    # Getting the c-index and adding it to the list so we can calculate the mean c-index later
    training_set_concordance <- training_statistics[["concordance"]]["concordance"]
    concordance_training_list[i] <- training_set_concordance[1]

    # Testing
    # -------------------------------------------------

    testing_features <- testing_fold[,3:ncol(testing_fold)]


    # Building a distance matrix
    testing_dist <- as.matrix(dist(testing_features))

    cat("-- Running KNN on testing fold ",i," with ",k," neighbors -- \n")
    for (j in 1:nrow(testing_fold)) {
      indices <- knn_for_one(j, testing_dist, k = k)
      # The following line worked during testing, but I think it's grabbing the wrong data
      # neighbor_pheno <- testing_fold[indices, c(1,2)]
      neighbor_pheno <- feature_data[indices,c(1,2)]
      testing_neighbor_times[j,] <- mean(neighbor_pheno$time)
      testing_neighbor_times$status <- testing_fold$status
    }


    cat(" -- Calculating the risk groups from survival times for fold ",i," -- \n")
    testing_risk_groups <- data.frame(matrix(nrow = nrow(testing_features), ncol = 1))
    rownames(testing_risk_groups) <- rownames(testing_features)
    colnames(testing_risk_groups) <- "Risk Group"
    testing_risk_groups$`Risk Group` <- risk_groups[match(rownames(testing_fold), rownames(risk_groups)), ]


    cat(" -- Running a cox regression for fold ",i," -- \n")
    testing_statistics <- coxph(Surv(testing_neighbor_times$time, testing_neighbor_times$status) ~ testing_risk_groups$`Risk Group`, data = testing_features)

    # Getting the c-index and adding it to the list so we can calculate the mean c-index later
    testing_set_concordance <- testing_statistics[["concordance"]]["concordance"]
    concordance_testing_list[i] <- testing_set_concordance[1]

    # if(KM.Plot == TRUE) {
    #   cat("-- Writing KM Plot for Fold ",i," -- \n")
    #   SurvTimes <- Surv(testing_neighbor_times$time, testing_neighbor_times$status)
    #   log1 = survdiff(SurvTimes ~ testing_risk_groups$`Risk Group`)
    #   p = pchisq(log1$chisq, 1, lower.tail=FALSE)
    #
    #
    #   fit <- survfit(SurvTimes ~ testing_risk_groups$`Risk Group`)
    #   group <- testing_risk_groups$`Risk Group`
    #   n1 <- sum(group == 1)
    #   n2 <- sum(group == 2)
    #   n3 <- sum(group == 3)
    #   leg1 <- paste("Low Risk(" ,n1, ")", sep = "")
    #   leg2 <- paste("Medium Risk(" ,n2, ")", sep = "")
    #   leg3 <- paste("High Risk(", n3, ")", sep = "")
    #   cat("-- Writing the Kaplan Meier Plot to KM-plot.png -- \n")
    #   png(filename = paste("KM-Plot_Fold",i,".png"), width = 5.5, height = 5.5, units = "in", res = 300, pointsize = 7)
    #   plot(fit, mark.time=TRUE, xlab = "Months", ylab = "Survival", main=paste("Survival Curve",sep=""),lty = 1:3, col = 1:3, cex = 0.5)
    #   grid()
    #   legend(x = "bottomright", legend = c(leg1, leg2, leg3), lty = 1:3,
    #          col = 1:3, cex = 0.65)
    #   text(10,0.1,paste("Log-Rank \n p = ", formatC(p, format="g", digits = 3), sep = ""),
    #        pos = 4, cex = 1)
    #   dev.off()
    # }
  }
  # Calculate average c-index for training and report it
  training_concordance <- mean(concordance_training_list)
  cat("Training concordance is: ",training_concordance,"\n")

  # Calculate average c-index for testing and report it
  testing_concordance <- mean(concordance_testing_list)
  cat("Testing concordance is: ",testing_concordance,"\n")
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

knn_for_one <- function(i, distance_matrix, k) {
  ordered_neighbors <<- order(distance_matrix[i, ])
  return(ordered_neighbors[2:(k + 1)])
}
