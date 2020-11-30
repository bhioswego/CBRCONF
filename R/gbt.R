library(survival)
library(caret)
library(gbm)
library(survcomp)
library(rms)
library(RegParallel)
library(ranger)
library(randomForestSRC)


#' Gradient Boosting Tree Module
#'
#' This module will use gradient boosting tree to predict survival.
#'
#'
#' @param dataset The dataset to use
#' @param n.folds The number of folds
#' @param ntree The number of trees in the random survival forest
#' @param blocksize The number of features to run through cox regression at once.
#' @import caret
#' @examples
#' @export
#'

gbt <- function(dataset, n.folds = 10, n.tree = 500, blocksize = 1000, Cores = 7) {
  feature_data <- dataset
  concordance_values <- list()
  #risk_groups <- CBRContainer$risk_groups


  concordance_training_list <- rep(NA, n.folds)
  concordance_testing_list <- rep(NA, n.folds)

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
    cores = Cores,
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
    testing_features <- testing_fold[,3:ncol(testing_fold)]


    cat(" -- Calculating the risk groups from survival times for fold ",i," -- \n")
    training_risk_groups <- data.frame(matrix(nrow = nrow(training_features), ncol = 1))
    rownames(training_risk_groups) <- rownames(training_features)
    colnames(training_risk_groups) <- "riskgroup"
    training_risk_groups$`Risk Group` <- risk_groups[match(rownames(training_fold), rownames(risk_groups)), ]
    training_fold <- cbind(training_risk_groups$`Risk Group`, training_fold)
    colnames(training_fold)[1] <- "riskgroup"

    testing_risk_groups <- data.frame(matrix(nrow = nrow(testing_features), ncol = 1))
    rownames(testing_risk_groups) <- rownames(testing_features)
    colnames(testing_risk_groups) <- "riskgroup"
    testing_risk_groups$`Risk Group` <- risk_groups[match(rownames(testing_fold), rownames(risk_groups)), ]
    testing_fold <- cbind(testing_risk_groups$`Risk Group`, testing_fold)
    colnames(testing_fold)[1] <- "riskgroup"


    cat(" -- Running a cox regression for fold ",i," -- \n")

    train <- gbm(Surv(time, status)~.,
               training_fold,
               interaction.depth=2,
               shrinkage=0.01,
               n.trees=n.tree,
               distribution="coxph")

    testing_model <- predict.gbm(train,
                            n.trees=n.tree,
                            newdata=testing_fold,
                            type="response")

    testing_concordance <- concordance.index(testing_model, testing_fold$time, testing_fold$status)


    concordance_values[i] <- testing_concordance$c.index

  }

  test_concordance <- mean(unlist(concordance_values))
  print(test_concordance)

}
