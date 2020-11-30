library(survival)
library(caret)
library(randomForestSRC)
library(survcomp)
library(rms)
library(RegParallel)
library(glmnet)

baseline <- function(dataset, n.folds = 10, blocksize = 1000, Cores = 7, KM.Plot = TRUE) {
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
    training_neighbor_times <- data.frame(matrix(nrow = nrow(training_fold), ncol = 3))
    rownames(training_neighbor_times) <- rownames(training_fold)
    colnames(training_neighbor_times)[1] <- "time"
    colnames(training_neighbor_times)[2] <- "status"
    colnames(training_neighbor_times)[3] <- "risk group"

    # Similar dataframe as above, but for the testing data
    testing_neighbor_times <- data.frame(matrix(nrow = nrow(testing_fold), ncol = 3))
    rownames(testing_neighbor_times) <- rownames(testing_fold)
    colnames(testing_neighbor_times)[1] <- "time"
    colnames(testing_neighbor_times)[2] <- "status"
    colnames(testing_neighbor_times)[3] <- "risk group"




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


    training_neighbor_times$time <- training_fold$time
    training_neighbor_times$status <- training_fold$status
    training_neighbor_times$`risk group`<- risk_groups[match(rownames(training_fold), rownames(risk_groups)),]

    testing_neighbor_times$time <- testing_fold$time
    testing_neighbor_times$status <- testing_fold$status
    testing_neighbor_times$`risk group` <- risk_groups[match(rownames(testing_fold), rownames(risk_groups)), ]
    testing_data <- as.data.frame(cbind(testing_neighbor_times, testing_features))


    train <- coxph(Surv(training_fold$time, training_fold$status) ~ training_risk_groups$`Risk Group`, data = training_features)
    test <- coxph(Surv(testing_fold$time, testing_fold$status) ~ testing_neighbor_times$`risk group`, data = testing_features)

    testing_set_concordance <- test[["concordance"]]["concordance"]

    concordance_values[i] <- testing_set_concordance
  }

  if(KM.Plot == TRUE) {
    km.file = readline(prompt = "Enter the full filename (Ex. kmplot.png): ")
    fit <- survfit(Surv(feature_data$time, feature_data$status) ~ 1, data=feature_data)
    t_group <- risk_groups$`Risk Group`
    t_n1 <- sum(t_group == 1)
    t_n2 <- sum(t_group == 2)
    t_n3 <- sum(t_group == 3)
    leg1 <- paste("Medium Risk(" ,t_n2, ")", sep = "")
    leg2 <- paste("High Risk(",t_n1, ")", sep = "")
    leg3 <- paste("Low Risk(" ,t_n3, ")", sep = "")
    cat("-- Writing the Kaplan Meier Plot -- \n")
    png(filename = paste(km.file), width = 7.5, height = 7.5, units = "in", res = 300, pointsize = 7)
    plot(fit, mark.time=TRUE, xlab = "Months", ylab = "Survival", main=paste("Survival Curve",sep=""),lty = 1:3, col = 1:3, cex = 0.5)
    grid()
    legend(x = "bottomright", legend = c(leg1, leg2, leg3), lty = 1:3,
           col = 1:3, cex = 0.65)
    #text(10,0.1,paste("Log-Rank \n p = ", formatC(test_p, format="g", digits = 3), sep = ""),
         #pos = 4, cex = 1)
    dev.off()

  }

  test_concordance <- mean(unlist(concordance_values))
  print(test_concordance)

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

isEmpty <- function(x) {
  return(identical(x, numeric(0)))
}
