#' Cox Module
#' This module performs the feature selection step, in which it runs for Num iterations. At each iteration, a training set
#' and testing set will be created and every feature will run through a multivariate cox regression. After Num runs, features found to be
#' significant in the percent of training sets specified with Percent are kept (log ranking < 0.01 is used for significance). Cores
#' and blocksize are used for a parallel process (uses RegParallel). Cores is the number of computing cores assigned to the task, and
#' blocksize is the number of features to run at once. Blocksize should not exceed the number of features in your dataset.
#'
#' @param CBRCONF A CBRCONF object
#' @param Num The number of iterations to run
#' @param Percent Features recurring in this percent of training sets will be kept. Needs to be a decimal value between 0 and 1 (Ex. 50% is 0.5)
#' @param Cores The number of cores to use for parallel processing
#' @import stats survival caret survcomp RegParallel
#' @export

cox_module <- function(CBRCONF, Num, Percent, Cores) {
  library(tictoc)
  set.seed(25)
  tic()
  cat("-- Entering the Cox Module -- \n")

  # Dataset is structured where first 2 columns are time and status, and the rest of the columns are features.
  feature_data <- CBRCONF$dataset

  training_sets <- list()
  testing_sets <- list()

  sig_feature_list <- list()
  results_list <- list()

  # This removes the time and status and leaves us with just the features.
  features <- feature_data[,3:ncol(feature_data)]

  retain_threshold <- Percent * Num
  feature_count <- data.frame(matrix(nrow = ncol(features), ncol = 1))
  rownames(feature_count) <- colnames(features)
  colnames(feature_count)[1] <- "Count"

  blocksize = ncol(feature_data) / 4

  for(i in 1:Num) {
    cat("-- Running a multivariate cox regression on training set ",i," -- \n")
    # Extracting testing and training fold
    sample <- sample.int(n = nrow(feature_data), size = floor(.75*nrow(feature_data)), replace = F)
    training_fold <- feature_data[sample, ]
    # Testing folds are presently not used during feature selection.
    testing_fold  <- feature_data[-sample, ]



    # Running a multivariate cox regression to test each feature for significance
    res <- RegParallel(
      data = training_fold,
      formula = 'Surv(time, status) ~ [*]',
      FUN = function(formula, data)
        coxph(formula = formula,
              data = data,
              ties = 'efron',
              singular.ok = TRUE),
      FUNtype = 'coxph',
      variables = colnames(training_fold)[3:ncol(training_fold)],
      blocksize = blocksize,
      cores = Cores,
      nestedParallel = FALSE,
      conflevel = 95)

    res <- res[!is.na(res$P),]
    res <- res[order(res$LogRank, decreasing = FALSE),]
    results_list[[i]] <- res
    # Getting only the features with a log rank less than 0.01
    final <- subset(res, LogRank < 0.01)
    sig_features <- gsub('^X', '',final$Variable)
    sig_feature_list[[i]] <- sig_features

  }

  # Initiating a dataframe to contain the features we want to keep
  features_to_keep <- rep(NA, nrow(feature_count))

  cat("-- Getting the counts of how many times each feature shows up as significant --\n")
  for(f in 1:nrow(feature_count)) {
    query <- rownames(feature_count)[f]
    current_count <- 0
    for(i in 1:Num) {
      if(query %in% sig_feature_list[[i]]) {
        current_count <- current_count + 1
      }
    }
    feature_count[f,] <- current_count
    if(feature_count[f,] > retain_threshold) {
      features_to_keep[f] <- rownames(feature_count)[f]
    }
  }

  # What remains in this list are the features that occurred in the specified percentage of training sets
  features_to_keep <- na.omit(features_to_keep)

  beta_results <- list()
  # Now we want just the Beta coefficient results
  for(r in 1:length(results_list)){
    res <- results_list[[r]]
    beta_table <- data.frame(matrix(nrow = ncol(features), ncol = 1))
    rownames(beta_table) <- res$Variable
    colnames(beta_table)[1] <- "Beta"
    beta_table$Beta <- res$Beta
    beta_results[[r]] <- beta_table
  }

  mean_beta_results <- list()
  # Now we want to determine the average Beta coefficient for each feature.
  for(b in 1:length(beta_results)) {
    beta_table1 <- beta_results[[b]]
    beta_table1 <- as.data.frame(beta_table1[match(features_to_keep, rownames(beta_table1)), ])
    rownames(beta_table1) <- features_to_keep
    colnames(beta_table1) <- "Beta"
    mean_beta_results[[b]] <- beta_table1
  }

  cat(" -- Calculating the mean beta coefficient of each feature -- \n")
  Beta_Means <- as.data.frame(rowMeans(do.call(cbind, lapply(mean_beta_results, "[", "Beta"))))
  colnames(Beta_Means) <- "Means"
  # Reduce the feature data to just those features that were significant
  sig_feature_data <- features[,match(features_to_keep, colnames(features))]

  # Setting up a table to hold the samples and prognostic scores
  prognostic_scores <- data.frame(matrix(nrow = nrow(sig_feature_data), ncol = 1))
  rownames(prognostic_scores) <- rownames(sig_feature_data)
  colnames(prognostic_scores) <- "prognostic score"

  # Setting up a dataframe to contain all of the individual feature scores. This is the feature's expression value multiplied by the mean beta coefficient
  feature_scores <- data.frame(matrix(nrow = nrow(sig_feature_data), ncol = ncol(sig_feature_data)))
  rownames(feature_scores) <- rownames(sig_feature_data)
  colnames(feature_scores) <- colnames(sig_feature_data)

  # Iterating through the significant feature data, and multiplying by the mean beta coefficient
  for(row in 1:nrow(sig_feature_data)) {
    for(column in 1:ncol(sig_feature_data)) {
      feature_score <- sig_feature_data[row, column] * Beta_Means[column,]
      feature_scores[row, column] <- feature_score
    }
  }

  cat(" -- Calculating the prognostic score for each sample -- \n")
  for(sample in 1:nrow(prognostic_scores)) {
    prognostic_scores[sample,] <- rowSums(feature_scores[sample,])
  }

  sig_feature_data <- cbind(feature_data$status, feature_data$time, sig_feature_data)
  colnames(sig_feature_data)[1] <- "status"
  colnames(sig_feature_data)[2] <- "time"

  cat(" -- Completed -- \n")
  cat(" -- The next module to run is the Risk Groups Module -- \n")

  toc()
  CBRCONF$sig_feature_data <- sig_feature_data
  CBRCONF$sig_feature_list <- sig_feature_list
  CBRCONF$feature_scores <- feature_scores
  CBRCONF$prognostic_scores <- prognostic_scores
  CBRCONF$results_list <- results_list

  return(CBRCONF)
}
