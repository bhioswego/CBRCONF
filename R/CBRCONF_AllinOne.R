#' CBRCONF_AllinOne
#'
#' This function will take a dataset, perform Num runs of cox regression to compute prognostic scores and split into Risk_Group_Num groups.
#' It will then use a case based reasoning method that has a confidence-based retrieval that takes into account each training case's proximity to
#' a prototype constructed from each risk group, and the average distance from that training case to testing cases. It continues to retrieve cases for each
#' testing case until a confidence threshold is reached. In this manner, a different number of cases are retrieved for each testing case.
#' @author Christopher Bartlett
#' @param Dataframe A dataframe where first column is time, second column is status and remaining columns are features.
#' @param Num The number of training sets to create. This is only for determining average beta coefficients in cox regression, reducing features, and determining risk groups.
#' @param Percent Features recurring in this percent of training sets will be kept. Needs to be a decimal value between 0 and 1 (Ex. 50% is 0.5)
#' @param Cores The number of cores to use for parallel processing
#' @param blocksize The number of features to test at once.
#' @param n.folds The number of folds to use for testing
#' @param KM.Plot If this is true, it will save a Kaplan-Meier plot
#' @import stats survival caret survcomp RegParallel
#' @export

library(RegParallel)
library(survival)
library(caret)
library(survcomp)
library(stats)
library(textshape)

CBRCONF_AllinOne <- function(dataset, Num = 10, Percent = 0.5, Cores = 2, blocksize = 4000, n.folds = 10, KM.Plot = TRUE) {
  set.seed(25)
  cat("-- Running the Cox Module -- \n")

  # Dataset is structured where first 2 columns are time and status, and the rest of the columns are features.
  feature_data <- dataset

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

  cat("-- Finding the risk groups -- \n")
  risk_groups <- data.frame(matrix(nrow = nrow(prognostic_scores), ncol = 1))
  rownames(risk_groups) <- rownames(prognostic_scores)
  colnames(risk_groups) <- "Risk Group"
  risk_groups$`Risk Group` <- ile(prognostic_scores$`prognostic score`, 3)


  cat("-- Finding the prototypical representation of each risk group -- \n")
  low_risk_group <- subset(risk_groups, risk_groups$`Risk Group` == 1)
  medium_risk_group <- subset(risk_groups, risk_groups$`Risk Group` == 2)
  high_risk_group <- subset(risk_groups, risk_groups$`Risk Group` == 3)

  low_risk_features <- sig_feature_data[match(rownames(low_risk_group), rownames(sig_feature_data)), ]
  medium_risk_features <- sig_feature_data[match(rownames(medium_risk_group), rownames(sig_feature_data)), ]
  high_risk_features <- sig_feature_data[match(rownames(high_risk_group), rownames(sig_feature_data)), ]

  low_risk_prototype <- as.data.frame(colMeans(low_risk_features))
  colnames(low_risk_prototype) <- "Low Risk Prototype"

  medium_risk_prototype <- as.data.frame(colMeans(medium_risk_features))
  colnames(medium_risk_prototype) <- "Medium Risk Prototype"

  high_risk_prototype <- as.data.frame(colMeans(high_risk_features))
  colnames(high_risk_prototype) <- "High Risk Prototype"

  group_means <- cbind(low_risk_prototype, medium_risk_prototype, high_risk_prototype)

  low_risk_prototype <- as.data.frame(t(low_risk_prototype))
  medium_risk_prototype <- as.data.frame(t(medium_risk_prototype))
  high_risk_prototype <- as.data.frame(t(high_risk_prototype))

  low_risk <- rbind(low_risk_prototype, low_risk_features)
  medium_risk <- rbind(medium_risk_prototype, medium_risk_features)
  high_risk <- rbind(high_risk_prototype, high_risk_features)

  cat("-- Calculating the distance between each sample and the prototype for each risk group -- \n")
  low_risk_distance <- as.matrix(dist(low_risk))
  medium_risk_distance <- as.matrix(dist(medium_risk))
  high_risk_distance <- as.matrix(dist(high_risk))

  distance_from_prototype_low <- as.data.frame(low_risk_distance[,1])
  distance_from_prototype_medium <- as.data.frame(medium_risk_distance[,1])
  distance_from_prototype_high <- as.data.frame(high_risk_distance[,1])

  distance_from_prototype_low <- as.data.frame(normalize(distance_from_prototype_low))
  distance_from_prototype_medium <- as.data.frame(normalize(distance_from_prototype_medium))
  distance_from_prototype_high <- as.data.frame(normalize(distance_from_prototype_high))

  confidence_weights_low <- as.data.frame(1 - distance_from_prototype_low)
  confidence_weights_medium <- as.data.frame(1 - distance_from_prototype_medium)
  confidence_weights_high <- as.data.frame(1 - distance_from_prototype_high)

  colnames(confidence_weights_low) <- "Confidence Weight"
  colnames(confidence_weights_medium) <- "Confidence Weight"
  colnames(confidence_weights_high) <- "Confidence Weight"

  confidence_weights <- rbind(confidence_weights_low, confidence_weights_medium, confidence_weights_high)


  confidence_weights <- confidence_weights * 0.5

  data.folds <- createFolds(y=sig_feature_data$time, k = n.folds)


  concordance_training_list <- rep(NA, n.folds)
  concordance_testing_list <- rep(NA, n.folds)

  number_retrieved <- data.frame(matrix(nrow = nrow(sig_feature_data), ncol = 1))
  colnames(number_retrieved)[1] <- "Number of Retrieved Neighbors"
  rownames(number_retrieved) <- rownames(sig_feature_data)

  testing_sample_information <- list()
  case_retrieval_information <- list()
  training_folds <- list()
  testing_folds <- list()

  cat(" -- Doing ",n.folds," fold testing using Case Based Reasoning with Confidence-Based Retrieval -- \n")
  for(i in 1:n.folds) {
    fold_case_information <- list()
    # Extracting testing and training fold
    testing_fold <- sig_feature_data[data.folds[[i]],]
    training_fold <- sig_feature_data[-data.folds[[i]],]
    training_folds[[i]] <- training_fold
    testing_folds[[i]] <- testing_fold

    testing_features <- testing_fold[,3:ncol(testing_fold)]
    training_features <- training_fold[,3:ncol(training_fold)]

    testing_neighbor_times <- data.frame(matrix(nrow = nrow(testing_fold), ncol = 4))
    rownames(testing_neighbor_times) <- rownames(testing_fold)
    colnames(testing_neighbor_times)[1] <- "time"
    colnames(testing_neighbor_times)[2] <- "status"
    colnames(testing_neighbor_times)[3] <- "prognostic scores"
    colnames(testing_neighbor_times)[4] <- "risk group"

    testing_neighbor_times$status <- testing_fold$status

    testing_distance <- as.matrix(distance_between(testing_features, training_features))

    case_base_average_distances <- as.data.frame(colMeans(testing_distance))
    colnames(case_base_average_distances)[1] <- "Average Distance"
    confidence_temp <- merge(confidence_weights, case_base_average_distances, by = "row.names")
    rownames(confidence_temp) <- confidence_temp[,1]
    confidence_temp <- confidence_temp[,-1]
    confidence_temp$`New Score` <- confidence_temp$`Confidence Weight` + confidence_temp$`Average Distance`
    confidence_weights_average <- as.data.frame(confidence_temp$`New Score`)
    rownames(confidence_weights_average) <- rownames(confidence_temp)
    colnames(confidence_weights_average)[1] <- "Confidence Weight"
    confidence_weights_average <- as.data.frame(normalize(confidence_weights))


    for(k in 1:nrow(testing_features)) {
      query_case <- rownames(testing_features)[k]

      # Creating lists to store survival times and prognostic scores for each query case
      retrieved_neighbors <- list()
      survival_times_for_case <- list()
      prognostic_scores_for_case <- list()
      risk_groups_for_case <- list()
      single_case_information <- list(retrieved_neighbors, survival_times_for_case, prognostic_scores_for_case, risk_groups_for_case)
      names(single_case_information) <- c("retrieved neighbors", "survival times", "prognostic scores", "risk groups")

      ordered_neighbors <- order(testing_distance[k, ])
      name_of_neighbors <- colnames(testing_distance)[ordered_neighbors]
      current_confidence <- 0


      for(n in 1:length(name_of_neighbors)) {
        neighbor <- name_of_neighbors[n]
        confidence <- confidence_weights_average[which(rownames(confidence_weights_average) == neighbor),]
        survival_time <- sig_feature_data[which(rownames(sig_feature_data) == neighbor), ]
        survival_time <- survival_time$time
        single_case_information[[1]] <- append(single_case_information[[1]], neighbor)
        single_case_information[[2]] <- append(single_case_information[[2]], survival_time)
        prognostic_score <- prognostic_scores[which(rownames(prognostic_scores) == neighbor), ]
        single_case_information[[3]] <- append(single_case_information[[3]], prognostic_score)
        risk_group <- risk_groups[which(rownames(risk_groups) == neighbor), ]
        single_case_information[[4]] <- append(single_case_information[[4]], risk_group)
        if(isEmpty(confidence)) {
          confidence <- 0
        }
        current_confidence <- current_confidence + confidence
        if(current_confidence > 1.0 && ! isEmpty(current_confidence)) {
          break
        }
        testing_neighbor_times[which(rownames(testing_neighbor_times) == query_case), 1] <- mean(unlist(single_case_information$`survival times`))
        testing_neighbor_times[which(rownames(testing_neighbor_times) == query_case), 3] <- prognostic_scores[which(rownames(prognostic_scores) == query_case),]

      }
      fold_case_information[[k]] <- single_case_information
      names(fold_case_information)[[k]] <- query_case
      number_retrieved[which(rownames(number_retrieved) == query_case),] <- length(single_case_information[[1]])
    }
    testing_neighbor_times$`risk group` <- risk_groups[match(rownames(testing_fold), rownames(risk_groups)), ]
    testing_statistics <- suppressWarnings(coxph(Surv(testing_neighbor_times$time, testing_neighbor_times$status) ~ testing_neighbor_times$`risk group`, data = testing_features))



    testing_sample_information[[i]] <- testing_neighbor_times

    testing_dataframe <- data.table::rbindlist(lapply(testing_sample_information, setDT, keep.rownames = TRUE))
    testing_dataframe <- column_to_rownames(testing_dataframe, loc = 1)

    testing_set_concordance <- testing_statistics[["concordance"]]["concordance"]
    concordance_testing_list[i] <- testing_set_concordance[1]

    case_retrieval_information[[i]] <- fold_case_information
    names(case_retrieval_information)[[i]] <- paste("fold", i)

  }


  if(KM.Plot == TRUE) {
    km.file = readline(prompt = "Enter the full filename (Ex. kmplot.png): ")
    SurvTimes_Test <- Surv(testing_dataframe$time, testing_dataframe$status)
    log2 = survdiff(SurvTimes_Test ~ testing_dataframe$`risk group`)
    test_p = pchisq(log2$chisq, 1, lower.tail = FALSE)

    fit_test <- survfit(SurvTimes_Test ~ testing_dataframe$`risk group`)
    t_group <- testing_dataframe$`risk group`
    t_n1 <- sum(t_group == 1)
    t_n2 <- sum(t_group == 2)
    t_n3 <- sum(t_group == 3)
    leg1 <- paste("Low Risk(" ,t_n1, ")", sep = "")
    leg2 <- paste("Medium Risk(" ,t_n2, ")", sep = "")
    leg3 <- paste("High Risk(",t_n3, ")", sep = "")
    cat("-- Writing the Kaplan Meier Plot for the test cases to KM-plot.png -- \n")
    png(filename = paste(km.file), width = 7.5, height = 7.5, units = "in", res = 300, pointsize = 7)
    plot(fit_test, mark.time=TRUE, xlab = "Months", ylab = "Survival", main=paste("Survival Curve",sep=""),lty = 1:3, col = 1:3, cex = 0.5)
    grid()
    legend(x = "bottomright", legend = c(leg1, leg2, leg3), lty = 1:3,
           col = 1:3, cex = 0.65)
    text(10,0.1,paste("Log-Rank \n p = ", formatC(test_p, format="g", digits = 3), sep = ""),
         pos = 4, cex = 1)
    dev.off()
  }

  testing_concordance <- mean(concordance_testing_list)
  cat("Testing concordance is: ",testing_concordance,"\n")
}

normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
}

distance_between <- function(df1, df2) {
  distance_frame <- sqrt(matrix(rowSums(expand.grid(rowSums(df1*df1),rowSums(df2*df2))),
                                nrow=nrow(df1)) - (2 * as.matrix(df1) %*% t(as.matrix(df2))))
  return(as.data.frame(distance_frame))
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
