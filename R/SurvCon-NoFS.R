#' SurvCon-CBR without Feature Selection Module
#'
#' This module will use SurvCon-CBR to predict survival.
#'
#'
#' @param dataset The dataset must have time and status as the first two columns, followed by the features.
#' @param n.folds The number of folds
#' @param KM.Plot If KM.Plot is true, a kaplan-meier plot will be constructed from 3 risk groups (low, medium and high)
#' @import caret dplyr data.table textshape
#' @examples
#' @export

library(RegParallel)
library(survival)
library(caret)
library(survcomp)
library(stats)
library(textshape)

SurvCon_NoFS <- function(dataset, n.folds = 10, blocksize = 100, Cores = 1, KM.Plot = TRUE) {
  feature_data <- dataset
  concordance_values <- list()
  #risk_groups <- CBRContainer$risk_groups


  concordance_training_list <- rep(NA, n.folds)
  concordance_testing_list <- rep(NA, n.folds)

  set.seed(25)
  cat("-- Setting up the training and testing folds -- \n")

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

  low_risk_group <- subset(risk_groups, risk_groups$`Risk Group` == 1)
  medium_risk_group <- subset(risk_groups, risk_groups$`Risk Group` == 2)
  high_risk_group <- subset(risk_groups, risk_groups$`Risk Group` == 3)

  low_risk_features <- feature_data[match(rownames(low_risk_group), rownames(feature_data)), ]
  medium_risk_features <- feature_data[match(rownames(medium_risk_group), rownames(feature_data)), ]
  high_risk_features <- feature_data[match(rownames(high_risk_group), rownames(feature_data)), ]

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
  data.folds <- createFolds(y=feature_data$time, k = n.folds)


  concordance_training_list <- rep(NA, n.folds)
  concordance_testing_list <- rep(NA, n.folds)

  number_retrieved <- data.frame(matrix(nrow = nrow(feature_data), ncol = 1))
  colnames(number_retrieved)[1] <- "Number of Retrieved Neighbors"
  rownames(number_retrieved) <- rownames(feature_data)

  testing_sample_information <- list()
  case_retrieval_information <- list()
  training_folds <- list()
  testing_folds <- list()

  cat(" -- Doing ",n.folds," fold testing using Case Based Reasoning with Confidence-Based Retrieval -- \n")
  for(i in 1:n.folds) {
    fold_case_information <- list()
    # Extracting testing and training fold
    testing_fold <- feature_data[data.folds[[i]],]
    training_fold <- feature_data[-data.folds[[i]],]
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
      neighbors <- list()
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


distance_between <- function(df1, df2) {
  distance_frame <- sqrt(matrix(rowSums(expand.grid(rowSums(df1*df1),rowSums(df2*df2))),
                                nrow=nrow(df1)) - (2 * as.matrix(df1) %*% t(as.matrix(df2))))
  return(as.data.frame(distance_frame))
}

isEmpty <- function(x) {
  return(identical(x, numeric(0)))
}

normalize <- function(x) {
  return ((x - min(x)) / (max(x) - min(x)))
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
