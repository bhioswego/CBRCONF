#' KNN Module
#'
#' This module will use K-nearest neighbor to predict survival.
#'
#'
#' @param CBRContainer A CBRContainer object
#' @param n.folds The number of folds
#' @param KM.Plot If KM.Plot is true, a kaplan-meier plot will be constructed from 3 risk groups (low, medium and high)
#' @import caret dplyr data.table textshape
#' @examples
#' @export

protoKNN <- function(CBRContainer, n.folds = 10, KM.Plot = TRUE) {
  confidence_weights <- CBRContainer$confidence_weights

  sig_feature_data <- CBRContainer$sig_feature_data
  prognostic_scores <- CBRContainer$prognostic_scores
  risk_groups <- CBRContainer$risk_groups


  set.seed(25)
  data.folds <- createFolds(y=sig_feature_data$time, k = n.folds)
  #CBRContainer$data.folds <- data.folds


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
        current_confidence <- current_confidence + confidence
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

  CBRContainer$case_retrieval_information <- case_retrieval_information
  CBRContainer$testing_number_retrieved <- number_retrieved
  CBRContainer$testing_concordance <- testing_concordance
  CBRContainer$testing_sample_information <- testing_sample_information
  CBRContainer$training_folds <- training_folds
  CBRContainer$testing_folds <- testing_folds
  CBRContainer$testing_results <- testing_dataframe

  return(CBRContainer)
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
