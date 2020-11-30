#' KNN Module
#'
#' This module will use K-nearest neighbor to predict survival.
#'
#'
#' @param CBRContainer A CBRContainer object
#' @param n.folds The number of folds
#' @param KM.Plot If KM.Plot is true, a kaplan-meier plot will be constructed from 3 risk groups (low, medium and high)
#' @import caret
#' @examples
#' @export

protoKNN <- function(CBRContainer, n.folds = 10, KM.Plot = TRUE) {
  confidence_weights <- CBRContainer$confidence_weights

  sig_feature_data <- CBRContainer$sig_feature_data
  prognostic_scores <- CBRContainer$prognostic_scores
  risk_groups <- CBRContainer$risk_groups

  set.seed(25)
  data.folds <- createFolds(y=sig_feature_data$time, k = n.folds)

  concordance_training_list <- rep(NA, n.folds)
  concordance_testing_list <- rep(NA, n.folds)

  number_retrieved <- data.frame(matrix(nrow = nrow(sig_feature_data), ncol = 1))
  colnames(number_retrieved)[1] <- "Number of Retrieved Neighbors"
  rownames(number_retrieved) <- rownames(sig_feature_data)


  testing_sample_information <- list()


  for(i in 1:n.folds) {
    # Extracting testing and training fold
    testing_fold <- sig_feature_data[data.folds[[i]],]
    training_fold <- sig_feature_data[-data.folds[[i]],]

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
      current_confidence <- 0.0


        for(n in 1:length(name_of_neighbors)) {
          neighbor <- name_of_neighbors[n]
          confidence <- confidence_weights[which(rownames(confidence_weights) == neighbor),]
          survival_time <- sig_feature_data[which(rownames(sig_feature_data) == neighbor), ]
          survival_time <- survival_time$time
          single_case_information[[1]] <- append(single_case_information[[1]], neighbor)
          single_case_information[[2]] <- append(single_case_information[[2]], survival_time)
          prognostic_score <- prognostic_scores[which(rownames(prognostic_scores) == neighbor), ]
          single_case_information[[3]] <- append(single_case_information[[3]], prognostic_score)
          risk_group <- risk_groups[which(rownames(risk_groups) == neighbor), ]
          single_case_information[[4]] <- append(single_case_information[[4]], risk_group)
          current_confidence <- current_confidence + confidence
          if(current_confidence > 1.0) {
            break
          }
          testing_neighbor_times[which(rownames(testing_neighbor_times) == query_case), 1] <- mean(unlist(single_case_information$`survival times`))
          testing_neighbor_times[which(rownames(testing_neighbor_times) == query_case), 3] <- prognostic_scores[which(rownames(prognostic_scores) == query_case),]
          #retrieved_risk_groups <- unlist(single_case_information[[4]])
          #risk_table <- table(retrieved_risk_groups)
          #group <- sample(names(which(risk_table == max(risk_table))), size = 1)
          #testing_neighbor_times[which(rownames(testing_neighbor_times) == query_case), 4] <- group
        }
      number_retrieved[which(rownames(number_retrieved) == query_case),] <- length(single_case_information[[1]])
    }

    # Use this for creating 3 groups from averaged prognostic scores
    #testing_neighbor_times$`risk group` <- ile(testing_neighbor_times$`prognostic scores`, 3)
    # Using the assigned risk group after cox regression for the 100 training sets
    testing_neighbor_times$`risk group` <- risk_groups[match(rownames(testing_fold), rownames(risk_groups)), ]
    #testing_neighbor_times$`risk group` <- ile(testing_neighbor_times$time)
    testing_statistics <- coxph(Surv(testing_neighbor_times$time, testing_neighbor_times$status) ~ testing_neighbor_times$`risk group`, data = testing_features)

    testing_sample_information[[i]] <- testing_neighbor_times

    testing_set_concordance <- testing_statistics[["concordance"]]["concordance"]
    concordance_testing_list[i] <- testing_set_concordance[1]


}


  if(KM.Plot == TRUE) {
    cat("-- Writing KM Plot for Fold ",i," -- \n")
    SurvTimes <- Surv(sig_feature_data$time, sig_feature_data$status)
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
    cat("-- Writing the Kaplan Meier Plot to KM-plot.png -- \n")
    png(filename = paste("KM-Plot.png"), width = 5.5, height = 5.5, units = "in", res = 300, pointsize = 7)
    plot(fit, mark.time=TRUE, xlab = "Months", ylab = "Survival", main=paste("Survival Curve",sep=""),lty = 1:3, col = 1:3, cex = 0.5)
    grid()
    legend(x = "bottomright", legend = c(leg1, leg2, leg3), lty = 1:3,
           col = 1:3, cex = 0.65)
    text(10,0.1,paste("Log-Rank \n p = ", formatC(p, format="g", digits = 3), sep = ""),
         pos = 4, cex = 1)
    dev.off()
  }

  testing_concordance <- mean(concordance_testing_list)
  cat("Testing concordance is: ",testing_concordance,"\n")
}
