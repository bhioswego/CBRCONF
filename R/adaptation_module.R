#' Adaptation Module
#'
#'
#'
#' @param CBRContainer A CBRContainer object
#' @import caret dplyr data.table textshape
#' @export

library(survival)
library(caret)
library(dplyr)
library(data.table)
library(textshape)

adaptation_module <- function(CBRContainer) {
  training_folds <- CBRContainer$training_folds
  testing_folds <- CBRContainer$testing_folds
  case_retrieval_information <- CBRContainer$case_retrieval_information
  testing_sample_information <- CBRContainer$testing_sample_information
  testing_dataframe <- CBRContainer$testing_dataframe

  for(fold in 1:length(training_folds)) {
    training_fold <- training_folds[[fold]]
    testing_fold <- testing_folds[[fold]]
    fold_information <- case_retrieval_information[[fold]]
    sample_information_for_fold <- testing_sample_information[[fold]]



    neighbors <- data.frame(matrix(nrow = nrow(training_fold), ncol = nrow(testing_fold)))
    rownames(neighbors) <- rownames(training_fold)
    colnames(neighbors) <- rownames(testing_fold)

    # This will create a dataframe where the neighbors for a testing case can be easily seen.
    for(training_case_index in 1:nrow(neighbors)) {
      for(testing_case_index in 1:ncol(neighbors)) {
        testing_case <- colnames(neighbors[testing_case_index])
        training_case <- rownames(neighbors)[training_case_index]
        index = which(names(fold_information) == testing_case)
        if(training_case %in% fold_information[[index]]$`retrieved neighbors`) {
          neighbors[training_case_index,testing_case_index] <- "NEIGHBOR"
        } else {
          neighbors[training_case_index,testing_case_index] <- "-"
        }
      }
    } # ADDED
  } # ADDED
      for(fold in 1:length(training_folds)) {
        training_fold <- training_folds[[fold]]
        testing_fold <- testing_folds[[fold]]
        fold_information <- case_retrieval_information[[fold]]
        sample_information_for_fold <- testing_sample_information[[fold]]



        neighbors <- data.frame(matrix(nrow = nrow(training_fold), ncol = nrow(testing_fold)))
        rownames(neighbors) <- rownames(training_fold)
        colnames(neighbors) <- rownames(testing_fold)

        # This will create a dataframe where the neighbors for a testing case can be easily seen.
        for(training_case_index in 1:nrow(neighbors)) {
          for(testing_case_index in 1:ncol(neighbors)) {
            testing_case <- colnames(neighbors[testing_case_index])
            training_case <- rownames(neighbors)[training_case_index]
            index = which(names(fold_information) == testing_case)
            if(training_case %in% fold_information[[index]]$`retrieved neighbors`) {
              neighbors[training_case_index,testing_case_index] <- "NEIGHBOR"
            } else {
              neighbors[training_case_index,testing_case_index] <- "-"
            }
          }
        }

        neighbor_survival_differences <- data.frame(matrix(nrow = nrow(neighbors), ncol = ncol(neighbors)))
        rownames(neighbor_survival_differences) <- rownames(neighbors)
        colnames(neighbor_survival_differences) <- colnames(neighbors)

        # This will create a dataframe of the differences in training survival and the retrieved survival times for the test cases
        for(training_case_index in 1:nrow(neighbor_survival_differences)) {
          for(testing_case_index in 1:ncol(neighbor_survival_differences)) {
            testing_case <- colnames(neighbor_survival_differences[testing_case_index])
            training_case <- rownames(neighbor_survival_differences)[training_case_index]
            training_index <- which(rownames(training_fold) == training_case)
            testing_index <- which(rownames(testing_dataframe) == testing_case)
            if(neighbors[training_case_index, testing_case_index] == "NEIGHBOR") {
              training_survival <- training_fold[training_index,]$time
              testing_survival <- testing_dataframe[testing_index,]$time
              time_difference <- training_survival - testing_survival
              neighbor_survival_differences[training_case_index, testing_case_index] <- time_difference
            } else {
              neighbor_survival_differences[training_case_index, testing_case_index] <- 0
            }
          }
        }

        # This will create a table of the average differences in survival time
        average_differences <- as.data.frame(rowMeans(neighbor_survival_differences))
        colnames(average_differences)[1] <- "Average Difference"

        new_training_fold <- training_fold
        new_training_fold$time <- training_fold$time + average_differences$`Average Difference`



        testing_features <- testing_fold[,3:ncol(testing_fold)]
        training_features <- new_training_fold[,3:ncol(new_training_fold)]

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
        new_confidence_weights <- as.data.frame(confidence_temp$`New Score`)
        rownames(new_confidence_weights) <- rownames(confidence_temp)
        colnames(new_confidence_weights)[1] <- "Confidence Weight"
        new_confidence_weights <- as.data.frame(normalize(new_confidence_weights))


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
            confidence <- new_confidence_weights[which(rownames(new_confidence_weights) == neighbor),]
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
          fold_information[[k]] <- single_case_information
          names(fold_information)[[k]] <- query_case
          number_retrieved[which(rownames(number_retrieved) == query_case),] <- length(single_case_information[[1]])
        }

        testing_neighbor_times$`risk group` <- risk_groups[match(rownames(testing_fold), rownames(risk_groups)), ]
        testing_statistics <- coxph(Surv(testing_neighbor_times$time, testing_neighbor_times$status) ~ testing_neighbor_times$`risk group`, data = testing_features)



        testing_sample_information[[i]] <- testing_neighbor_times


        testing_set_concordance <- testing_statistics[["concordance"]]["concordance"]
        concordance_testing_list[i] <- testing_set_concordance[1]

        case_retrieval_information[[i]] <- fold_information
        names(case_retrieval_information)[[i]] <- paste("fold", i)

      }

      testing_concordance <- mean(concordance_testing_list)
      cat("Testing concordance is: ",testing_concordance,"\n")
}
