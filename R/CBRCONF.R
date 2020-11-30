#' CBRCONF Class
#'
#' This is a class object to store all of the results performed with the modules of this package
#' @import R6

CBRCONF <- R6Class(classname = "CBRCONF",
                  public = list(
                    dataset = NULL,
                    #' @field dataset The dataset where time and status are the first and second columns and subsequent columns are features
                    sig_feature_list = NULL,
                    #' @field sig_feature_list List of all significant features found in each iteration of feature selection
                    sig_feature_data = NULL,
                    #' @field sig_feature_data Significant features that are retained after all iterations of feature selection
                    feature_scores = NULL,
                    #' @field feature_scores A dataframe to contain all of the individual feature scores. This is the feature's expression value multiplied by the mean beta coefficient
                    prognostic_scores = NULL,
                    #' @field prognostic_scores Dataframe of prognostic scores for all samples
                    results_list = NULL,
                    #' @field results_list Holds the results for feature selection
                    risk_groups = NULL,
                    #' @field risk_groups Holds the risk groups, which is assigned by the prognostic risk scores
                    group_means = NULL,
                    #' @field group_means The features of each sample within each risk group is averaged to create one prototypical representation of each group. This stores that information
                    low_risk_features = NULL,
                    #' @field low_risk_features A subset of the original dataset that is just the low risk group, and their corresponding features
                    medium_risk_features = NULL,
                    #' @field medium_risk_features A subset of the original dataset that is just the medium risk group, and their corresponding feature
                    high_risk_features = NULL,
                    #' @field high_risk_features A subset of the original dataset that is just the high risk group, and their corresponding feature
                    low_risk_prototype = NULL,
                    #' @field low_risk_prototype A vector of the averaged features of the low risk group
                    medium_risk_prototype = NULL,
                    #' @field medium_risk_prototype A vector of the averaged features of the medium risk group
                    high_risk_prototype = NULL,
                    #' @field high_risk_prototype A vector of the averaged features of the high risk group
                    confidence_weights_low = NULL,
                    #' @field confidence_weights_low The confidence weights for just the low risk group
                    confidence_weights_medium = NULL,
                    #' @field confidence_weights_medium The confidence weights for just the medium risk group
                    confidence_weights_high = NULL,
                    #' @field confidence_weights_high The confidence weights for just the high risk group
                    confidence_weights = NULL,
                    #' @field confidence_weights A dataframe that combines the confidence weights of the three risk groups
                    data.folds = NULL,
                    #' @field data.folds Each of the training and testing folds
                    training_folds = NULL,
                    #' @field training_folds Each of the training folds
                    testing_folds = NULL,
                    #' @field testing_folds Each of the testing folds
                    case_retrieval_information = NULL,
                    #' @field case_retrieval_information Which cases were retrieved by which case
                    testing_number_retrieved = NULL,
                    #' @field testing_number_retrieved How many cases each testing sample retrieved
                    testing_concordance = NULL,
                    #' @field testing_concordance The concordance values for each testing fold
                    testing_sample_information = NULL,
                    #' @field testing_sample_information The risk group, predicted survival times, status, and testing_results
                    #' @description
                    #' Create a new CBRContainer object.
                    #' @param dataset
                    #' @return A new 'CBRMSR' object.
                    initialize = function(dataset) {
                      self$dataset <- dataset
                    }
                  )
)
