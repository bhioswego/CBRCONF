#' Create a CBRCONF object
#'
#' CBRCONF is an R6 class container that holds the dataset and relevant results at each stage.
#' @param dataset A dataset of features, where the first column is the time and the second column is the event statuses. Remaining columns are features.
#' @export

create_CBRCONF <- function(dataset) {
  object <- CBRCONF$new(dataset)
  cat(" -- Created -- \n")
  cat(" -- The next module to run is the Cox Module -- \n")
  return(object)
}
