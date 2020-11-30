#' Distance between Two Dataframes
#'
#' This function takes two dataframes and calculates the Euclidean distances between their columns.
#' @param df1 The first dataframe
#' @param df2 The second dataframe
#' @export

distance_between <- function(df1, df2) {
  distance_frame <- sqrt(matrix(rowSums(expand.grid(rowSums(df1*df1),rowSums(df2*df2))),
                 nrow=nrow(df1)) - (2 * as.matrix(df1) %*% t(as.matrix(df2))))
  return(as.data.frame(distance_frame))
}

distance_between2 <- function(df1, df2) {
  distance_frame <- sqrt(matrix(rowSums(expand.grid(rowSums(df1*df1),rowSums(df2*df2))),
                                nrow=nrow(df1)) - (2 * t(as.matrix(df2) %*% as.matrix(df1))))
  return(as.data.frame(distance_frame))
}
