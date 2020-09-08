#' @title Normalization
#'
#' @description Normalization of the performance matrix.
#'
#' @usage normalize(x, method = c("norm", "gauss", "minmax","none"))
#'
#' @param x Matrix with the performances of each alternative at each criterion
#' @param method Character string. Normalization method. Either "norm", "gauss", "minmax" or "none".
#' If "none" is selected, data matrix given in x is assumed to be already normalized and no further operations are carried out.
#' @author
#'
#' \strong{Rafael Benítez} (\email{rafael.suarez@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' \strong{Vicente Liern} (\email{vicente.liern@@uv.es}).
#' \emph{Department of Business Mathematics}
#'
#' University of Valencia (Spain)
#' @examples
#'
#' x <- matrix(1:16, nrow = 4)
#' normalize(x)
#'
#' @export


normalize <- function(x, method = c( "norm", "gauss", "minmax")){
  method <- match.arg(method)
  if (method == "norm"){
    norm_x <-  apply(x, MARGIN = 2, function(y) {
      res <- y/sqrt(sum(y^2))
      return(unname(res))
    }
    )
  } else if (method == "minmax"){
    norm_x <- apply(x, MARGIN = 2, function (y){
      res <- (y - min(y)) / diff(range(y))
      return(unname(res))
    })
  } else if (method == "gauss") {
    norm_x <- apply(x, MARGIN = 2, function (y){
      res <- (y - mean(y)) / sd(y)
      return(unname(res))
    })
  } else {
    norm_x <- x
  }
  return(norm_x)
}
