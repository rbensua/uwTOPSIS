#' @title Unweighted TOPSIS
#'
#' @description Computes the unweighted TOPSIS given the performance table.
#'
#' @usage uwTOPSIS(x, directions, norm.method = c("norm", "gauss", "minmax"), L = NULL, U = NULL, w0 = NULL)
#'
#' @param x Dataframe with the performances of each alternative at each criterion. The first column should be the alternatives definition, the subsequent columns correspond to the different criteria.
#' @param norm.method Character string. Normalization method. Either "norm", "gauss" or "minmax".
#' @param L numeric. Vector containing the lower bound for the weights of the criteria. If NULL (default) it will be zero.
#' @param U numeric. Vector containing the upper bound for the weights of the criteria. If NULL (default) it will be one.
#' @param w0 numeric. Vector containing the initial guess for the optimal weights of the criteria. If NULL (default) it will be (L+ U) / 2.
#' @param ordered. If TRUE the resulting table is ordered with respect the average TOPSIS score (descendent order). If FALSE (default) the resulting table is given in the same order as the input performance table.
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
#' @import nloptr
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot geom_line scale_x_continuous theme_classic xlab ylab aes  geom_point geom_ribbon scale_linetype_manual scale_color_manual
#' @importFrom methods show
#' @export
uwTOPSIS <- function(x,
                     directions,
                     norm.method = c("norm", "gauss", "minmax"),
                     L = NULL,
                     U = NULL,
                     w0 = NULL,
                     ordered = FALSE,
                     makefigure = TRUE){

  # Checking arguments

  # Checking input data

  if(!is.data.frame(x)){
    stop("Input data should be a data frame")
  }else {
    x <- as.data.frame(x)
  }

  # Checking directions

  if (length(directions) != ncol(x)-1){
    stop("Number of optimal directions should be equal to the number of criteria!")
  }

  if ( any( !directions %in% c("min","max"))){
    stop("Directions should be either 'max' or 'min'!")
  }

  # Checking normalization method
  norm.method <- match.arg(norm.method)

  # Checking weight bounds

  if(is.null(L)){
    L <- rep(0,ncol(x)-1)
  } else{
    if (length(L) != ncol(x) -1 ){
      stop("Number of lower bounds should be equal to the number of criteria!")
    } else if (any(L < 0) | any(L > 1) ){
      stop("Lower bounds should be between 0 and 1!")
    }
  }

  if(is.null(U)){
    U <- rep(1,ncol(x)-1)
  } else{
    if (length(U) != ncol(x) -1 ){
      stop("Number of upper bounds should be equal to the number of criteria!")
    } else if (any(U < 0) | any(U > 1) ){
      stop("Upper bounds should be between 0 and 1!")
    }
  }

  # Checking initial weights

  if(is.null(w0)){
    w0 <- (L + U) / 2
  } else if (length(w0) != ncol(x) -1){
    stop("Number of upper bounds should be equal to the number of criteria!")
  } else if (any(w0 < L) | any(w0 > U)) {
    stop("Initial bounds should be between the lower and upper bounds!")
  }


  # Normalizing data

  NormMat <- normalize(x[,-c(1)], method = norm.method)
  rownames(NormMat) <- x[,c(1)]


  #
  # DEFINING OBJECTIVE DIRECTIONS
  #



  obj_pos <- directions
  obj_neg <- ifelse(obj_pos == "max", "min","max")

  mask_pos_max <- obj_pos == "max"
  mask_pos_min <- obj_pos == "min"
  mask_neg_max <- obj_neg == "max"
  mask_neg_min <- obj_neg == "min"

  idx_max <- max.col(t(NormMat))
  idx_min <- max.col(-t(NormMat))

  #
  # DETERMINING THE IDEAL AND ANTI-IDEAL SOLUTIONS
  #

  maximums <- diag(NormMat[idx_max,])
  minimums <- diag(NormMat[idx_min,])

  opt_pos <- maximums*mask_pos_max + minimums*mask_pos_min
  opt_neg <- maximums*mask_neg_max + minimums*mask_neg_min


  #
  # DEFINING THE SCORE FUNCTION AND ITS GRADIENT
  #



  D <- function(w,direction = c("pos","neg"),i){
    direction <- match.arg(direction)
    if(direction == "pos"){
      opt <- opt_pos
    } else{
      opt <- opt_neg
    }
    Dist <- sqrt(sum((w*NormMat[i,]-w*opt)^2))
    return(Dist)
  }

  gradD <- function(w,direction = c("pos","neg"),i){
    direction <- match.arg(direction)
    opt <- ifelse(direction == "pos", opt_pos, opt_neg)
    num <- w * (NormMat[i,] - opt)^2
    denom <- D(w, direction = direction)
    return(num/denom)
  }

  score <- function(w,i,case = c("lower","upper")) {
    case <- match.arg(case)
    r <- D(w, direction = "neg",i) / (D(w, direction = "neg",i) + D(w, direction = "pos",i))
    r <- ifelse(case == "upper", -r,r)
    return(r)
  }

  gradscore <- function(w,i,case = c("lower","upper")) {
    case <- match.arg(case)
    grad <- numeric(length = length(w))

    prefactor <- 1/(D(w, direction = "pos", i) + D(w, direction = "neg", i))^2

    matrow1 <- c(D(w, direction = "pos", i) ,  D(w, direction = "neg", i))

    gD_pos <- gradD(w,direction = "pos")
    gD_neg <- gradD(w,direction = "neg")

    for (k in seq_along(w)){
      matrow2 <- c(gD_pos[k], gD_neg[k])
      mat <- matrix(c(matrow1,matrow2), nrow = 2, byrow = TRUE)
      grad[k] <- det(mat)
    }
    grad <- prefactor * grad

    sign_case <- ifelse(case == "upper", -1, 1)

    return(sign_case*grad)
  }


  #
  # RESTRICTIONS
  #
  eval_g0_ineq <- function(w,i,case){
    constr <- c(sum(w)-1,
                1-sum(w))
    return(constr)
  }

  grad_eval_g0_ineq <- function(w,i,case){
    res1 <- rep(1,length(w))
    res2 <- rep(-1,length(w))
    return(rbind(res1,res2))
  }









  #
  # OPTIMIZING
  #

  N <- nrow(NormMat)
  m <- length(obj_pos)


  # In this matrices we will store the weights
  solutions_min <- matrix(nrow = N, ncol = m) #
  solutions_max <- matrix(nrow = N, ncol = m) #

  # In this vectors we will store the scores
  score_min <- numeric(N)
  score_max <- numeric(N)
  names(score_min) <- rownames(NormMat)
  names(score_max) <- rownames(NormMat)
  for(i in 1:N){
    sols <- nloptr(x0 = w0,
                   eval_f =  score,
                   #eval_grad_f = gradscore, # COBYLA IS GRADIENT-FREE
                   eval_g_ineq =  eval_g0_ineq,
                   #eval_jac_g_ineq = grad_eval_g0_ineq,
                   lb = L,
                   ub = U,
                   opts = list("algorithm" = "NLOPT_LN_COBYLA",
                               "xtol_rel" = 1e-27,
                               "xtol_abs" = 1e-27,
                               "maxeval" = 2000),
                   case = "lower", # lower = minimizing score, upper = maximizing score
                   i = i) # i is the number of alternative
    solutions_min[i,] <- sols$solution
    score_min[i] <- sols$objective
  }

  for(i in 1:N){
    sols <- nloptr(x0 = w0,
                   eval_f =  score,
                   # val_grad_f = gradscore, # COBYLA IS GRADIENT-FREE
                   eval_g_ineq =  eval_g0_ineq,
                   #eval_jac_g_ineq = grad_eval_g0_ineq,
                   lb = L,
                   ub = U,
                   opts = list("algorithm" = "NLOPT_LN_COBYLA",
                               "xtol_rel" = 1e-15,
                               "xtol_abs" = 1e-15,
                               "maxeval" = 2000),
                   case = "upper", # lower = minimizing score, upper = maximizing score
                   i = i) # i is the number of alternative
    solutions_max[i,] <- sols$solution
    score_max[i] <- -sols$objective
  }

#
# MERGING THE RESULTS IN A DATA.FRAME
#

  scores_DF <- cbind(x[,1],
                     data.frame(Min = score_min,
                                Max = score_max,
                                uwTOPSIS = 0.5*(score_min+score_max),
                                row.names = NULL
                                )
                     )
  colnames(scores_DF)[1] <- colnames(x)[1]

  if(ordered){
    scores_DF <- scores_DF[sort(scores_DF$uwTOPSIS, decreasing = TRUE, index.return = TRUE)$ix,]
  }

#
# PLOT FIGURE
#


  if(makefigure){
    scores_DF$Idx <- 1:nrow(scores_DF)

    df <- pivot_longer(scores_DF, -c(1,ncol(scores_DF)), names_to = "Score", values_to = "Value")

   p <-  ggplot(data = df, aes(x = Idx)) +
      geom_ribbon(data = scores_DF, aes(ymin = Min, ymax = Max), fill = "gray", alpha =0.5) +
      geom_line(aes( y = Value, col = Score, group = Score, lty = Score)) +
      scale_linetype_manual(values=c("solid","dotted","dashed")) +
      scale_color_manual(values = c("red", "black","black")) +
      geom_point(aes( y = Value, colour = Score, group = Score)) +
      theme_classic() +
      scale_x_continuous(breaks = 1:nrow(scores_DF), labels = scores_DF[,1]) +
      xlab(colnames(scores_DF)[1]) + ylab("Score")
   show(p)
  }


return(list(scores = scores_DF, weights_min = solutions_min, weights_max = solutions_max) )
}

