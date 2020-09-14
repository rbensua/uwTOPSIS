TOPSIS <- function(x,
                   directions,
                   w,
                   norm.method = c("norm", "gauss", "minmax"),
                   ordered = FALSE,
                   makefigure = TRUE){

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

  # Cheking weights

  if(any(w<0) | abs(sum(w) - 1) > 1e-8 ){
    stop("Weights should be positive and add up to 1!")
  }

  # Checking normalization method
  norm.method <- match.arg(norm.method)

  ##
  ## DEFINING THE NORMALIZATION FUNCTION
  ##


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
    } else {
      norm_x <- apply(x, MARGIN = 2, function (y){
        res <- (y - mean(y)) / sd(y)
        return(unname(res))
      })
    }
    return(norm_x)
  }


  # Normalizing data

  altnames <- apply(x, MARGIN = 1, FUN = function(y) y[1])
  NormMat <- normalize(x[,-c(1)], method = norm.method)

  rownames(NormMat) <-   altnames



  #
  # DEFINING OBJECTIVE DIRECTIONS
  #



  obj_pos <- directions
  obj_neg <- ifelse(obj_pos == "max", "min","max")

  mask_pos_max <- obj_pos == "max"
  mask_pos_min <- obj_pos == "min"
  mask_neg_max <- obj_neg == "max"
  mask_neg_min <- obj_neg == "min"

  #
  #  DETERMINING THE IDEAL AND ANTI-IDEAL SOLUTIONS
  #
  #

    idx_max <- max.col(t(NormMat))
    idx_min <- max.col(-t(NormMat))
    maximums <- diag(NormMat[idx_max,])
    minimums <- diag(NormMat[idx_min,])


  opt_pos <- maximums*mask_pos_max + minimums*mask_pos_min
  opt_neg <- maximums*mask_neg_max + minimums*mask_neg_min

  ##
  ## DEFINING THE DISTANCE FUNCTION (EUCLIDEAN)
  ##


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

  ##
  ## DEFINING THE SCORE FUNCTION
  ##

  score <- function(w,i) {
    r <- D(w, direction = "neg",i) / (D(w, direction = "neg",i) + D(w, direction = "pos",i))
    return(r)
  }

  ##
  ## COMPUTING ALL THE SCORES
  ##
  ##

  TOPSIS_scores <- sapply(seq_along(altnames), function(k) score(w,k))
  names(TOPSIS_scores) <- altnames

  #
  # MERGING THE RESULTS IN A DATA.FRAME
  #

  scores_DF <- cbind(x[,1],
                     data.frame(TOPSIS = TOPSIS_scores,
                                row.names = NULL
                     )
  )
  colnames(scores_DF)[1] <- colnames(x)[1]

  if(ordered){
    scores_DF <- scores_DF[sort(scores_DF$TOPSIS, decreasing = TRUE, index.return = TRUE)$ix,]
  }


  #
  # PLOT FIGURE
  #


  if(makefigure){
    scores_DF$Idx <- 1:nrow(scores_DF)



    p <-  ggplot(data = scores_DF, aes(x = Idx)) +
      geom_line(aes( y = TOPSIS)) +
      geom_point(aes( y = TOPSIS)) +
      theme_classic() +
      scale_x_continuous(breaks = 1:nrow(scores_DF), labels = scores_DF[,1]) +
      xlab(colnames(scores_DF)[1]) + ylab("Score")
    show(p)
  }
}
