% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/uwTOPSIS.R
\name{uwTOPSIS}
\alias{uwTOPSIS}
\title{Unweighted TOPSIS}
\usage{
uwTOPSIS(x, directions, norm.method = c("norm", "gauss", "minmax"), L = NULL, U = NULL, w0 = NULL)
}
\arguments{
\item{x}{Dataframe with the performances of each alternative at each criterion. The first column should be the alternatives definition, the subsequent columns correspond to the different criteria.}

\item{norm.method}{Character string. Normalization method. Either "norm", "gauss" or "minmax".}

\item{L}{numeric. Vector containing the lower bound for the weights of the criteria. If NULL (default) it will be zero.}

\item{U}{numeric. Vector containing the upper bound for the weights of the criteria. If NULL (default) it will be one.}

\item{w0}{numeric. Vector containing the initial guess for the optimal weights of the criteria. If NULL (default) it will be (L+ U) / 2.}

\item{forceideal}{logical. Forcing the ideal solution. If this parameter is TRUE, the normalized ideal solution is forced to be 1s for the maximizing criteria and 0 for the minimizing criteria. The corresponding antiideal solution would be its oposite. Defaults to FALSE.}

\item{ordered.}{If TRUE the resulting table is ordered with respect the average TOPSIS score (descendent order). If FALSE (default) the resulting table is given in the same order as the input performance table.}
}
\description{
Computes the unweighted TOPSIS given the performance table.
}
\examples{

x <- matrix(1:16, nrow = 4)
normalize(x)

}
\author{
\strong{Rafael Benítez} (\email{rafael.suarez@uv.es}).
\emph{Department of Business Mathematics}

\strong{Vicente Liern} (\email{vicente.liern@uv.es}).
\emph{Department of Business Mathematics}

University of Valencia (Spain)
}
