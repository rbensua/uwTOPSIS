% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{GShotels}
\alias{GShotels}
\title{Data: A. Oztel, A. Karakas, A. Ahmet (2018)}
\format{
Data frame (tibble) with 35 rows and 11 columns. The variables are:
\describe{
  \item{HOTEL}{Code of the hotel.}
  \item{dSubway}{Distance to the closest subway station (km).}
  \item{dSubway}{Distance to the Airport (km).}
  \item{dCoachStation}{Distance to Coach Station (km).}
  \item{dSutanSquare}{Distance to Sultan Ahmet Square (km).}
  \item{Traffic}{Traffic intensity (1-4).}
  \item{PriceTA.}{Price given in Trip Advisor platform }
  \item{scTrivago}{Score given in Trivago platform (0-10)}
  \item{scTripAdvisor}{Score given in Trip Advisor platform (1-5)}
  \item{scHotels.com}{Score given in Hotels.com platform (0-10)}
  \item{scGoogle}{Score given in Google.com (1-5).}
}
}
\source{
Oztel A, Karakas A, Ahmet A. 2018. Green Star Hotel Selection with Fuzzy TOPSIS Method:Case of Istanbul. In: Dorczak R, Arslan H, Musialik R, editors. Recent researches on socialsciences. First edit ed. Krakow: Institute of Public Affairs - Jagiellonian University; chap.Green Star; p. 473–485.
}
\usage{
data("GShotels")
}
\description{
Decision matrix for Green Star hotel selection problem in Turkey.
}
\author{
\strong{Vicente Liern Carrion} (\email{vicente.liern@uv.es}).
\emph{Department of Business Mathematics}

\strong{Rafael Benitez} (\email{rafael.suarez@uv.es}).
\emph{Department of Business Mathematics}

University of Valencia (Spain)
}
\keyword{datasets}
