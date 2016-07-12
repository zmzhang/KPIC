\name{PICrefine}
\alias{PICrefine}
\title{
Refine the PICs
}
\description{
Remove unsatisfactory PICs by calculating the Gaussian similarity, sharpness and SNR
}
\usage{
PICrefine(PICs, n = 1)
}

\arguments{
  \item{PICs}{
the result of the extracting through getPIC function
}
  \item{n}{
a parameter for calculating the threshold i.e. the lambda referred in the article
}}

\value{
  \item{PICs}{the list of all extracted PICs}
  \item{Info}{the matrix of the properties of the features}
  \item{rt}{the scan time}
  \item{scores}{the score of each PIC}
}
