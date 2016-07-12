\name{getTIC}
\alias{getTIC}
\title{get TICs of samples}
\description{Return a list which is comprised of TICs}

\usage{getTIC(files,method='TIC')}

\arguments{
  \item{files}{path of the LC-MS data file to read, can be.mzXML or .CDF}
  \item{method}{'TIC'or'BPC'}
}

\value{
  \item{TICs}{the list of TICs}
}