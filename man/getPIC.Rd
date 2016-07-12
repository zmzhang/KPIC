\name{getPIC}
\alias{getPIC}
\title{Extracted the pure ion chromatograms via Optimal K-means Clustering}
\description{Return a list which is comprised of PICs and the main information of them}

\usage{getPIC(filename, range, level, alpha = 0.3, Iterate = 4000, n = 3)}

\arguments{
  \item{filename}{path of the LC-MS data file to read, can be.mzXML or .CDF}
  \item{range}{the range of retention time of ROI to be clursted}
  \item{level}{the minimum intensity of a feature, i.e. the minimum intensity of the landmark ion to be considered}
  \item{alpha}{the adjustment factor of ROI in m/z direction}
  \item{Iterate}{the maximum number of iterate computing each part}
  \item{n}{parameter for multi-core computing}
}

\value{
  \item{PICs}{the list of all extracted PICs}
  \item{Info}{the matrix of the properties of the features}
  \item{rt}{the scan time}
}