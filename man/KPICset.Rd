\name{KPICset}
\alias{KPICset}
\title{Extracted the pure ion chromatograms for a set of samples}
\description{Return a KPICset object}
\usage{
KPICset(files, range = 10, level = 500, alpha = 0.2, Iterate = 4000, n = 3, refine = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{files}{path names of the NetCDF/mzXML files to read}
  \item{range}{the range of retention time of ROI to be clursted}
  \item{level}{the minimum intensity of a feature, i.e. the minimum intensity of the landmark ion to be considered}
  \item{alpha}{the adjustment factor of ROI in m/z direction}
  \item{Iterate}{the maximum number of iterate computing each part}
  \item{n}{parameter for multi-core computing}
  \item{refine}{the parameter of the PICrefine function, if refine = NA the procedure of refinement is ignored}
}

\value{
   A KPICSet object.
   \item{@PICset}{the list of extracted PICs of each sample}
   \item{@peakmat}{the matrix of the properties of the features}
   \item{@sample}{the names of each sample}
}