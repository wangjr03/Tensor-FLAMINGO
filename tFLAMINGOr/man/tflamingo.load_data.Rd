% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tflamingo.load_data.R
\name{tflamingo.load_data}
\alias{tflamingo.load_data}
\title{tflamingo.load_data}
\usage{
tflamingo.load_data(input_PD, input_IF, chr_name)
}
\arguments{
\item{input_PD}{Path to the completed pair-wise distance matrix.}

\item{input_IF}{Path to the completed interaction frequency matrix.}

\item{chr_name}{Name of the chromosome}
}
\value{
A flamingo data object for later use
}
\description{
load the completed scHi-C matrices and make flamingo data obj
}
\keyword{tFLAMINGO}
