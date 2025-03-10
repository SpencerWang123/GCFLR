% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distribution_parameters.R
\name{distribution_parameters}
\alias{distribution_parameters}
\title{multinomial parameters for the full condition distribution at the mask position}
\usage{
distribution_parameters(
  reg_dt,
  estimated_parameters,
  mask_index,
  K,
  lambda = 1,
  method = "ls"
)
}
\arguments{
\item{reg_dt}{N*M principle components scores of the observations}

\item{estimated_parameters}{coefficient estimation of the clustering regression results}

\item{mask_index}{int in 1:N used to mask several specific label (for full condition distribution estimation)}

\item{K}{number of clusters}

\item{lambda}{parameters used to construct a multinomial distribution}

\item{method}{loss function used in the regression}
}
\value{
multinomial parameters for the full condition distribution at the mask position
}
\description{
multinomial parameters for the full condition distribution at the mask position
}
