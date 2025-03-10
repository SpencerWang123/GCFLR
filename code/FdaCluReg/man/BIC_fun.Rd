% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BIC_fun.R
\name{BIC_fun}
\alias{BIC_fun}
\title{BIC loss result to determine the number of clusters}
\usage{
BIC_fun(reg_dt, label, temp_K_seq, method = "ls")
}
\arguments{
\item{reg_dt}{N*M principle components scores of the observations}

\item{label}{N*1 label vector that conveys clustering results}

\item{method}{loss function used in the regression}

\item{K}{number of clusters}
}
\value{
calculated BIC results
}
\description{
BIC loss result to determine the number of clusters
}
