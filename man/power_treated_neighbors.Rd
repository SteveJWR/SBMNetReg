% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network_features.R
\name{power_treated_neighbors}
\alias{power_treated_neighbors}
\title{Power of treated neighbors}
\usage{
power_treated_neighbors(data, G, treatment_column = "A", power = 2)
}
\arguments{
\item{data}{data frame}

\item{G}{Adjacency matrix}

\item{treatment_column}{Which column name indicates the treatment}

\item{power}{Which power to raise the adjacency matrix to}
}
\value{
number of treated neighbours in the power graph
}
\description{
Power of treated neighbors
}
