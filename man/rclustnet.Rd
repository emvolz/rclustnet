% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rclustnet.R
\name{rclustnet}
\alias{rclustnet}
\title{Generate a random graph with desired size, clustering coefficient, and degree distribution}
\usage{
rclustnet(n, CC, rdegdist = "pois", ...)
}
\arguments{
\item{n}{Integer size of random graph}

\item{CC}{Desired clustering coefficient. Value in (0,1).}

\item{rdegdist}{This is either a character string specifiying the desired degree distribution or an R function which generates the desired distribution.}
}
\value{
A graph in igraph format
}
\description{
Generate a random graph with desired size, clustering coefficient, and degree distribution
}
\examples{
rclustnet( n=500, CC = .25, rdegdist='pois', lambda=5 )
}
