% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gutils.r
\name{gintervals.neighbors1}
\alias{gintervals.neighbors1}
\title{Finds neighbors between two sets of intervals (and does not return conflicting column names)}
\usage{
gintervals.neighbors1(intervals1 = NULL, intervals2 = NULL,
  maxneighbors = 1, mindist = -1e+09, maxdist = 1e+09,
  na.if.notfound = TRUE)
}
\arguments{
\item{intervals1}{intervals}

\item{intervals2}{intervals}

\item{maxneighbors}{maximal number of neighbors}

\item{mindist}{distance range for 1D intervals}

\item{maxdist}{distance range for 1D intervals}

\item{na.if.notfound}{if 'TRUE' return 'NA' interval if no matching neighbors were found, otherwise omit the interval in the answer}
}
\description{
Finds neighbors between two sets of intervals (and does not return conflicting column names)
}
