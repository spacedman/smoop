\docType{package}
\name{smoop}
\alias{smoop}
\alias{smoop-package}
\title{smoop}
\usage{
  smoop(y, n, spdata, bounds, M, nx = 64, ny = 64)
}
\arguments{
  \item{y}{formula for variable of interest}

  \item{n}{formula for population}

  \item{spdata}{spatial data frame}

  \item{region}{for bounds and clipping}

  \item{M}{minimal population size}

  \item{nx}{grid size x}

  \item{ny}{grid size y}
}
\value{
  a raster stack with values and SEs and something
}
\description{
  smoop

  Do a spatial smoothing where the baseline ratio is over a
  minimum population size.
}

