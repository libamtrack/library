\name{AT.SPC.spectrum.at.depth.g.cm2.Rd}
\alias{AT.SPC.tapply}
\title{AT.SPC.tapply}
\description{Similar to R's \code{tapply} this applies a function to cells defined by indices in a spc data object}
\usage{AT.SPC.tapply( spc, INDEX, FUN, names.results = NULL, ...)
}
\arguments{
  \item{spc}{ spc data}
  \item{INDEX}{ vector of column names in spc data that should be used as indices (indicating the cells).}
  \item{FUN}{ Function to be applied - preferably this is a libamtrack function obeying the standard naming of mixed field variables \code{\link{number.of.field.components}} and \code{\link{material.no}}. Additional arguments to \code{FUN} can be passed using \code{...})}
  \item{names.results}{ optional vector with names for the returned values of \code{FUN}}
  \item{...}{ Optional additional arguments to \code{FUN}}
}
\value{
  A data frame with the following columns:
  \item{index.columns}{ Columns for indices given in \code{INDEX} to be looped over}
  \item{results[]}{ Additional columns containing the returned values from \code{FUN}}
}

\examples{
# Compute the dose with depth from spc data 
#
# AT.SPC.tapply( spc           = spc, 
#                INDEX         = "depth.g.cm2", 
#                FUN           = AT.total.D.Gy, 
#                names.results = "D.Gy")
#
#
# Compute the dose with depth, but differentiate contribution
# from individual particle species 
#
# AT.SPC.tapply( spc           = spc, 
#                INDEX         = c("depth.g.cm2", "particle.no"), 
#                FUN           = AT.total.D.Gy, 
#                names.results = "D.Gy")
}
