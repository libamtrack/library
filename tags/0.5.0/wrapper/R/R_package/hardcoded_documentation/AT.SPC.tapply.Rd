\name{AT.SPC.tapply}
\alias{AT.SPC.tapply}
\title{AT.SPC.tapply}
\description{Similar to R's \code{tapply} this applies a function to cells defined by indices in a spc data object}
\usage{AT.SPC.tapply( spc, INDEX, FUN, additional.arguments = NULL, names.results = NULL)
}
\arguments{
  \item{spc}{ spc data as returned by \code{\link{AT.SPC.read}}}
  \item{INDEX}{ vector of column names in spc data that should be used as indices (indicating the cells).}
  \item{FUN}{ Function to be applied - preferably this is a libamtrack function obeying the standard naming of mixed field variables \code{\link{number.of.field.components}} and \code{\link{material.no}}. Additional arguments to \code{FUN} can be passed using \code{...})}
  \item{additional.arguments}{ Optional additional arguments to \code{FUN}. This should be a list with a three-entry vector for each argument. The first entry is the argument name, the second the value, the third indicated if it should be applied cell-wise (if TRUE). 
  For example: the primitive R function to evaluate the mean energy per cell: \code{mean(x)}: \code{additional.arguments = list(c("x", "E.MeV.u", TRUE))}
  Or to pass the material.no (same for all cells): \code{additional.arguments = list(c("material.no", "1", FALSE))}}
  \item{names.results}{ optional vector with names for the returned values of \code{FUN}}
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
# Compute the fluence with depth, but differentiate contribution
# from individual particle species 
# As 'mean' does not take any arguments obeying the mixed field naming
# convention in libamtrack, but a vector "x",
# one must specify this as additional argument
# AT.SPC.tapply( spc                  = spc, 
#                INDEX                = c("depth.g.cm2", "particle.no"), 
#                FUN                  = mean,
#                additional.arguments = list(c("x", "fluence.cm2")),
#                names.results        = "fluence.cm2")
}
