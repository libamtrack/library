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
  \item{mixed.field.arguments}{ Named list containing the variables necessary to describe a mixed field (energy, fluence, particle type). The defaults are set to correspond to the column names used in \code{\link{AT.SPC.read}}. The user can change the list refering to other columns (of same length), e.g. when having multiplied the fluence per primary in the spc to get a realistic dose.}
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

# Load required libraries
require(libamtrack)
require(lattice)

# Use example data set
file.name 			   <- system.file("extdata", "libamtrack.12C.H2O.active3.MeV27000.spc", package = "libamtrack")
endian                 <- c("big", "little")[2]

# Read in spc data, we use old style R reader, so endianess has to be given
spc                    <- AT.SPC.read ( file.name = file.name,
                                        endian    = endian,
                                        flavour   = "vanilla")$spc
										
# Compute and plot dose with depth from spc data 
df  <- AT.SPC.tapply( spc           = spc, 
                      INDEX         = "depth.g.cm2", 
                      FUN           = AT.total.D.Gy, 
                      names.results = "D.Gy")

xyplot( D.Gy ~ depth.g.cm2,
        df,
		type = "o",
		ylab = "dose per primary / Gy"
		xlab = depth / (g/cm2)")

# Compute the dose with depth, but differentiate contribution
# from individual particle species 
AT.SPC.tapply( spc           = spc, 
               INDEX         = c("depth.g.cm2", "particle.no"), 
               FUN           = AT.total.D.Gy, 
               names.results = "D.Gy")

# Compute the fluence with depth, but differentiate contribution
# from individual particle species 
# As 'mean' does not take any arguments obeying the mixed field naming
# convention in libamtrack, but a vector "x",
# one must specify this as additional argument
# Also, for the fluence data, it is referred to a non-standard column (not fluence.cm2 = "N.per.primary"),
# which contains the fluence scaled to 


spc$fluence.cm2 <- spc$N.per.primary

AT.SPC.tapply( spc                   = spc, 
               INDEX                 = c("depth.g.cm2", "particle.no"), 
               FUN                   = mean,
			   mixed.field.arguments = list(E.MeV.u = "E.mid.MeV.u", fluence.cm2 = "fluence.cm2", particle.no = "particle.no")
               additional.arguments  = list(c("x", "fluence.cm2")),
               names.results         = "fluence.cm2")
}
}
