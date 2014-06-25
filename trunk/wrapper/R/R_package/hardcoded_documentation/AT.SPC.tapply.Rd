\name{AT.SPC.tapply}
\alias{AT.SPC.tapply}
\title{AT.SPC.tapply}
\description{Similar to R's \code{tapply} this applies a function to cells defined by indices in a spc data object}
\usage{AT.SPC.tapply( spc, INDEX, FUN, mixed.field.arguments = list(E.MeV.u = "E.mid.MeV.u", fluence.cm2 = "N.per.primary", particle.no = "particle.no"), additional.arguments = NULL, names.results = NULL)
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
\dontrun{
# Load required libraries
require(libamtrack)
require(lattice)

# Use example data set
file.name 	       <- system.file("extdata", "libamtrack.12C.H2O.active3.MeV27000.zip", package = "libamtrack")
endian                 <- c("big", "little")[2]

# Read in spc data, we use old style R reader, so endianess has to be given
spc                    <- AT.SPC.read ( file.name = file.name,
                                        endian    = endian,
                                        flavour   = "vanilla")$spc

# Translate particle numbers in particle names (looks better)
spc$particle.name      <- AT.particle.name.from.particle.no(particle.no = spc$particle.no)
spc$particle.name      <- factor(spc$particle.name)
										
# Compute and plot dose per primary with depth from spc data 
df1 <- AT.SPC.tapply( spc                   = spc, 
                      INDEX                 = "depth.g.cm2", 
                      FUN                   = AT.total.D.Gy, 
                      additional.arguments  = list( c("material.no", "1", FALSE),                   # Water
                                                    c("stopping.power.source.no", "0", FALSE)),     # PSTAR
                      names.results         = "D.Gy")

xyplot( D.Gy ~ depth.g.cm2,
        df1,
        type = "o",
		grid = TRUE,
        ylab = "dose per primary / Gy",
        xlab = "depth / (g/cm2)")

# Compute and plot dose (1 Gy entrance) with depth from spc data 
# To do so, first the fluence per primary has to be scaled and then
# AT.SPC.tapply has to be referred to the new non-standard column (default: fluence.cm2 = "N.per.primary")

fluence.factor  <- 1.0 / df1$D.Gy[1]					# factor between dose per primary (from above) and 2 Gy
spc$fluence.cm2 <- spc$N.per.primary * fluence.factor

df2 <- AT.SPC.tapply( spc                   = spc, 
                      INDEX                 = "depth.g.cm2", 
                      FUN                   = AT.total.D.Gy,
					  mixed.field.arguments = list( fluence.cm2 = "fluence.cm2", E.MeV.u = "E.mid.MeV.u", particle.no = "particle.no"),
                      additional.arguments  = list( c("material.no", "1", FALSE),                   # Water
                                                    c("stopping.power.source.no", "0", FALSE)),     # PSTAR
                      names.results         = "D.Gy")

xyplot( D.Gy ~ depth.g.cm2,
        df2,
        type = "o",
		col  = "red",
        grid = TRUE,
		ylab = "dose / Gy",
        xlab = "depth / (g/cm2)")
			   
			   
# Compute and plot the dose with depth, but differentiate contribution
# from individual particle species
# Also: use nicer (human-readable) code for material
df3  <- AT.SPC.tapply( spc                  = spc, 
                       INDEX                = c("depth.g.cm2", "particle.name"),  
                       FUN                  = AT.total.D.Gy, 
                       additional.arguments = list(c("material.no", "AT.material.no.from.material.name('Water, Liquid')", FALSE),
                                                   c("stopping.power.source.no", "0", FALSE)),
                       names.results        = "D.Gy")

xyplot( log10(D.Gy) ~ depth.g.cm2,
        df3,
        groups      = particle.name,
        type     	= 'o',
        grid        = TRUE,
        auto.key    = list(space = 'right'),
        ylab        = "log10( dose per primary / Gy )",
        xlab        = 'depth / (g/cm2)')
}
}
