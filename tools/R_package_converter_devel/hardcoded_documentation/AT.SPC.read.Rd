\name{AT.SPC.read}
\alias{AT.SPC.read}
\title{AT.SPC.read}
\description{Read a spc-formatted data file with energy-fluence in n depth steps. For original TRiP format definition by M. Kraemer, please see \url{http://bio.gsi.de/DOCS/TRiP98BEAM/DOCS/trip98fmtspc.html}.}
\usage{AT.SPC.read( file.name, endian, mean, raw, compress)
}
\arguments{
  \item{file.name}{ name of spc-file to be read (with extension).}
  \item{endian}{ byte-order in spc file, "big" or "little" depending on the machine of origin.}
  \item{mean}{ 0 if geometrical bin-mid should be used, 1 for arithmetic (default).}
  \item{raw}{ if TRUE all data from spc file will be kept, if FALSE (default) the data frame will be reformatted for use with libamtrack, including using \code{\link{particle.no}} and \code{\link{material.name}} conventions.}
  \item{compress}{ if TRUE all entries with zero fluence will be removed to save memory.}
}
\value{
  A list with the following elements:
  \item{file.type}{ Type of spc file}
  \item{file.version}{ Version}
  \item{file.date}{ Date of version}
  \item{target.name}{ Target material, transformed to \code{\link{material.name}} if \code{raw = T})}
  \item{projectile.name}{ Projectile, according \code{\link{particle.name}})}
  \item{beam.energy.MeV.u}{ Initial projectile energy}
  \item{peak.position.g.cm2}{ Position (mass density) of Bragg peak}
  \item{normalization}{ Normalization constant (not used)}
  \item{n.depth.steps}{ Number of depth steps in data set}
  \item{data}{ Actual data, organized in a data-frame}
}
\examples{
# Read in TRiP data from 142 MeV protons and
# reformat/compress for use with libamtrack
# functions
# 
# spc <- AT.SPC.read( file.name = "TRIP_GSI_1H.H2O.MeV14260.spc",
#                     endian    = "little",
#                     mean      = 1,
#                     raw       = FALSE,
#                     compress  = TRUE)
#
#
# Read in same data but give complete overview on data content
# 
# spc.raw <- AT.SPC.read( file.name = "TRIP_GSI_1H.H2O.MeV14260.spc",
#                         endian    = "little",
#                         mean      = 1,
#                         raw       = TRUE,
#                         compress  = FALSE)
}
