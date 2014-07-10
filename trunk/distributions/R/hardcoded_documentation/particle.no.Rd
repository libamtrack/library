% Documentation for 'particle.no' variable
\name{particle.no}

\alias{particle.no}
\alias{particle number}
\alias{particle index}
\alias{particle numbering}

\title{Particle indexing convention}

\description{The particle indexing system in libamtrack works as follows: for a particle with atomic number Z and mass number A the index is \code{particle.no = 1000 * Z + A}. This ensures uniqueness for all particle types used in libamtrack. Also, an extension to anti-protons or ions could be made using negative numbers.}

\seealso{You can use \code{AT.particle.name.from.particle.no} and \code{AT.particle.no.from.particle.name} to convert particle names into particle indices and vice versa.}

\details{
	\itemize{
		\item{}{For a proton (Z=1,A=1) \code{particle.no} is 1001.}
		\item{}{For deuterium (Z=1, A=2) \code{particle.no} is 1002.}
		\item{}{For He-4 (Z=2, A=4) \code{particle.no} is 2004.}
		\item{}{For C-12 (Z=6, A=12) \code{particle.no} is 6012.}
		\item{}{For a U-238 (Z=92, A=238) \code{particle.no} is 92238.}
	}
}

\keyword{particle indexing}
