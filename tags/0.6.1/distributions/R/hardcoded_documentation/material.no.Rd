% Documentation for 'material.no' variable
\name{material.no}

\alias{material.no}
\alias{material number}
\alias{material index}
\alias{material numbering}

\title{Material indexing convention}

\description{Currently seven predefined materials exist in libamtrack. They are referred to using sequential positive integer numbers (see details). The material definition follows largely the pstar and astar definitions of NIST (http://physics.nist.gov/cgi-bin/Star/compos.pl?ap).

\code{material.no = 0} points in contrast to a material defined by the user during run-time. This feature is still under development and not available yet.}

\seealso{You can use \code{AT.material.name.from.material.no} and \code{AT.material.no.from.material.name} to convert material names into material indices and vice versa. \code{AT.get.materials.data} will return all available data for a predefined material.}

\details{
	The predefined materials in libamtrack are:
	\tabular{rl}{
		\bold{\code{material.no}} \tab \bold{material}\cr
		1 \tab Water, Liquid \cr
		2 \tab Aluminum Oxide (Al2O3) \cr
		3 \tab Aluminum \cr
		4 \tab PMMA (Polymethylmethacrylate, Perspex, Plexiglass, [C5H8O2]n) \cr
		5 \tab Alanine \cr
		6 \tab Lithium Fluoride (LiF) \cr
		7 \tab Air (dry, at sea level) \cr
		8 \tab Silicon \cr
		9 \tab Copper
	}
}

\keyword{material indexing}
