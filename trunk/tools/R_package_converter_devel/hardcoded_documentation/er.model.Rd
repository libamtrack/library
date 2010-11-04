% Documentation for 'er.model' variable
\name{er.model}

\alias{er.model}
\alias{electron-range models}
\alias{electron-range model index}
\alias{electron-range model numbering}

\title{Electron-range models}

\description{In libamtrack, there are currently seven parametrizations of the maximal distance electrons travel from the center of a particle track as a function of the primary particle's energy. They are referred to using sequential positive integer numbers (see details).}

\seealso{More information on the electron-range models, especially the explicit formulas used, valid energy ranges and references to literature are found in the libamtrack reference manual (\url{http://libamtrack.dkfz.org/libamtrack/images/3/31/LibamtrackReferenceManual.pdf}).}

\details{
	\tabular{rl}{
		\bold{\code{er.model}} \tab \bold{description}\cr
		1 \tab A simple test function (constant range)\cr
		2 \tab Parametrization proposed by Butts and Katz (Butts & Katz, 1972)\cr
		3 \tab Parametrization proposed by M. Waligorski (Waligorski, 1985)\cr
		4 \tab Parametrization proposed by O. Geiss (Geiss, 1997)\cr
		5 \tab Parametrization proposed by M. Scholz (Scholz, 1997)\cr
		6 \tab Parametrization proposed by J. Edmund (Edmund et al., 2007)\cr
		7 \tab Parametrization proposed by Tabata (Tabata, 1972)
	}
}

\keyword{electron-range models}
