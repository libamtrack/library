% Documentation for 'stopping.power.source.no' variable
\name{stopping.power.source.no}

\alias{stopping.power.source.no}
\alias{Sources of stopping power data}
\alias{Source index of stopping power data}

\title{Sources of stopping power data}

\description{In libamtrack, there are currently four stopping power data sources. They are referred to using sequential positive integer numbers (see details). Please note, that not all sources can be used for every material. The most versatile, but also inaccurate source in the context is the Bethe formular.}

\seealso{More information, especially the explicit formulas used, valid energy ranges and references to literature are found in the libamtrack reference manual (\url{http://libamtrack.dkfz.org/libamtrack/images/3/31/LibamtrackReferenceManual.pdf}).}

\details{
	\bold{The stopping power sources are:}
	\tabular{rl}{
		\bold{\code{stopping.power.source.no}} \tab \bold{description}\cr
		0 \tab Tabulated data from NIST's PSTAR code, for other ions scales by the square of the effective charge (default until v0.5.0r959, 'Black Wombat')\cr
		1 \tab Analytical Bethe formular, including density effect\cr
		2 \tab Tabulated data from the SHIELD-HIT code, generate with an extended Bethe-Bloch approach, for liquid water only\cr
		3 \tab Tabulated data from ICRU49 (H, He) and ICRU73, for liquid water only. N.B. the I value inconsistency between the two reports
	}
}

\keyword{Sources of stopping power data}
