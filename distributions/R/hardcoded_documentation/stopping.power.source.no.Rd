% Documentation for 'stopping.power.source.no' variable
\name{stopping.power.source.no}

\alias{stopping.power.source.no}
\alias{Sources of stopping power data}
\alias{Source index of stopping power data}

\title{Sources of stopping power data}

\description{Three sources of stopping power data are implemented in libamtrack. They are referred to using keywords (see below). In addition, the user can provide a filename point to custom data in text format. The built-in sources provide electronic mass stopping power as a function of specific energy (details see below).\cr
Scaling with density (or use of material parameters by Bethe formular) is done using the material specified. In case of custom stopping power data, the user is responsible to match the correct material, e.g. by use of a custom material (\code{\link{material.no}=0}). \cr


Usually not relevant, but for high-accuracy application, e/amu, I
\cr
Stopping power data handling underwent a major rework in v0.6.0.}

\seealso{More information, especially the explicit formulas used, valid energy ranges and references to literature are found in the libamtrack reference manual (\url{http://libamtrack.dkfz.org/libamtrack/images/3/31/LibamtrackReferenceManual.pdf}).}

\details{
	\bold{The stopping power sources are:}
	\tabular{rl}{
		\bold{\code{stopping.power.source}} \tab \bold{description}\cr
  	Bethe \tab Analytical Bethe formular, including density effect, for any material, if composition is given.\cr
		PSTAR \tab Tabulated data from NIST's PSTAR code, for other ions scales by the square of the effective charge. !Scaling by E/n\cr
		ICRU \tab Tabulated data from ICRU49 (H, He) and ICRU73, for liquid water only. N.B. the I value inconsistency between the two reports
    libdEdx_SOURCE \tab Reserved keyword for later use with libdEdx and data source SOURCE
	}
}

\keyword{Sources of stopping power data}
