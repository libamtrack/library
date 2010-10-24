% Documentation for 'gamma.model' and 'gamma.parameters' variable
\name{gamma.model}

\alias{gamma.model}
\alias{gamma response models}
\alias{X ray response models}
\alias{gamma response model index}
\alias{gamma response model numbering}

\title{Gamma (X ray) response models}

\description{In libamtrack, there are currently six parametrizations of the gamma / X ray response of a system. They are referred to using sequential positive integer numbers (see details).}

\seealso{More information on models, especially the explicit formulas used, valid energy ranges and references to literature are found in the libamtrack reference manual (\url{http://libamtrack.dkfz.org/libamtrack/images/3/31/LibamtrackReferenceManual.pdf}).}

\details{
	The implemented gamma models are:
	\tabular{rl}{
		\bold{\code{gamma.model}} \tab \bold{description}\cr
		1 \tab A simple test function (constant over entire track width/electron range)\cr
		Rest \tab TODO \cr
	}

	The corresponding gamma parameters and their default value are:
	\tabular{rlll}{
		\bold{\code{gamma.model}} \tab \bold{\code{gamma.parameters[1]}} \tab \bold{\code{gamma.parameters[2]}} \tab \bold{\code{gamma.parameters[3]}}\cr
		1 \tab NA \tab NA \tab NA \cr
	}
}

\keyword{gamma response models}
