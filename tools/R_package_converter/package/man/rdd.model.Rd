% Documentation for 'rdd.model' and 'rdd.parameters' variable
\name{rdd.model}

\alias{rdd.model}
\alias{radial dose distribution models}
\alias{radial dose distribution model index}
\alias{radial dose distribution model numbering}

\title{Radial dose distribution (rdd) models}

\description{In libamtrack, there are currently seven parametrizations of distribution of local dose around particle track as a function of the primary particle's energy and type. They are referred to using sequential positive integer numbers (see details).}

\seealso{More information on the radial dose distribution models, especially the explicit formulas used, valid energy ranges and references to literature are found in the libamtrack reference manual (\url{http://libamtrack.dkfz.org/libamtrack/images/3/31/LibamtrackReferenceManual.pdf}).}

\details{
	\bold{The implemented rdd models are:}
	\tabular{rl}{
		\bold{\code{rdd.model}} \tab \bold{description}\cr
		1 \tab A simple test function (constant over entire track width/electron range)\cr
		2 \tab Point target rdd proposed by Butts and Katz (Butts & Katz, 1969)\cr
		3 \tab Simplified extended target rdd proposed by O. Geiss (Geiss, 1997)\cr
		4 \tab Simplified extended target rdd proposed by Katz (Katz et al., 1972)\cr
		5 \tab Point target rdd proposed by Cucinotta (Cucinotta, XXXX)\cr
		6 \tab Extended target rdd proposed by Katz (Katz, XXXX)\cr
		7 \tab Extended target rdd proposed by Cucinotta (Cucinotta, XXXX)
	}

	\bold{The corresponding rdd parameters are:}
	\tabular{rlll}{
		\bold{\code{rdd.model}} \tab \bold{\code{rdd.parameters[1]}} \tab \bold{\code{rdd.parameters[2]}} \tab \bold{\code{rdd.parameters[3]}}\cr
		1 \tab NA \tab NA \tab NA \cr
		2 \tab minimal radius (integration limit) [m]\tab lower dose cut-off [Gy]\tab NA\cr
		3 \tab core diameter/target size [m]\tab NA\tab NA\cr
		4 \tab core diameter/target size [m]\tab lower dose cut-off [Gy]\tab NA\cr
		5 \tab minimal radius (integration limit) [m]\tab lower dose cut-off [Gy]\tab NA\cr
		6 \tab minimal radius (integration limit) [m]\tab core diameter/target size [m]\tab lower dose cut-off [Gy]\cr
		7 \tab minimal radius (integration limit) [m]\tab core diameter/target size [m]\tab lower dose cut-off [Gy]
	}

	\bold{The corresponding default values are:}
	\tabular{rlll}{
		\bold{\code{rdd.model}} \tab \bold{\code{rdd.parameters[1]}} \tab \bold{\code{rdd.parameters[2]}} \tab \bold{\code{rdd.parameters[3]}}\cr
		1 \tab NA \tab NA \tab NA \cr
		2 \tab 1e-10\tab 1e-10\tab NA\cr
		3 \tab 5e-8\tab NA\tab NA\cr
		4 \tab 5e-8\tab 1e-10\tab NA\cr
		5 \tab 5e-11\tab 1e-10\tab NA\cr
		6 \tab 1e-10\tab 1e-8\tab 1e-10\cr
		7 \tab 5e-11\tab 1e-8\tab 1e-10
	}
}

\keyword{radial dose distribution models}
