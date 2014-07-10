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
	\bold{The implemented gamma models are:}
	\tabular{rl}{
		\bold{\code{gamma.model}} \tab \bold{description}\cr
		1 \tab A simple test function (constant over entire track width/electron range)\cr
		2 \tab Generalized hit/target model\cr
		3 \tab Radioluminescence response\cr
		4 \tab Simple exponential saturation (special case of 2: one hit / one target)\cr
		5 \tab Linear-quadratic model\cr
		6 \tab Mejdahl's TLD response for LiF as used by Geiss et al. (1998)		
	}

	\bold{The corresponding gamma parameters are:}
	\tabular{rlllll}{
		\bold{\code{gamma.model}} \tab \bold{\code{gamma.parameters[1]}} \tab \bold{\code{gamma.parameters[2]}} \tab \bold{\code{gamma.parameters[3]}} \tab \bold{\code{gamma.parameters[4]}} \tab \bold{\code{gamma.parameters[5]}}\cr
		1 \tab NA \tab NA \tab NA \tab NA \tab NA\cr
		2 \tab rel. response\tab characteristic dose [Gy]\tab hittedness \tab multiplicity \tab 0 to finish or another set of 4 parameters [1]-[4]\cr
		3 \tab max. response\tab transition dose [Gy]\tab dynamic \tab NA \tab NA\cr
		4 \tab max. response\tab characteristic dose [Gy] \tab NA \tab NA \tab NA\cr
		5 \tab alpha [1/Gy]\tab beta[1/Gy2]\tab transition dose [Gy]\tab NA \tab NA\cr
		6 \tab const \tab k1 \tab a1 \tab k2 \tab a2\
	}

	\bold{The corresponding default/example values are:}
	\tabular{rlllll}{
		\bold{\code{gamma.model}} \tab \bold{\code{gamma.parameters[1]}} \tab \bold{\code{gamma.parameters[2]}} \tab \bold{\code{gamma.parameters[3]}} \tab \bold{\code{gamma.parameters[4]}} \tab \bold{\code{gamma.parameters[5]}}\cr
		1 \tab NA \tab NA \tab NA \tab NA \tab NA\cr
		2 \tab 1.0\tab 10.0\tab 1 \tab 1 \tab 0\cr
		3 \tab 1e6\tab 20.0\tab 6.0 \tab NA \tab NA\cr
		4 \tab 1e6\tab 10.0 \tab NA \tab NA \tab NA\cr
		5 \tab 0.2\tab 0.018\tab 30.0\tab NA \tab NA\cr
		6 \tab 6e3 \tab 0.59 \tab 5e-4 \tab 0.41 \tab 2e-6
	}	
}

\keyword{gamma response models}
