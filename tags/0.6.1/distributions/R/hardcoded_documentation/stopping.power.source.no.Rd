% Documentation for 'stopping.power.source' variable
\name{stopping.power.source}

\alias{stopping.power.source}
\alias{stopping.power.source.no}
\alias{Sources of stopping power data}
\alias{Source index of stopping power data}

\title{Sources of stopping power data}

\description{Three stopping power data sources are implemented in libamtrack. They are referred to using keywords or integer identifiers (see below). In addition, the user can provide a filename point to custom data in text format. The built-in sources provide electronic mass stopping power as a function of specific energy (details see below).\cr

Scaling with density (or use of material parameters by Bethe formular) is done using the material specified. In case of custom stopping power data, the user has either to make sure that the right data (i.e. electronic mass stopping power) are provided with the correct units (E/u and MeV*cm2/g) OR use any kind of data / units if they know what they are doing... (\code{\link{material.no} = 0}). \cr

Built-in data are given in specific energy per amu for the most frequent isotope (i.e. H-1, He-4 etc.) or directly in specific energy per nucleon (i.e. in case of the Bethe formula). Scaling to other isotopes happens via the number of nucleons. Thus, effects of changes in binding energy per nucleon are neglected. These are usually small (0.7 % for H-1) but might be important for high-accuracy application. The same is true for the material parameters, e.g. the I-value. In that case, it is adviced to use external tabulated data.\cr

These tabulated data have to be plain ASCII with three columns (separated by space): charge, energy, and stopping power. The datahave to be sorted in ascending order by first charge than energy. In addition, any alphanumeric comment can be inserted (in separate lines).

Stopping power data handling underwent a major rework in v0.6.0 and identifiers were changed. To avoid confusion, downward compatibility might be not maintained and integer identifiers depreceated.}

\seealso{More information, especially the explicit formulas used, valid energy ranges and references to literature are found in the libamtrack reference manual (\url{http://libamtrack.dkfz.org/libamtrack/images/3/31/LibamtrackReferenceManual.pdf}).}

\details{
	\bold{The stopping power sources are:}
	\tabular{rcl}{
		\bold{\code{stopping.power.source.no}} \tab \bold{\code{stopping.power.source}} \tab \bold{description}\cr
 		--- \tab <FILENAME> \tab data are read from external file (see details)\cr		
 		1 \tab Bethe \tab Analytical Bethe formular, including density effect, for any material, if composition is given.\cr
		2 \tab PSTAR \tab Tabulated data from NIST's PSTAR code, for other ions scales by the square of the effective charge. CAVE: Scaling by E/n!\cr
		3 \tab ICRU \tab Tabulated data from ICRU49 (H, He) and ICRU73, for liquid water only. N.B. the I value inconsistency between the two reports\cr
		--- \tab libdEdx_SOURCE \tab Reserved keyword for later use with libdEdx and data source SOURCE
	}
}

\keyword{Sources of stopping power data}

