\name{AT.gamma.response}
\alias{AT.gamma.response}
\title{AT.gamma.response}
\description{Returns a system (detector or cells) response for given doses
 according to the chosen gamma response model
}
\usage{AT.gamma.response(d.Gy, gamma.model, gamma.parameter, lethal.event.mode)
}
\arguments{
  \item{d.Gy}{ doses in Gy (array of size number.of.doses).}
  \item{gamma.model}{ gamma response model index.}
  \item{gamma.parameter}{ vector holding necessary parameters for the chose gamma response model (array of size 9).}
  \item{lethal.event.mode}{ if true computation is done in lethal event mode.}
}
\value{
% TODO proper return definition of lists!!! ADD NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)
  \item{S}{ gamma responses (array of size number.of.doses)}
}
\examples{
# Show the gamma response of two Al2O3 detectors (A & B) and two protocols ('peak' and 'total')
# parametrized in two components (single hit/single target and two hit/single target)
# as measured and published by Edmund et al., NIM B 262 (2007), 261-275
require(lattice)
d.Gy						<-	10^seq(from = log10(0.1), to = log10(25), length.out = 100)     # Compute 100 points between 0.1 and 25 Gy
gamma.model					<-	2                                                               # General hit/target model
R							<-	1                                                               # Probe A, 'peak'
Smax						<-	0.81e6
k1							<-	Smax * (R / 100)
k2							<-	Smax * (1 - R / 100)
gamma.parameter.peak.A		<-	c(	k1 = k1, D01 = 0.36, c1 = 1, m1 = 1,
								    k2 = k2, D02 = 3.06, c2 = 2, m2 = 1,
								    0)
R							<-	33                                                              # Probe A, 'total'
Smax						<-	6.2e6
k1							<-	Smax * (R / 100)
k2							<-	Smax * (1 - R / 100)
gamma.parameter.total.A	    <-	c(	k1 = k1, D01 = 1.13, c1 = 1, m1 = 1,
								    k2 = k2, D02 = 1.77, c2 = 2, m2 = 1,
								    0)
R							<-	13                                                              # Probe B, 'peak'
Smax						<-	2.84e6
k1							<-	Smax * (R / 100)
k2							<-	Smax * (1 - R / 100)
gamma.parameter.peak.B	    <-	c(	k1 = k1, D01 = 4.15, c1 = 1, m1 = 1,
								    k2 = k2, D02 = 5.14, c2 = 2, m2 = 1,
								    0)
R							<-	44                                                              # Probe B, 'total'
Smax						<-	27.6e6
k1							<-	Smax * (R / 100)
k2							<-	Smax * (1 - R / 100)
gamma.parameter.total.B	    <-	c(	k1 = k1, D01 = 2.90, c1 = 1, m1 = 1,
								    k2 = k2, D02 = 4.66, c2 = 2, m2 = 1,
								    0)
vecA <- AT.gamma.response(	d.Gy              = d.Gy,
																	gamma.model       = gamma.model,
																	gamma.parameter   = gamma.parameter.peak.A,
                                                                    lethal.event.mode = FALSE)$response
vecB <- AT.gamma.response(	d.Gy              = d.Gy,
																	gamma.model       = gamma.model,
																	gamma.parameter   = gamma.parameter.total.A,
                                                                    lethal.event.mode = FALSE)$response
vecC <- AT.gamma.response(	d.Gy              = d.Gy,
																	gamma.model       = gamma.model,
																	gamma.parameter   = gamma.parameter.peak.B,
                                                                    lethal.event.mode = FALSE)$response
vecD <- AT.gamma.response(	d.Gy              = d.Gy,
																	gamma.model       = gamma.model,
																	gamma.parameter   = gamma.parameter.total.B,
                                                                    lethal.event.mode = FALSE)$response
df			<-	data.frame(	d.Gy 		= rep(	d.Gy, 4),                                       # Compose data frame
							S 			= c(	vecA, vecB, vecC, vecD ),
								which 		= rep(	c(	rep("peak", length(d.Gy)),
														rep("total", length(d.Gy))), 2),
								probe		= c(	rep("probe A", 2 * length(d.Gy)),
													rep("probe B", 2 * length(d.Gy))))
													
xyplot(	log10(S) ~ log10(d.Gy)|probe,
            df,	groups	= which,
			type	= 'l',
			lwd		= 2,
			ylim	= log10(c(1e3, 4e7)),
			ylab	= list(	"OSL response", cex = 1.2),
			xlim	= log10(c(0.1, 25)),
			xlab	= list(	"dose / Gy", cex = 1.2),
			scales	= list(	x = list( at = log10(c(1,10,20)), labels = as.character(c(1,10,20))), y =list(	at = c(4,5,6,7), labels	= 10^(c(4,5,6,7)))),
			aspect	= 2.5)
}
