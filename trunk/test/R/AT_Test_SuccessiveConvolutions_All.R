# Testing script for AT_SuccessiveConvolutions workflow
# Created: 2010-10-06
# Creator: greilich

rm(list = ls())

try(dyn.load("..\..\lib\libamtrack.dll"))
try(dyn.load("../../lib/libamtrack.so"))
try(dyn.load("../../lib/libamtrack.dylib"))

source("../../wrapper/R/AmTrack.R")

library(lattice)

#################################################
# Set parameters

E.MeV.u <- 10
particle.no <- 1001
fluence.cm2.or.dose.Gy <- -1
material.no <- 1				# Liquid water
RDD.model <- 3				# Geiss RDD
RDD.parameters <- 5e-8			# a0 = 50 nm
ER.model <- 4				# Geiss ER
gamma.model <- 2				# General hit-target
gamma.parameters <- c(1,10,1,1,0)	# One single-hit-single-target (exp-sat) component, characteristic dose 10 Gy
N2 <- 20				# 20 bins per factor 2 in histograms
fluence.factor <- 1				# use fluence as given
write.output <- T				# no log file
shrink.tails <- T				# cut insignificant tails
shrink.tails.under <- 1e-30			# cut them in case contribution to first moment is lower than
adjust.N2 <- T				# adjust bin width during convolution
lethal.events.mode <- F				# use survival instead of activation

# Get histogram size for single-impact dose distribution
res.get.f1.array.size	<-	AT.SC.get.f1.array.size(E.MeV.u = E.MeV.u,
								particle.no = particle.no,
								fluence.cm2.or.dose.Gy = fluence.cm2.or.dose.Gy,
								material.no = material.no,
								RDD.model = RDD.model,
								RDD.parameters = RDD.parameters,
								ER.model = ER.model,
								N2 = N2)
								
print(res.get.f1.array.size$n.bins.f1)

# Get single-impact RDD parameters
res.f1.parameters	<-	AT.RDD.f1.parameters.mixed.field(E.MeV.u = E.MeV.u,
								particle.no = particle.no,
								material.no = material.no,
								ER.model = ER.model,
								RDD.model = RDD.model,
								RDD.parameters = RDD.parameters)
print(res.f1.parameters)

# Get single-impact dose distribution
res.get.f1	<-	AT.SC.get.f1(	E.MeV.u = E.MeV.u,
							particle.no = particle.no,
							fluence.cm2.or.dose.Gy = fluence.cm2.or.dose.Gy,
							material.no = material.no,
							RDD.model = RDD.model,
							RDD.parameters = RDD.parameters,
							ER.model = ER.model,
							N2 = N2,
							n.bins.f1 = res.get.f1.array.size$n.bins.f1,
							f1.parameters = res.f1.parameters)

xyplot(log10(f1)~log10(f1.d.Gy),
res.get.f1$f1)

xyplot(log10(f1*f1.dd.Gy)~log10(f1.d.Gy),
res.get.f1$f1)


# Get mean impact number u
res.u <- AT.u(E.MeV.u = E.MeV.u,
							particle.no = particle.no,
							fluence.cm2.or.dose.Gy = fluence.cm2.or.dose.Gy,
							material.no = material.no,
							er.model = ER.model)

print(res.u)

# Get histogram size for low-fluence dose distribution, low-fluence impact number (u.start) and number of convolutions
res.get.f.array.size <-	AT.SC.get.f.array.size(	u = res.u,
								fluence.factor = fluence.factor,
								N2 = N2,
								n.bins.f1 = res.get.f1.array.size$n.bins.f1,
								f1.d.Gy = res.get.f1$f1$f1.d.Gy,
								f1.dd.Gy = res.get.f1$f1$f1.dd.Gy,
								f1 = res.get.f1$f1$f1)

print(res.get.f.array.size)


# Get low-fluence dose distribution (computes again u.start -> remove from function above)
res.get.f.start	<-	AT.SC.get.f.start(	N2 = N2,
								 n.bins.f1 = res.get.f1.array.size$n.bins.f1,
								 f1.d.Gy = res.get.f1$f1$f1.d.Gy,
								 f1.dd.Gy = res.get.f1$f1$f1.dd.Gy,
								 f1 = res.get.f1$f1$f1,
								 n.bins.f = res.get.f.array.size$n.bins.f)

xyplot(log10(f.start)~log10(f.start.d.Gy),
res.get.f.start)

# Perform SC
res.SC	<-	AT.SC.SuccessiveConvolutions(	u = res.u,
									n.bins.f = res.get.f.array.size$n.bins.f,
									N2 = N2,
									n.bins.f.used = res.get.f1.array.size$n.bins.f1,
									f.d.Gy = res.get.f.start$f.start.d.Gy,
									f.dd.Gy = res.get.f.start$f.start.dd.Gy,
									f = res.get.f.start$f.start,
									write.output = write.output,
									shrink.tails = shrink.tails,
									shrink.tails.under = shrink.tails.under,
									adjust.N2 = adjust.N2)

print(res.SC$N2)
print(res.SC$n.bins.f.used)
print(res.SC$f0)
print(res.SC$d)

xyplot(log10(f)~log10(f.d.Gy),
res.SC$f)

#df <- read.table("SuccessiveConvolutions.log", header = TRUE, sep = ";")
#xyplot(xyplot(log10(f)~log10(f.d.Gy)|n.convolution,
#df))


