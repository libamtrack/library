################################################################################################
# Testing script for AT_SuccessiveConvolutions workflow
################################################################################################
# Copyright 2006, 2010 The libamtrack team
# 
# This file is part of the AmTrack program (libamtrack.sourceforge.net).
#
# AmTrack is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# AmTrack is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# long with AmTrack (file: copying.txt).
# If not, see <http://www.gnu.org/licenses/>
#
################################################################################################
rm(list = ls())

# Build latest version of libamtrack and load for direct access
recompile <- TRUE
source("AT_Test_PreRun.R")

library(lattice)

#################################################
# Set parameters

E.MeV.u <- 50
particle.no <- 6012
dose.Gy <- 10
fluence.cm2.or.dose.Gy <- -dose.Gy
material.no <- 1				# Liquid water
RDD.model <- 3				# Geiss RDD
RDD.parameters <- 5e-8			# a0 = 50 nm
ER.model <- 3				# Geiss ER
#gamma.model <- 5				# General hit-target
#gamma.parameters <- c(0.2,0.02,10,0,0)	# One single-hit-single-target (exp-sat) component, characteristic dose 10 Gy
gamma.model            <- 2                                                    # General hit/target
gamma.parameters       <- c(1,10,1,1,0)                                        # Exp.-sat. (one-hit/one-target) with 10 Gy sat.-dose
N2 <- 10				# 20 bins per factor 2 in histograms
fluence.factor <- 1				# use fluence as given
write.output <- F				# no log file
shrink.tails <- T				# cut insignificant tails
shrink.tails.under <- 1e-30			# cut them in case contribution to first moment is lower than
adjust.N2 <- T				# adjust bin width during convolution
lethal.events.mode <- F				# use survival instead of activation
stopping.power.source.no   <- 0    # PSTAR

print(dose.Gy)

fluence.cm2 <- AT.fluence.cm2.from.dose.Gy(
			E.MeV.u,
			particle.no,
			dose.Gy,
			material.no,
			stopping.power.source.no)

print(fluence.cm2)

# Get histogram size for single-impact dose distribution
res.get.f1.array.size	<-	AT.n.bins.for.single.impact.local.dose.distrib(E.MeV.u = E.MeV.u,
								particle.no = particle.no,
                                material.no = material.no,
								rdd.model = RDD.model,
								rdd.parameter = RDD.parameters,
								er.model = ER.model,
								N2 = N2,
								stopping.power.source.no = stopping.power.source.no)[[1]]
								
print(res.get.f1.array.size)

# Get single-impact RDD parameters
res.f1.parameters	<-	AT.RDD.f1.parameters.mixed.field(E.MeV.u = E.MeV.u,
								particle.no = particle.no,
								material.no = material.no,
								rdd.model = RDD.model,
								rdd.parameter = RDD.parameters,
								er.model = ER.model,
								stopping.power.source.no = stopping.power.source.no)[[1]]
print(res.f1.parameters)

# Get single-impact dose distribution
res.get.f1	<-	AT.single.impact.local.dose.distrib(	E.MeV.u = E.MeV.u,
							particle.no = particle.no,
							fluence.cm2.or.dose.Gy = fluence.cm2.or.dose.Gy,
							material.no = material.no,
							rdd.model = RDD.model,
							rdd.parameter = RDD.parameters,
							er.model = ER.model,
							N2 = N2,
							n.bins.f1 = res.get.f1.array.size,
							f1.parameters = res.f1.parameters,							
							stopping.power.source.no = stopping.power.source.no)

print(res.get.f1)

p1 <- plot(log10(res.get.f1$f1)~log10(res.get.f1$f1.d.Gy))


# Get mean impact number u
res.u <- AT.mean.number.of.tracks.contrib(E.MeV.u = E.MeV.u,
							particle.no = particle.no,
							fluence.cm2 = fluence.cm2,
							material.no = material.no,
							er.model = ER.model,
							stopping.power.source.no = stopping.power.source.no)

print(res.u)

# Get histogram size for low-fluence dose distribution, low-fluence impact number (u.start) and number of convolutions
res.get.f.array.size <-	AT.n.bins.for.low.fluence.local.dose.distribution(	u = res.u$returnValue,
								fluence.factor = fluence.factor,
								N2 = N2,
								f1.d.Gy = res.get.f1$f1.d.Gy,
								f1.dd.Gy = res.get.f1$f1.dd.Gy,
								f1 = res.get.f1$f1)

print(res.get.f.array.size)


# Get low-fluence dose distribution (computes again u.start -> remove from function above)
res.get.f.start	<-	AT.low.fluence.local.dose.distribution(	N2 = N2,
								 f1.d.Gy = res.get.f1$f1.d.Gy,
								 f1.dd.Gy = res.get.f1$f1.dd.Gy,
								 f1 = res.get.f1$f1,
								 n.bins.f = res.get.f.array.size$n.bins.f)

xyplot(log10(res.get.f.start$f.start) ~ log10(res.get.f.start$f.d.Gy))

# Perform SC
res.SC	<-	AT.SuccessiveConvolutions(	final.mean.number.of.tracks.contrib = res.u,
									N2 = N2,
									n.bins.f.used = res.get.f1.array.size,
									f.d.Gy = res.get.f.start$f.d.Gy,
									f.dd.Gy = res.get.f.start$f.dd.Gy,
									f = res.get.f.start$f.start,
									write.output = write.output,
									shrink.tails = shrink.tails,
									shrink.tails.under = shrink.tails.under,
									adjust.N2 = adjust.N2)

#print(res.SC$N2)
#print(res.SC$n.bins.f.used)
#print(res.SC$f0)
#print(res.SC$d)

p2 <- xyplot(log10(res.SC$f)~log10(res.SC$f.d.Gy))

pdf("SC.pdf")

p1
p2

dev.off()
