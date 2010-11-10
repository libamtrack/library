################################################################################################
# R test script for implemented models
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

dev.off

# clear workspace
rm( list = ls() )

# load libAmTrack library
try(dyn.load("../../lib/libamtrack.dll"))
try(dyn.load("../../lib/libamtrack.so"))
try(dyn.load("../../lib/libamtrack.dylib"))

# load wrapping scripts
source("../../wrapper/R/AmTrack.R")

# necessary library for plotting
library("lattice")

####################### GSM and CPPSC function test #########################################

# Mixed field of 3 particles:

E.MeV.u <- c( 1, 10, 100)
particle.no <- c( 1001, 6012, 1001)
fluence.cm2.or.dose.Gy <- c( -1, -5, -4)  # in Gy

factor <- 10^seq( -4, -2, length = 30)

material.no <- 1

# electron range models definition
ER.model.names <- c("simple test",  "Butts & Katz' (linear)",  "Waligorski's (power-law wmax)",  "Geiss' (power-law E)", "Scholz' (power-law E)", "Edmund' (power-law wmax)","Tabata")
ER.model <- c(4)

# RDD models definition
RDD.model.names <- c("Simple step test function",  "Katz' point target", "Geiss'", "Site", "Cucinotta", "KatzExtTarget", "CucinottaExtTarget")
RDD.model <- c(3)

# RDD parameters
RDD.parameters <- list(c(1),c(1e-10,1e-10), c(5e-8),c(1e-8,1e-10),c(5e-11,1e-10),c(1e-10,1e-8,1e-10),c(5e-11,1e-8,1e-10))

# data frame setup
df1 <- expand.grid( factor = factor, ER.model = ER.model, RDD.model = RDD.model )

df1$ER.model.name	<- as.character(ER.model.names[df1$ER.model])
df1$RDD.model.name	<- as.character(RDD.model.names[df1$RDD.model])

df1$CPPSC.ion	<- numeric(nrow(df1))
df1$CPPSC.gamma	<- numeric(nrow(df1))
df1$CPPSC.check	<- numeric(nrow(df1))

df1$GSM.ion   <- numeric(nrow(df1))
df1$GSM.gamma <- numeric(nrow(df1))
df1$GSM.check <- numeric(nrow(df1))

# 1, D0, c, m, 0
GR.parameters <- c(1, 3.0, 1, 2, 0)

for( i in ER.model) {
 ii 			<- df1$ER.model == i
	for( j in RDD.model) {
 	jj 			<- ((df1$RDD.model == j) & (df1$ER.model == i))
 	 for( k in factor ){   
 	    kk 			<- ((df1$RDD.model == j) & (df1$ER.model == i) & (df1$factor == k))
 	    
 	    current.fluence.cm2.or.dose.Gy <- df1$factor[kk] * fluence.cm2.or.dose.Gy
 	    
 	    cat( "Fluences ", current.fluence.cm2.or.dose.Gy , "\n")
 	        
 	    CPPSC.res <- AT.run.CPPSC.method( E.MeV.u = E.MeV.u,
									particle.no = particle.no,
									fluence.cm2.or.dose.Gy = current.fluence.cm2.or.dose.Gy,
									material.no = material.no,
									RDD.model = j,
									RDD.parameters = RDD.parameters[[j]],
									ER.model = i,
									gamma.model = 2,
									gamma.parameters = GR.parameters,
									N2 = 3,
									fluence.factor = 1.0,
									write.output = F,
									shrink.tails = T,
									shrink.tails.under = 1e-30,
									adjust.N2 = T,
									lethal.events.mode = F)
         cat( "CPPSC results ", CPPSC.res , "\n")
         df1$CPPSC.ion[kk] <- CPPSC.res[3]
         df1$CPPSC.gamma[kk] <- CPPSC.res[4]
         df1$CPPSC.check[kk] <- CPPSC.res[2]
         
         GSM.res <- AT.run.GSM.method( E.MeV.u = E.MeV.u,
                                    particle.no = particle.no,
                                    fluence.cm2.or.dose.Gy = current.fluence.cm2.or.dose.Gy,
                                    material.no = material.no,
                                    RDD.model = j,
                                    RDD.parameters = RDD.parameters[[j]],
                                    ER.model = i,
                                    gamma.model = 2,
                                    gamma.parameters = GR.parameters,
                                    N.runs = 1,
                                    write.output = F,
                                    nX = 5,
                                    voxel.size.m = 1e-7,
                                    lethal.events.mode = F)
         cat( "GSM results ", GSM.res , "\n")
         df1$GSM.ion[kk] <- GSM.res[3]
         df1$GSM.gamma[kk] <- GSM.res[4]
         df1$GSM.check[kk] <- GSM.res[2]
			}
		}
}

df1

# plots...

p1 <- xyplot( CPPSC.ion + GSM.ion ~ factor | ER.model.name , groups = RDD.model.name, ref = TRUE, data=df1, pch = ".", lty = 1, type = "l", xlab = "Dose factor", ylab = "Response", auto.key = list(title = "ion response",points = FALSE, lines = TRUE), scales = list(log = 10))
p2 <- xyplot( CPPSC.gamma + GSM.gamma ~ factor | ER.model.name , groups = RDD.model.name, ref = TRUE, data=df1, pch = ".", lty = 1, type = "l", xlab = "Dose factor", ylab = "Response", auto.key = list(title = "gamma response",points = FALSE, lines = TRUE), scales = list(log = 10))
p3 <- xyplot( CPPSC.check + GSM.check ~ factor | ER.model.name , groups = RDD.model.name, ref = TRUE, data=df1, pch = ".", lty = 1, type = "l", xlab = "Dose factor", ylab = "Dose check", auto.key = list(title = "CPPSC dose check",points = FALSE, lines = TRUE), scales = list(log = 10))


pdf("Models.pdf")

p1
p2
p3

dev.off()
