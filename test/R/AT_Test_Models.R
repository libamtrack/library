################################################################################################
# R test script for implemented models
################################################################################################
# Copyright 2006, 2009 Steffen Greilich / the libamtrack team
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
dyn.load("../../Release/libamtrack.dll")

# load wrapping scripts
source("../../wrapper/R/AmTrack.R")

# necessary library for plotting
library("lattice")

####################### INACTIVATION PROBABILITY function test #########################################

# Mixed field of 3 particles:

E.MeV.u <- c( 1, 10, 100)
particle.no <- c( 1001, 6012, 1001)
fluence.cm2 <- c( -1, -5, -4)  # in Gy

factor <- 10^seq( -4, -1, length = 5)

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

df1$SPIFF.ion	<- numeric(nrow(df1))
df1$SPIFF.gamma	<- numeric(nrow(df1))
df1$SPIFF.check	<- numeric(nrow(df1))

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
 	    
 	    current.fluence.cm2 <- df1$factor[kk] * fluence.cm2
 	    
 	    cat( "Fluences ", current.fluence.cm2 , "\n")
 	        
 	    SPIFF.res <- AT.run.SPIFF.method( E.MeV.u = E.MeV.u,
									particle.no = particle.no,
									fluence.cm2 = current.fluence.cm2,
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
         cat( "SPIFF results ", SPIFF.res , "\n")
         df1$SPIFF.ion[kk] <- SPIFF.res[2]
         df1$SPIFF.gamma[kk] <- SPIFF.res[3]
         df1$SPIFF.check[kk] <- SPIFF.res[1]
         
         GSM.res <- AT.run.GSM.method( E.MeV.u = E.MeV.u,
                                    particle.no = particle.no,
                                    fluence.cm2 = current.fluence.cm2,
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
         df1$GSM.ion[kk] <- GSM.res[2]
         df1$GSM.gamma[kk] <- GSM.res[3]
         df1$GSM.check[kk] <- GSM.res[1]
			}
		}
}

df1

# plots...

p1 <- xyplot( SPIFF.ion ~ factor | ER.model.name , groups = RDD.model.name, ref = TRUE, data=df1, pch = ".", lty = 1, type = "l", xlab = "Distance [m]", ylab = "Inactivation prob.", auto.key = list(title = "Protons in liquid water",points = FALSE, lines = TRUE), scales = list(log = 10))
p2 <- xyplot( GSM.ion ~ factor | ER.model.name , groups = RDD.model.name, ref = TRUE, data=df1, pch = ".", lty = 1, type = "l", xlab = "Distance [m]", ylab = "Inactivation prob.", auto.key = list(title = "Protons in liquid water",points = FALSE, lines = TRUE), scales = list(log = 10))


pdf("Models.pdf")

p1
p2

dev.off()
