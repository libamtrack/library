################################################################################################
# R test script for GSM model
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

# clear workspace
rm( list = ls() )

# Build latest version of libamtrack and load for direct access
source("AT_Test_PreRun.R")

require(lattice)

####################### GSM algorithm function test #########################################

# Mixed field of 3 particles :
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

df1$GSM.ion   <- numeric(nrow(df1))
df1$GSM.gamma <- numeric(nrow(df1))
df1$GSM.check <- numeric(nrow(df1))

# 1, D0, c, m, 0
GR.model      <- 2
GR.parameters <- c(1, 3.0, 1, 2, 0)

N.runs        <- 10
nX            <- 5
voxel.size.m  <- 1e-7

for( i in ER.model) {
 # i <- ER.model[1]
   ii 			<- df1$ER.model == i
	for( j in RDD.model) {
      # j <- RDD.model[1]
 	jj 			<- ((df1$RDD.model == j) & (df1$ER.model == i))
 	 for( k in factor ){   
          # k <- factor[1]
 	    kk 			<- ((df1$RDD.model == j) & (df1$ER.model == i) & (df1$factor == k))
 	    
 	    current.fluence.cm2.or.dose.Gy <- df1$factor[kk] * fluence.cm2.or.dose.Gy
 	    
 	    cat( "Fluences ", current.fluence.cm2.or.dose.Gy , "\n")
         
        GSM.res <- AT.run.GSM.method( E.MeV.u = E.MeV.u,
                                    particle.no = particle.no,
                                    fluence.cm2.or.dose.Gy = current.fluence.cm2.or.dose.Gy,
                                    material.no = material.no,
                                    rdd.model = j,
                                    rdd.parameters = RDD.parameters[[j]],
                                    er.model = i,
                                    gamma.model = GR.model,
                                    gamma.parameters = GR.parameters,
                                    N.runs = N.runs,
                                    write.output = F,
                                    nX = nX,
                                    voxel.size.m = voxel.size.m,
                                    lethal.events.mode = F)
         cat( "############\nGSM results\n############\n")
         print(GSM.res)
         df1$GSM.ion[kk] <- GSM.res$S.HCP
         df1$GSM.gamma[kk] <- GSM.res$S.gamma
         df1$GSM.check[kk] <- GSM.res$d.check
			}
		}
}

df1

# plots...

p1 <- xyplot( GSM.ion ~ factor | ER.model.name , groups = RDD.model.name, ref = TRUE, data=df1, pch = ".", lty = 1, type = "l", xlab = "Dose factor", ylab = "Response", auto.key = list(title = "ion response",points = FALSE, lines = TRUE), scales = list(log = 10))
p2 <- xyplot( GSM.gamma ~ factor | ER.model.name , groups = RDD.model.name, ref = TRUE, data=df1, pch = ".", lty = 1, type = "l", xlab = "Dose factor", ylab = "Response", auto.key = list(title = "gamma response",points = FALSE, lines = TRUE), scales = list(log = 10))
p3 <- xyplot( GSM.check ~ factor | ER.model.name , groups = RDD.model.name, ref = TRUE, data=df1, pch = ".", lty = 1, type = "l", xlab = "Dose factor", ylab = "Dose check", auto.key = list(title = "CPPSC dose check",points = FALSE, lines = TRUE), scales = list(log = 10))


pdf("GSM.pdf")

p1
p2
p3

dev.off()
