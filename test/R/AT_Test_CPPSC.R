################################################################################################
# R test script for CPPSC
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

####################################### THIS CAN BE LATER COLLECTED IN A SCRIPT SHARED BY ALL R TESTSCRIPTS ##########
# Trigger build of latest libamtrack version for direct R access 
# and load resulting dll / wrappers
start.dir           <- getwd()
if(try(setwd("..\\..\\wrapper\\R\\R_direct_access")) == FALSE){stop("Please start script from /test")}
try(system("create.direct.access.windows.bat"))

# load library and wrappers
try(dyn.load("libamtrack.dll"))
try(dyn.load("libamtrack.so"))
try(dyn.load("libamtrack.dylib"))

source("libamtrack.R")

setwd(cur.dir)

# necessary library for plotting
require("lattice")
####################################### THIS CAN BE LATER COLLECTED IN A SCRIPT SHARED BY ALL R TESTSCRIPTS ##########


####################### CPPSC function test #########################################

# Mixed field of 3 particles:

E.MeV.u                <- c( 1, 10, 100)
particle.no            <- c( 1001, 6012, 1001)
fluence.cm2.or.dose.Gy <- c( -1, -5, -4)  # in Gy
material.no            <- 1 # Liquid water

CPPSC.res              <- .C("AT_run_CPPSC_method",   E.MeV.u = E.MeV.u,
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
