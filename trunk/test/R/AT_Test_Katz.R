################################################################################################
# R test script for implemented electron range models
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
dyn.load("../../Release/libAmTrack.dll")

# load wrapping scripts
source("../../wrapper/R/AmTrack.R")

# necessary library for plotting
library("lattice")

####################### DOSE(DISTANCE) function test #########################################

# radius definitions:
r.m <- 10^seq (-14, -3, length.out=100)

# electron range models definition
ER.model.names <- c("simple test",  "Butts & Katz' (linear)",  "Waligorski's (power-law wmax)",  "Geiss' (power-law E)", "Scholz' (power-law E)", "Edmund' (power-law wmax)","Tabata")
ER.model <- c(2,3,4,5,6,7)

# RDD models definition
RDD.model.names <- c("Simple step test function",  "Katz' point target", "Geiss'", "Site", "Cucinotta", "KatzExtTarget", "CucinottaExtTarget")
RDD.model <- c(1,2,3,4,5,6,7)

# RDD parameters
RDD.parameters <- list(c(1),c(1e-10,1e-10), c(1e-8),c(1e-8,1e-10),c(5e-11,1e-10),c(1e-10,1e-8,1e-10),c(5e-11,1e-8,1e-10))

# data frame setup
df1 <- expand.grid( r.m = r.m, ER.model = ER.model, RDD.model = RDD.model )

df1$ER.model.name	<- as.character(ER.model.names[df1$ER.model])
df1$RDD.model.name	<- as.character(RDD.model.names[df1$RDD.model])

df1$inact.prob	<- numeric(nrow(df1))

df1$r.m.check	<- numeric(nrow(df1))


material.no  <-  c(1)
E.MeV.u      <-  c(100.0)
particle.no  <-  c(1001)

# 1, D0, c, m, 0
GR.parameters <- c(1, 3.0, 1, 2, 0)

for( i in ER.model) {
 ii 			<- df1$ER.model == i
	for( j in RDD.model) {
 	jj 			<- ((df1$RDD.model == j) & (df1$ER.model == i)) 
			df1$inact.prob[jj] <- AT.Katz.inactivation.probability(r.m = df1$r.m[jj], E.MeV.u, particle.no, material.no, rdd.model = j, rdd.parameters = RDD.parameters[[j]], er.model = i, GR.parameters)
		}
}

#df1

# plots...

p1 <- xyplot( inact.prob ~ r.m | ER.model.name , groups = RDD.model.name, ref = TRUE, data=df1, pch = ".", lty = 1, type = "l", xlab = "Distance [m]", ylab = "Inactivation prob.", auto.key = list(title = "Protons in liquid water",points = FALSE, lines = TRUE), scales = list(log = 10))
p2 <- xyplot( inact.prob ~ r.m | RDD.model.name, groups = ER.model.name, ref = TRUE, data=df1, pch = ".", lty = 1, type = "l", xlab = "Distance [m]", ylab = "Inactivation prob.", auto.key = list(title = "Protons in liquid water",points = FALSE, lines = TRUE), scales = list(log = 10))

p1logx <- xyplot( inact.prob ~ r.m | ER.model.name , groups = RDD.model.name, ref = TRUE, data=df1, pch = ".", lty = 1, type = "l", xlab = "Distance [m]", ylab = "Inactivation prob.", auto.key = list(title = "Protons in liquid water",points = FALSE, lines = TRUE), scales = list( x = list(log = 10)))
p2logx <- xyplot( inact.prob ~ r.m | RDD.model.name, groups = ER.model.name, ref = TRUE, data=df1, pch = ".", lty = 1, type = "l", xlab = "Distance [m]", ylab = "Inactivation prob.", auto.key = list(title = "Protons in liquid water",points = FALSE, lines = TRUE), scales = list( x = list(log = 10)))


pdf("Katz.pdf")

p1
p2
p1logx
p2logx

dev.off()
