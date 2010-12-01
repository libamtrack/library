if(FALSE){
################################################################################################
# R test script for implemented RDD models
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

# load libAmTrack library
try(dyn.load("../../lib/libamtrack.dll"))
try(dyn.load("../../lib/libamtrack.so"))
try(dyn.load("../../lib/libamtrack.dylib"))

# load wrapping scripts
source("../../wrapper/R/AmTrack.R")

# necessary library for plotting
library("lattice")

####################### DOSE(DISTANCE) function test #########################################

# radius definitions:
r.m <- 10^seq (-14, -3, length.out=100)
#r.m <- seq( 1e-9,1.5*1e-8, length.out=10)
#r.m <- c(1e-6, 1e-7, 1e-8, 1e-9, 1e-12)

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

df1$D.Gy	<- numeric(nrow(df1))

df1$r.m.check	<- numeric(nrow(df1))


material.no <- c(1)
E.MeV.u = rep( 100 , nrow(df1))
particle.no = rep( 1001, nrow(df1))


for( i in ER.model) {
 ii 			<- df1$ER.model == i
	for( j in RDD.model) {
 	jj 			<- ((df1$RDD.model == j) & (df1$ER.model == i)) 
			df1$D.Gy[jj] <- AT.D.RDD.Gy(r.m = df1$r.m[jj], E.MeV.u, particle.no, material.no, ER.model = i, RDD.model = j, RDD.parameters=RDD.parameters[[j]] )
 		df1$r.m.check[jj] <- AT.r.RDD.m(D.Gy = df1$D.Gy[jj], E.MeV.u, particle.no, material.no, ER.model = i, RDD.model = j, RDD.parameters=RDD.parameters[[j]] ) 	
		}
}

#df1

# plots...

p1 <- xyplot( D.Gy ~ r.m | ER.model.name , groups = RDD.model.name, ref = TRUE, data=df1, pch = ".", lty = 1, type = "l", xlab = "Distance [m]", ylab = "Dose [Gy]", auto.key = list(title = "Protons in liquid water, RDD",points = FALSE, lines = TRUE), scales = list(log = 10))
p2 <- xyplot( D.Gy ~ r.m | RDD.model.name, groups = ER.model.name, ref = TRUE, data=df1, pch = ".", lty = 1, type = "l", xlab = "Distance [m]", ylab = "Dose [Gy]", auto.key = list(title = "Protons in liquid water, RDD",points = FALSE, lines = TRUE), scales = list(log = 10))

p1
p2

p1.check <- xyplot( r.m.check ~ r.m | ER.model.name , groups = RDD.model.name, ref = TRUE, data=df1, pch = ".", lty = 1, type = "l", xlab = "Distance [m]", ylab = "Distance [m]", auto.key = list(title = "Protons in liquid water, test RDD",points = FALSE, lines = TRUE), scales = list(log = 10))
p2.check <- xyplot( r.m.check ~ r.m | RDD.model.name, groups = ER.model.name, ref = TRUE, data=df1, pch = ".", lty = 1, type = "l", xlab = "Distance [m]", ylab = "Distance [m]", auto.key = list(title = "Protons in liquid water, test RDD",points = FALSE, lines = TRUE), scales = list(log = 10))


####################### DISTANCE(DOSE) function test #########################################

# energy range definitions:
#D.Gy <- c( 1e-8,1e-6,1e-4)
#D.Gy <- c( 7.115367e-02, 7.116197e+02)
#D.expn <- seq (-2, 7, by=0.1)
D.Gy <- 10^seq(-10,6,length.out=50)


df2 <- expand.grid( D.Gy = D.Gy, ER.model = ER.model, RDD.model = RDD.model )

df2$ER.model.name	<- as.character(ER.model.names[df2$ER.model])
df2$RDD.model.name	<- as.character(RDD.model.names[df2$RDD.model])
df2$r.m	<- numeric(nrow(df2))

E.MeV.u = rep( 100, nrow(df2))
particle.no = rep( 1001, nrow(df2))

#df2

#RDD.parameters

for( i in ER.model) {
 ii 			<- df2$ER.model == i
 for( j in RDD.model) {
 	jj 			<- ((df2$RDD.model == j) & (df2$ER.model == i)) 
  		df2$r.m[jj] <- AT.r.RDD.m(D.Gy = df2$D.Gy[jj], E.MeV.u, particle.no, material.no, ER.model = i, RDD.model = j, RDD.parameters=RDD.parameters[[j]] ) 	
		}
}

#df2

# plots...

p3 <- xyplot( r.m ~ D.Gy | ER.model.name , groups = RDD.model.name, ref = TRUE, data=df2, pch = ".", lty = 1, type = "l", xlab = "Dose [Gy]", ylab = "Distance [m]", auto.key = list(title = "Protons in liquid water, inverse RDD",points = FALSE, lines = TRUE), scales = list(log = 10))
p4 <- xyplot( r.m ~ D.Gy | RDD.model.name, groups = ER.model.name, ref = TRUE, data=df2, pch = ".", lty = 1, type = "l", xlab = "Dose [Gy]", ylab = "Distance [m]", auto.key = list(title = "Protons in liquid water, inverse RDD",points = FALSE, lines = TRUE), scales = list(log = 10))

# saving plots to the file

pdf("RDD.pdf")

p1
p2
p3
p4
p1.check
p2.check

dev.off()
}