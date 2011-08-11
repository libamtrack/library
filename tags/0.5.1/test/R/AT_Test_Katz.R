################################################################################################
# R test script for implemented Katz model
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
recompile <- TRUE
source("AT_Test_PreRun.R")

require(lattice)

####################### INACTIVATION PROBABILITY function test #########################################


# radius definitions:
r.m <- 10^seq (-10, -3, length.out=3)

# electron range models definition
ER.model.names <- c("simple test",  "Butts & Katz' (linear)",  "Waligorski's (power-law wmax)",  "Geiss' (power-law E)", "Scholz' (power-law E)", "Edmund' (power-law wmax)","Tabata")
ER.model <- c(2,3,7)

# RDD models definition
RDD.model.names <- c("Simple step test function",  "Katz' point target", "Geiss'", "Site", "Cucinotta", "KatzExtTarget", "CucinottaExtTarget")
RDD.model <- c(6,7)

# RDD parameters
RDD.parameters <- list(c(1),c(1e-10,1e-10), c(1e-8),c(1e-8,1e-10),c(5e-11,1e-10),c(1e-10,1e-8,1e-10),c(5e-11,1e-8,1e-10))

# data frame setup
df1 <- expand.grid( r.m = r.m, ER.model = ER.model, RDD.model = RDD.model )

df1$ER.model.name	<- as.character(ER.model.names[df1$ER.model])
df1$RDD.model.name	<- as.character(RDD.model.names[df1$RDD.model])

df1$inact.prob	<- numeric(nrow(df1))

material.no  <-  c(1)
E.MeV.u      <-  c(100.0)
particle.no  <-  c(1001)

# 1, D0, c, m, 0
GR.parameters <- c(1, 3.0, 1, 2, 0)
SP.source <- 0 # (0 -> PSTAR)

for( i in ER.model) {
 ii 			<- df1$ER.model == i
	for( j in RDD.model) {
 	jj 			<- ((df1$RDD.model == j) & (df1$ER.model == i)) 
			df1$inact.prob[jj] <- AT.KatzModel.inactivation.probability(r.m = df1$r.m[jj], E.MeV.u, particle.no, material.no, rdd.model = j, rdd.parameters = RDD.parameters[[j]], er.model = i, GR.parameters, stop.power.source = SP.source)[[1]]
		}
}

df1

# plots...

p1 <- xyplot( inact.prob ~ r.m | ER.model.name , groups = RDD.model.name, ref = TRUE, data=df1, pch = ".", lty = 1, type = "l", xlab = "Distance [m]", ylab = "Inactivation prob.", auto.key = list(title = "Protons in liquid water",points = FALSE, lines = TRUE), scales = list(log = 10))
p2 <- xyplot( inact.prob ~ r.m | RDD.model.name, groups = ER.model.name, ref = TRUE, data=df1, pch = ".", lty = 1, type = "l", xlab = "Distance [m]", ylab = "Inactivation prob.", auto.key = list(title = "Protons in liquid water",points = FALSE, lines = TRUE), scales = list(log = 10))

p1logx <- xyplot( inact.prob ~ r.m | ER.model.name , groups = RDD.model.name, ref = TRUE, data=df1, pch = ".", lty = 1, type = "l", xlab = "Distance [m]", ylab = "Inactivation prob.", auto.key = list(title = "Protons in liquid water",points = FALSE, lines = TRUE), scales = list( x = list(log = 10)))
p2logx <- xyplot( inact.prob ~ r.m | RDD.model.name, groups = ER.model.name, ref = TRUE, data=df1, pch = ".", lty = 1, type = "l", xlab = "Distance [m]", ylab = "Inactivation prob.", auto.key = list(title = "Protons in liquid water",points = FALSE, lines = TRUE), scales = list( x = list(log = 10)))

####################### INACTIVATION CROSS SECTION function test #########################################

E.MeV.u <- 10^seq (0, 3, length.out=10)

# data frame setup
df1 <- expand.grid( E.MeV.u = E.MeV.u, ER.model = ER.model, RDD.model = RDD.model )

df1$ER.model.name	<- as.character(ER.model.names[df1$ER.model])
df1$RDD.model.name	<- as.character(RDD.model.names[df1$RDD.model])

df1$inact.cross.sect.m2	<- numeric(nrow(df1))

material.no  <-  c(1)
particle.no  <-  c(1001)

# 1, D0, c, m, 0
GR.parameters <- c(1, 3.0, 1, 2, 0)

for( i in ER.model) {
 ii 			<- df1$ER.model == i
	for( j in RDD.model) {
 	jj 			<- ((df1$RDD.model == j) & (df1$ER.model == i)) 
			df1$inact.cross.sect.m2[jj] <- AT.KatzModel.inactivation.cross.section.m2(E.MeV.u = df1$E.MeV.u[jj], particle.no, material.no, rdd.model = j, rdd.parameters = RDD.parameters[[j]], er.model = i, GR.parameters, stop.power.source = SP.source)[[1]]
		}
}

#df1

p3 <- xyplot( inact.cross.sect.m2 ~ E.MeV.u | ER.model.name , groups = RDD.model.name, ref = TRUE, data=df1, pch = ".", lty = 1, type = "l", xlab = "Energy [MeV/u]", ylab = "Inactivation cross section [m2]", auto.key = list(title = "Protons in liquid water",points = FALSE, lines = TRUE), scales = list(log = 10))
p4 <- xyplot( inact.cross.sect.m2 ~ E.MeV.u | RDD.model.name, groups = ER.model.name, ref = TRUE, data=df1, pch = ".", lty = 1, type = "l", xlab = "Energy [MeV/u]", ylab = "Inactivation cross section [m2]", auto.key = list(title = "Protons in liquid water",points = FALSE, lines = TRUE), scales = list(log = 10))

p3logx <- xyplot( inact.cross.sect.m2 ~ E.MeV.u | ER.model.name , groups = RDD.model.name, ref = TRUE, data=df1, pch = ".", lty = 1, type = "l", xlab = "Energy [MeV/u]", ylab = "Inactivation cross section [m2]", auto.key = list(title = "Protons in liquid water",points = FALSE, lines = TRUE), scales = list( x = list(log = 10)))
p4logx <- xyplot( inact.cross.sect.m2 ~ E.MeV.u | RDD.model.name, groups = ER.model.name, ref = TRUE, data=df1, pch = ".", lty = 1, type = "l", xlab = "Energy [MeV/u]", ylab = "Inactivation cross section [m2]", auto.key = list(title = "Protons in liquid water",points = FALSE, lines = TRUE), scales = list( x = list(log = 10)))

####################### SURVIVAL function test #########################################

D.Gy <- 10^seq (-2, 1, length.out=10)
E.MeV.u <- 60

# data frame setup
df1 <- expand.grid( D.Gy = D.Gy, ER.model = ER.model, RDD.model = RDD.model )

df1$ER.model.name   <- as.character(ER.model.names[df1$ER.model])
df1$RDD.model.name  <- as.character(RDD.model.names[df1$RDD.model])

df1$survival <- numeric(nrow(df1))

material.no  <-  c(1)
particle.no  <-  c(6012)
sigma0.m2 <- 10^(-10)
E.array <- rep( E.MeV.u, length(D.Gy) )
particle.array <- rep( particle.no, length(D.Gy) )

# 1, D0, c, m, 0
GR.parameters <- c(1, 3.0, 1, 2, 0)

for( i in ER.model) {
 ii             <- df1$ER.model == i
    for( j in RDD.model) {
    jj          <- ((df1$RDD.model == j) & (df1$ER.model == i))
            for( k in D.Gy ) {
                kk   <- ((df1$RDD.model == j) & (df1$ER.model == i) & (df1$D.Gy == k) )
                fluence.cm2 <- AT.fluence.cm2.from.dose.Gy( E.MeV.u, D.Gy = df1$D.Gy[kk], particle.no, material.no, SP.source)[[1]]                 
                df1$survival[kk] <- AT.KatzModel.single.field.survival( fluence.cm2, E.MeV.u, particle.no, material.no,
                rdd.model = j,
                rdd.parameters = RDD.parameters[[j]],
                er.model = i,
                D0.characteristic.dose.Gy = GR.parameters[[2]],
                m.number.of.targets = GR.parameters[[4]],
                sigma0.m2 = sigma0.m2,
                stopping.power.source.no = SP.source)[[1]]
            }
        }
}

df1

p5logy <- xyplot( survival ~ D.Gy | ER.model.name , groups = RDD.model.name, ref = TRUE, data=df1, pch = ".", lty = 1, type = "l", xlab = "Dose [Gy]", ylab = "Survival", auto.key = list(title = "Carbons in liquid water",points = FALSE, lines = TRUE), scales = list( y = list(log = 10)))
p6logy <- xyplot( survival ~ D.Gy | RDD.model.name, groups = ER.model.name, ref = TRUE, data=df1, pch = ".", lty = 1, type = "l", xlab = "Dose [Gy]", ylab = "Survival", auto.key = list(title = "Carbons in liquid water",points = FALSE, lines = TRUE), scales = list( y = list(log = 10)))



pdf("Katz.pdf")

#p1
#p2
#p1logx
p2logx

#p3
p4
#p3logx
p4logx

p5logy
p6logy


dev.off()
