################################################################################################
# R test script for implemented electron range models
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
dyn.load("../../AT_Release/libamtrack.dll")

# load wrapping scripts
source("../../wrapper/R/AmTrack.R")

# necessary library for plotting
library("lattice")

# energy range definitions:
expn <- seq (-3, 4, by=1e-1)
E.MeV.u <- 10^expn

# models definition
er.models.names <- c("simple test ER model",  "Butts & Katz' ER model (linear)",  "Waligorski's ER model (power-law wmax)",  "Geiss' ER model (power-law E)", "Scholz' ER model (power-law E)", "Edmund' ER model (power-law wmax)","Tabata  ER model")
er.models <- c(1,2,3,4,5,6,7)

# other parameters
n <- length(E.MeV.u) 
material.number <- rep(1, 1) # Water

# data frame setup
df <- expand.grid( E.MeV.u = E.MeV.u , er.models = er.models )

df$er.models.name  <-  as.character(er.models.names[df$er.models])
df$range.m         <-  numeric(nrow(df))
df$wmax            <-  0

# calculations
j <- 0
for( i in er.models ){
  j                <-  j+1
  ii			   <-  df$er.models == i
  df$er.models.name[ii] <- er.models.names[j]
  wmax_MeV         <-  AT.max.E.transfer.MeV( E.MeV.u = df$E.MeV.u[ii] )
  df$wmax[ii]      <-  wmax_MeV
  df$range.m[ii]   <-  AT.max.electron.range(	E.MeV.u = df$E.MeV.u[ii], material.number, i)
}

# plots...

logplot <- xyplot( 1e2*range.m ~ wmax, groups = er.models.name, ref = TRUE, data=df, pch = ".", lty = 1, type = "l", xlab = "wmax [MeV]", ylab = "Range [cm]", auto.key = list(title = "Range of delta electrons in liquid water",points = FALSE, lines = TRUE), scales = list(log = 10))
linplot <- xyplot( 1e2*range.m ~ wmax, groups = er.models.name, ref = TRUE, data=df, pch = ".", lty = 1, type = "l", xlab = "wmax [MeV]", ylab = "Range [cm]", auto.key = list(title = "Range of delta electrons in liquid water",points = FALSE, lines = TRUE))


pdf("ER.pdf")

logplot
linplot

dev.off()