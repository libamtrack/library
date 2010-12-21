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

# clear workspace
rm( list = ls() )

# Build latest version of libamtrack and load for direct access
recompile <- TRUE
#recompile <- FALSE
source("AT_Test_PreRun.R")

# necessary library for plotting
require("lattice")

# energy range definitions:
E.MeV.u           <- 10^seq(-1, 3, length.out = 50)

# models definition
er.models.names   <- c( "Simple test ER model",  
                        "Butts & Katz' ER model (linear)",  
                        "Waligorski's ER model (power-law)",  
                        "Geiss' ER model", 
                        "Scholz' ER model", 
                        "Edmund' ER model",
                        "Tabata  ER model")
er.models         <- 1:7

material.no       <- 1            # Water, Liquid

# data frame setup
df.range                <- expand.grid( E.MeV.u                 = E.MeV.u , 
                                        er.models               = er.models, 
                                        range.m                 = 0, 
                                        stringsAsFactors        = FALSE)
df.range$er.models.name <- er.models.names[df.range$er.models]

df.energy               <- expand.grid( E.MeV.u                 = E.MeV.u,
                                        which                   = c("non-relativistic", "relativistic"),
                                        max.electron.energy.keV = 0)


# calculations
ii                                     <-  df.energy$which == "relativistic"
df.energy$max.electron.energy.keV[ii]  <-  AT.max.E.transfer.MeV( E.MeV.u = df.energy$E.MeV.u[ii]  )$max.E.transfer.MeV * 1000
df.energy$max.electron.energy.keV[!ii] <-  AT.max.E.transfer.MeV( E.MeV.u = -df.energy$E.MeV.u[!ii])$max.E.transfer.MeV * 1000

for( i in er.models ){
  ii			               <-  df.range$er.models == i
  df.range$range.m[ii]           <-  AT.max.electron.ranges.m( E.MeV.u      = df.range$E.MeV.u[ii], 
                                                               material.no  = material.no, 
                                                               er.model     = i)$max.electron.range.m
}

# Check sanity of computation output
failures <- which(!(df.range$range.m >= 1e-15 & df.range$range.m < 10))
if( length(failures > 0)){
   stop("At least one range value outside range (1e-15 m to 10 m)!")
}

# Plots

plot1 <- xyplot( log10(max.electron.energy.keV) ~ log10(E.MeV.u), 
                 df.energy,
                 type         = "l", 
                 groups       = which,
                 xlab         = "particle energy [MeV/u]", 
                 ylab         = "log max. electron energy [keV]",
                 panel        = function(...){
                                    panel.grid( v = -1, h = -1)
                                    panel.xyplot(...)},
                 scales       = list( x = list( at     = -1:3,
                                                labels = c("0.1", "1", "10", "100", "1000")),
                                      y = list( at     = 0:3,
                                                labels = c("1", "10", "100", "1000"))),
                 aspect       = 1,
                 auto.key     = list(columns = 2))

plot2 <- xyplot( log10(range.m) ~ log10(E.MeV.u), 
                 df.range,
                 type         = "l", 
                 groups       = er.models.name,
                 xlab         = "particle energy [MeV/u]", 
                 ylab         = "log range [m]",
                 aspect       = 1,
                 panel        = function(...){
                                    panel.grid( v = -1, h = -1)
                                    panel.xyplot(...)},
                 scales       = list( x = list( at     = -1:3,
                                                labels = c("0.1", "1", "10", "100", "1000"))),
                 auto.key     = list(space  = 'right',
                                     cex    = 0.5,
                                     points = FALSE,
                                     lines  = TRUE))



pdf("ER.pdf")

plot1
plot2

dev.off()
