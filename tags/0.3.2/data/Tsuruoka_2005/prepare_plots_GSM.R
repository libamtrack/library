################################################################################################
# <file description>
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

# clear workspace
rm( list = ls() )

# load libAmTrack library
dyn.load("../../Release/libAmTrack.dll")

# load wrapping scripts
source("../../WrappingScripts/AmTrack_S.ssc")

# load data access script
source("data_operations.R")

# necessary libraries for plotting
library("lattice")
library("gplots")

############ useful functions ####################

AT.particle.no.from.Z <- function( Z ){
   if( Z == 6 )
   		return (18)
   if( Z == 10 )
   		return (33)
   if( Z == 14 )
   		return (34)
   if( Z == 26 )
   		return (35)
   return (-1)
}

AT.particle.name.from.Z <- function( Z ){
   if( Z == 6 )
   		return ("Carbon")
   if( Z == 10 )
   		return ("Neon")
   if( Z == 14 )
   		return ("Silicon")
   if( Z == 26 )
   		return ("Iron")
   return ("")
}


############ main script #########################

# read data into data frames
experimental.data = AT.survival.data.read()

cond  <- (experimental.data$radiation == "gamma")
gamma.data <- data.frame( D.Gy = experimental.data$D.Gy[ cond ], Survival = experimental.data$Survival[ cond ], Sur.min = experimental.data$Sur.min[ cond ] , Sur.max = experimental.data$Sur.max[ cond ] , type = "gamma")   
gamma.LQ.params <- AT.fit.linquad( gamma.data$D.Gy, gamma.data$Survival, gamma.data$Sur.min, gamma.data$Sur.max, 0., 0. )

# model setup - physics
material.no      <- 1                       # liquid water
RDD.model        <- 3                       # [ 3 - Geiss RDD ]
RDD.parameters   <- 1e-8                    # for Geiss RDD this means a0
ER.model         <- 2                       # [2 - Waligorski ER]
ER.parameters    <- 1
gamma.model      <- 5                       # LQ model
gamma.parameters <- c(gamma.LQ.params,20.0) # Gamma parameters + Dt = 20 Gy


# plotting
pdf("GSM_SurvivalCurves.pdf", width=15, height=15, pointsize=24)

simulation.results <- read.table("GSM_results.dat")

LQ.D.Gy = seq( 0, 6, by = 0.01)

Z_init.energy <- unique(paste( simulation.results$Z , simulation.results$Init.Energy, sep = "-"))
# loop over all datasets
for( z_en in Z_init.energy){
  cond.z_en <- (paste( simulation.results$Z , simulation.results$Init.Energy, sep = "-") == z_en)
  Z <- unique(simulation.results$Z[cond.z_en])
  Init.Energy <- unique(simulation.results$Init.Energy[cond.z_en])
  LET.vector <- unique(simulation.results$LET[cond.z_en])
  particle.no <- AT.particle.no.from.Z(Z) # ion number
#  if ((particle.no > 0) & (Z == 6) & (Init.Energy == 290)) {
  if (particle.no > 0) {
   cat("**************************************************\n")
   cat("Z = ",Z,"\n")
   cat("Init.Energy = ",Init.Energy,"\n")
   # loop over all possible LET values
   for( LET in LET.vector ){
   
      ###################### data preparation ###############################
     cat("******* LET = ",LET,"**********\n")
     E.MeV.u <- AT.E.MeV.u( LET , particle.no, material.no )[1] # energy from LET
     cond.LET.GSM  <- ((simulation.results$Z == Z) & (simulation.results$Init.Energy == Init.Energy) & (simulation.results$LET == LET))
     cond.LET.exp  <- ((experimental.data$Z == Z) & (experimental.data$Init.Energy == Init.Energy) & (experimental.data$LET == LET))     

     # read data from experimental dataset
     D.Gy.exp <- experimental.data$D.Gy[cond.LET.exp]
     Survival.exp.ion <- experimental.data$Survival[cond.LET.exp]
     Survival.min.exp.ion <- experimental.data$Sur.min[cond.LET.exp]
     Survival.max.exp.ion <- experimental.data$Sur.max[cond.LET.exp]

     # read dose values from model simulation
     D.Gy <- simulation.results$D.Gy[cond.LET.GSM]
     
     # read dose values from model simulation (gamma radiation)
     Survival.GSM.gamma <- simulation.results$Survival.gamma[cond.LET.GSM]
					
     # fit LQ model to simulated gamma response
     GSM.gamma.params <- AT.fit.linquad( D.Gy, Survival.GSM.gamma, rep(1,length(D.Gy)), 0., 0. )
     cat( "GSM gamma: alpha = ", GSM.gamma.params[1], " beta = ",  GSM.gamma.params[2],"\n")
     
     # prepare survival curve with LQ parameters for simulated gamma response
     Survival.GSM.gamma.LQ <- 100.*exp( - GSM.gamma.params[1] * LQ.D.Gy - GSM.gamma.params[2] * LQ.D.Gy * LQ.D.Gy)
     
     # read dose values from model simulation (ion radiation)
     Survival.GSM.ion <- simulation.results$Survival.ion[cond.LET.GSM]
     
     # fit LQ model to simulated ion response
     GSM.ion.params <- AT.fit.linquad( D.Gy, Survival.GSM.ion, rep(1,length(D.Gy)), 0., 0. )
     cat( "GSM ion: alpha = ", GSM.ion.params[1], " beta = ",  GSM.ion.params[2],"\n")

     # prepare survival curve with LQ parameters for simulated ion response
     Survival.GSM.ion.LQ <- 100.*exp( - GSM.ion.params[1] * LQ.D.Gy - GSM.ion.params[2] * LQ.D.Gy * LQ.D.Gy)
     
     # dose reproducibility factor
     Dose.Factor.GSM <- simulation.results$Dose.Factor[cond.LET.GSM]
     
     ###################### plotting ###############################
     
     # plot gamma experimental data
     plotCI( x = gamma.data$D.Gy , y = gamma.data$Survival/100., ui = gamma.data$Sur.max/100., li = gamma.data$Sur.min/100., gap = 0 , xlab = "Dose (Gy)" , ylab = "Surviving fraction" , log = "y" , axis = FALSE, pch = 1, xlim = c(0,6), ylim = c(0.001, 1.), col = "green")

     # plot ion experimental data     
     par(new=TRUE)
     plotCI( x = D.Gy.exp , y = Survival.exp.ion/100., ui = Survival.max.exp.ion/100., li = Survival.min.exp.ion/100., gap = 0 , xlab = "" , ylab = "" , log = "y" , axis = FALSE, pch = 2, xlim = c(0,6), ylim = c(0.001, 1.), col = "red")
     
     # plot gamma GSM simulated data and LQ fit
     points( x = D.Gy, y = Survival.GSM.gamma/100., pch = 3, col = "green")
     lines( x = LQ.D.Gy , y = Survival.GSM.gamma.LQ/100. , col = "green")
     
     # plot ion GSM simulated data and LQ fit
     points( x = D.Gy, y = Survival.GSM.ion/100., pch = 4, col = "red")     
     lines( x = LQ.D.Gy , y = Survival.GSM.ion.LQ/100. , col = "red")
     
     # plot of dose reproducibility factor
     points( x = D.Gy , y = Dose.Factor.GSM )
     abline()
     
     # refinement
     abline( h = c(1,0.1) ,lw = 2)
     axis( 1, at = c(0,1,2,3,4,5,6), lw = 4)
     axis( 2, at = c(0.001,0.01,0.1,1), lw = 4)
     legend( "topright", legend = c("exp. data gamma","model gamma","exp. data ions","model ions"), pch=c(1,3,2,4), col=c("green","green","red","red"), lw = 4)
     title(paste("Z = ",Z," ,",AT.particle.name.from.Z(Z), " ions , ",Init.Energy, " MeV/u, LET = ", LET/10, "keV/um" ))
     
   } # end LET loop
  } # end if
} # end dataset loop

dev.off()
