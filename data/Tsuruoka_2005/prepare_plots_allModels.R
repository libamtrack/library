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

####################################### PLOT 1 ######################################################

# plotting
pdf("models_RBE.pdf", width=15, height=15, pointsize=24)

SPIFF.results <- read.table("SPIFF_results.dat")
Katz.results <- read.table("Katz_results.dat")
GSM.results <- read.table("GSM_results.dat")

RBE.article <- read.table("RBE_article.dat")

plotCI( x = 1, y = 1, xlim = c(0,500), ylim = c(0, 5), axis = FALSE, xlab = "", ylab = "")

Z_init.energy <- unique(paste( SPIFF.results$Z , SPIFF.results$Init.Energy, sep = "-"))
# loop over all datasets
for( z_en in Z_init.energy){
  cond.z_en.SPIFF <- (paste( SPIFF.results$Z , SPIFF.results$Init.Energy, sep = "-") == z_en)
  Z <- unique(SPIFF.results$Z[cond.z_en.SPIFF])
  Init.Energy <- unique(SPIFF.results$Init.Energy[cond.z_en.SPIFF])
  LET.vector.SPIFF <- unique(SPIFF.results$LET[cond.z_en.SPIFF])
  particle.no <- AT.particle.no.from.Z(Z) # ion number
  if (particle.no > 0) {
   cat("**************************************************\n")
   cat("Z = ",Z,"\n")
   cat("Init.Energy = ",Init.Energy,"\n")
   # loop over all possible LET values
   
   cond.RBE.article <- ((RBE.article$Z == Z) & (RBE.article$Init.Energy == Init.Energy))
   
   LET <- RBE.article$LET[cond.RBE.article]
   LET.min <- RBE.article$LET.min[cond.RBE.article]
   LET.max <- RBE.article$LET.max[cond.RBE.article]
   RBE <- RBE.article$RBE[cond.RBE.article]
   RBE.err  <- RBE.article$RBE.err[cond.RBE.article]
   RBE.min  <- RBE - RBE.err
   RBE.max  <- RBE + RBE.err
   
   par(new=TRUE)
   plotCI( x = LET, y = RBE, ui = RBE.max, li = RBE.min, err = "y", gap = 0., xlim = c(0,500), ylim = c(0, 5), axis = FALSE, xlab = "LET [keV/um]", ylab = "RBE")
   par(new=TRUE)
   plotCI( x = LET, y = RBE, ui = LET.max, li = LET.min, err = "x", gap = 0., xlim = c(0,500), ylim = c(0, 5), axis = FALSE, xlab = NA, ylab = NA, add = TRUE)
   #title(paste("Z = ",Z," ,",AT.particle.name.from.Z(Z), " ions , ",Init.Energy, " MeV/u" ))
   legend( "topright", legend = c("exp. data - article","exp. data - LQ","model - SPIFF","model - Katz"), pch=c(1,1,2,3), col=c("black","green","red","blue"), lw = 4)

   RBE.vector.SPIFF <- numeric()
   LET.vector.SPIFF <- sort( LET.vector.SPIFF ) 

   for( LET in LET.vector.SPIFF ){
      ###################### data preparation ###############################
     cat("******* LET = ",LET,"**********\n")
     E.MeV.u <- AT.E.MeV.u( LET , particle.no, material.no )[1] # energy from LET
     cond.LET.SPIFF  <- ((SPIFF.results$Z == Z) & (SPIFF.results$Init.Energy == Init.Energy) & (SPIFF.results$LET == LET))

     # read dose values from model simulation
     D.Gy <- SPIFF.results$D.Gy[cond.LET.SPIFF]
     
     # read dose values from model simulation (gamma radiation)
     Survival.SPIFF.gamma <- SPIFF.results$Survival.gamma[cond.LET.SPIFF]
					
     # fit LQ model to simulated gamma response
     SPIFF.gamma.params <- AT.fit.linquad( D.Gy, Survival.SPIFF.gamma, rep(1,length(D.Gy)), 0., 0. )
     cat( "SPIFF gamma: alpha = ", SPIFF.gamma.params[1], " beta = ",  SPIFF.gamma.params[2],"\n")
     
     # read dose values from model simulation (ion radiation)
     Survival.SPIFF.ion <- SPIFF.results$Survival.ion[cond.LET.SPIFF]
     
     # fit LQ model to simulated ion response
     SPIFF.ion.params <- AT.fit.linquad( D.Gy, Survival.SPIFF.ion, rep(1,length(D.Gy)), 0., 0. )
     cat( "SPIFF ion: alpha = ", SPIFF.ion.params[1], " beta = ",  SPIFF.ion.params[2],"\n")
     
     # calculate simulated (SPIFF) RBE value 
     RBE.SPIFF <- AT.RBE( SPIFF.ion.params[1], SPIFF.ion.params[2], SPIFF.gamma.params[1], SPIFF.gamma.params[2], survival = 10)
     
     RBE.vector.SPIFF <- c( RBE.vector.SPIFF , RBE.SPIFF)
          
   } # end LET loop

   ###################### plotting ###############################
   lines( x = LET.vector.SPIFF/10., y = RBE.vector.SPIFF, col = "red")


   cond.Z.en.Katz  <- ((experimental.data$Z == Z) & (experimental.data$Init.Energy == Init.Energy))
   LET.vector.Katz <- unique(experimental.data$LET[cond.Z.en.Katz])

   for( LET in LET.vector.Katz ){
      ###################### data preparation ###############################
     cat("******* LET = ",LET,"**********\n")
     E.MeV.u <- AT.E.MeV.u( LET , particle.no, material.no )[1] # energy from LET
     cond.LET.Katz  <- ((Katz.results$Z == Z) & (Katz.results$Init.Energy == Init.Energy) & (Katz.results$LET == LET))

     # read dose values from model simulation
     D.Gy <- Katz.results$D.Gy[cond.LET.Katz]
     
     # read dose values from model simulation (gamma radiation)
     Survival.Katz.gamma <- Katz.results$Survival.gamma[cond.LET.Katz]
					
     # fit LQ model to simulated gamma response
     Katz.gamma.params <- AT.fit.linquad( D.Gy, Survival.Katz.gamma, rep(1,length(D.Gy)), 0., 0. )
     cat( "Katz gamma: alpha = ", Katz.gamma.params[1], " beta = ",  Katz.gamma.params[2],"\n")
     
     # read dose values from model simulation (ion radiation)
     Survival.Katz.ion <- Katz.results$Survival.ion[cond.LET.Katz]
     
     # fit LQ model to simulated ion response
     Katz.ion.params <- AT.fit.linquad( D.Gy, Survival.Katz.ion, rep(1,length(D.Gy)), 0., 0. )
     cat( "Katz ion: alpha = ", Katz.ion.params[1], " beta = ",  Katz.ion.params[2],"\n")
     
     # calculate simulated (Katz) RBE value 
     RBE.Katz <- AT.RBE( Katz.ion.params[1], Katz.ion.params[2], Katz.gamma.params[1], Katz.gamma.params[2], survival = 10)

     ###################### plotting ###############################     
     # plot RBE
     if( RBE.Katz > 1 ){
		     points( x = LET/10., y = RBE.Katz, pch = 3, col = "blue")
		   }     
   } # end LET loop

   cond.Z.en.GSM  <- ((GSM.results$Z == Z) & (GSM.results$Init.Energy == Init.Energy))
   LET.vector.GSM <- unique(GSM.results$LET[cond.Z.en.GSM])

   for( LET in LET.vector.GSM ){
      ###################### data preparation ###############################
     cat("******* LET = ",LET,"**********\n")
     E.MeV.u <- AT.E.MeV.u( LET , particle.no, material.no )[1] # energy from LET
     cond.LET.GSM  <- ((GSM.results$Z == Z) & (GSM.results$Init.Energy == Init.Energy) & (GSM.results$LET == LET))

     # read dose values from model simulation
     D.Gy <- GSM.results$D.Gy[cond.LET.GSM]
        
     # read dose values from model simulation (ion radiation)
     Survival.GSM.ion <- GSM.results$Survival.ion[cond.LET.GSM]
     
     # fit LQ model to simulated ion response
     GSM.ion.params <- AT.fit.linquad( D.Gy, Survival.GSM.ion, rep(1,length(D.Gy)), 0., 0. )
     cat( "GSM ion: alpha = ", GSM.ion.params[1], " beta = ",  GSM.ion.params[2],"\n")
     
     # calculate simulated (Katz) RBE value 
     RBE.GSM <- AT.RBE( GSM.ion.params[1], GSM.ion.params[2], gamma.LQ.params[1], gamma.LQ.params[2], survival = 10)

     ###################### plotting ###############################     
     # plot RBE
     points( x = LET/10., y = RBE.GSM, pch = 4, col = "pink")     
   } # end LET loop

   cond.Z.en.exp  <- ((experimental.data$Z == Z) & (experimental.data$Init.Energy == Init.Energy))
   LET.vector.exp <- unique(experimental.data$LET[cond.Z.en.exp])

   for( LET in LET.vector.exp ){
      ###################### data preparation ###############################
     cat("******* LET = ",LET,"**********\n")
     E.MeV.u <- AT.E.MeV.u( LET , particle.no, material.no )[1] # energy from LET
     cond.LET.exp  <- ((experimental.data$Z == Z) & (experimental.data$Init.Energy == Init.Energy) & (experimental.data$LET == LET))     

     # read data from experimental dataset
     D.Gy.exp <- experimental.data$D.Gy[cond.LET.exp]
     Survival.exp.ion <- experimental.data$Survival[cond.LET.exp]
     Survival.min.exp.ion <- experimental.data$Sur.min[cond.LET.exp]
     Survival.max.exp.ion <- experimental.data$Sur.max[cond.LET.exp]
     Survival.err.exp.ion <- experimental.data$Sur.err[cond.LET.exp]

     # fit LQ model to experimental ion response
     exp.ion.params <- AT.fit.linquad( D.Gy.exp, Survival.exp.ion, Survival.err.exp.ion , 0., 0. )
     cat( "Experiment ion: alpha = ", exp.ion.params[1], " beta = ",  exp.ion.params[2],"\n")
     
     # calculate experimental RBE value 
     RBE.exp <- AT.RBE( exp.ion.params[1], exp.ion.params[2], gamma.LQ.params[1], gamma.LQ.params[2], survival = 10)

     ###################### plotting ###############################     
     # plot RBE
     points( x = LET/10., y = RBE.exp, pch = 1, col = "green")     
   } # end LET loop
  } # end if
} # end dataset loop

dev.off()


####################################### PLOT 2 ######################################################

# plotting
pdf("models_RBE_details.pdf", width=15, height=15, pointsize=24)

SPIFF.results <- read.table("SPIFF_results.dat")
Katz.results <- read.table("Katz_results.dat")

RBE.article <- read.table("RBE_article.dat")

Z_init.energy <- unique(paste( SPIFF.results$Z , SPIFF.results$Init.Energy, sep = "-"))
# loop over all datasets
for( z_en in Z_init.energy){
  cond.z_en.SPIFF <- (paste( SPIFF.results$Z , SPIFF.results$Init.Energy, sep = "-") == z_en)
  Z <- unique(SPIFF.results$Z[cond.z_en.SPIFF])
  Init.Energy <- unique(SPIFF.results$Init.Energy[cond.z_en.SPIFF])
  LET.vector.SPIFF <- unique(SPIFF.results$LET[cond.z_en.SPIFF])
  particle.no <- AT.particle.no.from.Z(Z) # ion number
  if (particle.no > 0) {
   cat("**************************************************\n")
   cat("Z = ",Z,"\n")
   cat("Init.Energy = ",Init.Energy,"\n")
   # loop over all possible LET values
   
   cond.RBE.article <- ((RBE.article$Z == Z) & (RBE.article$Init.Energy == Init.Energy))
   
   LET <- RBE.article$LET[cond.RBE.article]
   LET.min <- RBE.article$LET.min[cond.RBE.article]
   LET.max <- RBE.article$LET.max[cond.RBE.article]
   RBE <- RBE.article$RBE[cond.RBE.article]
   RBE.err  <- RBE.article$RBE.err[cond.RBE.article]
   RBE.min  <- RBE - RBE.err
   RBE.max  <- RBE + RBE.err
   
   plotCI( x = LET, y = RBE, ui = RBE.max, li = RBE.min, err = "y", gap = 0., xlim = c(0,max(LET.max)), ylim = c(0, max(RBE.max)+1), axis = FALSE, xlab = "LET [keV/um]", ylab = "RBE")
   plotCI( x = LET, y = RBE, ui = LET.max, li = LET.min, err = "x", gap = 0., xlim = c(0,max(LET.max)), ylim = c(0, max(RBE.max)+1), axis = FALSE, xlab = NA, ylab = NA, add = TRUE)
   title(paste("Z = ",Z," ,",AT.particle.name.from.Z(Z), " ions , ",Init.Energy, " MeV/u" ))
   legend( "bottomright", legend = c("exp. data - article","exp. data - LQ","model - SPIFF","model - Katz","model - GSM"), pch=c(1,1,2,3,4), col=c("black","green","red","blue","pink"), lw = 4)

   RBE.vector.SPIFF <- numeric()
   LET.vector.SPIFF <- sort( LET.vector.SPIFF ) 

   for( LET in LET.vector.SPIFF ){
      ###################### data preparation ###############################
     cat("******* LET = ",LET,"**********\n")
     E.MeV.u <- AT.E.MeV.u( LET , particle.no, material.no )[1] # energy from LET
     cond.LET.SPIFF  <- ((SPIFF.results$Z == Z) & (SPIFF.results$Init.Energy == Init.Energy) & (SPIFF.results$LET == LET))

     # read dose values from model simulation
     D.Gy <- SPIFF.results$D.Gy[cond.LET.SPIFF]
     
     # read dose values from model simulation (gamma radiation)
     Survival.SPIFF.gamma <- SPIFF.results$Survival.gamma[cond.LET.SPIFF]
					
     # fit LQ model to simulated gamma response
     SPIFF.gamma.params <- AT.fit.linquad( D.Gy, Survival.SPIFF.gamma, rep(1,length(D.Gy)), 0., 0. )
     cat( "SPIFF gamma: alpha = ", SPIFF.gamma.params[1], " beta = ",  SPIFF.gamma.params[2],"\n")
     
     # read dose values from model simulation (ion radiation)
     Survival.SPIFF.ion <- SPIFF.results$Survival.ion[cond.LET.SPIFF]
     
     # fit LQ model to simulated ion response
     SPIFF.ion.params <- AT.fit.linquad( D.Gy, Survival.SPIFF.ion, rep(1,length(D.Gy)), 0., 0. )
     cat( "SPIFF ion: alpha = ", SPIFF.ion.params[1], " beta = ",  SPIFF.ion.params[2],"\n")
     
     # calculate simulated (SPIFF) RBE value 
     RBE.SPIFF <- AT.RBE( SPIFF.ion.params[1], SPIFF.ion.params[2], SPIFF.gamma.params[1], SPIFF.gamma.params[2], survival = 10)
     
     RBE.vector.SPIFF <- c( RBE.vector.SPIFF , RBE.SPIFF)
          
   } # end LET loop

   ###################### plotting ###############################
   lines( x = LET.vector.SPIFF/10., y = RBE.vector.SPIFF, col = "red", lw = 2)


   cond.Z.en.Katz  <- ((experimental.data$Z == Z) & (experimental.data$Init.Energy == Init.Energy))
   LET.vector.Katz <- unique(experimental.data$LET[cond.Z.en.Katz])
   RBE.vector.Katz <- numeric()

   for( LET in LET.vector.Katz ){
      ###################### data preparation ###############################
     cat("******* LET = ",LET,"**********\n")
     E.MeV.u <- AT.E.MeV.u( LET , particle.no, material.no )[1] # energy from LET
     cond.LET.Katz  <- ((Katz.results$Z == Z) & (Katz.results$Init.Energy == Init.Energy) & (Katz.results$LET == LET))

     # read dose values from model simulation
     D.Gy <- Katz.results$D.Gy[cond.LET.Katz]
     
     # read dose values from model simulation (gamma radiation)
     Survival.Katz.gamma <- Katz.results$Survival.gamma[cond.LET.Katz]
					
     # fit LQ model to simulated gamma response
     Katz.gamma.params <- AT.fit.linquad( D.Gy, Survival.Katz.gamma, rep(1,length(D.Gy)), 0., 0. )
     cat( "Katz gamma: alpha = ", Katz.gamma.params[1], " beta = ",  Katz.gamma.params[2],"\n")
     
     # read dose values from model simulation (ion radiation)
     Survival.Katz.ion <- Katz.results$Survival.ion[cond.LET.Katz]
     
     # fit LQ model to simulated ion response
     Katz.ion.params <- AT.fit.linquad( D.Gy, Survival.Katz.ion, rep(1,length(D.Gy)), 0., 0. )
     cat( "Katz ion: alpha = ", Katz.ion.params[1], " beta = ",  Katz.ion.params[2],"\n")
     
     # calculate simulated (Katz) RBE value 
     RBE.Katz <- AT.RBE( Katz.ion.params[1], Katz.ion.params[2], Katz.gamma.params[1], Katz.gamma.params[2], survival = 10)
     RBE.vector.Katz <- c( RBE.vector.Katz , RBE.Katz )

   } # end LET loop

   ###################### plotting ###############################     
   # plot RBE
   lines( x = LET.vector.Katz[RBE.vector.Katz > 1]/10., y = RBE.vector.Katz[RBE.vector.Katz > 1], pch = 3, col = "blue", lw = 2)     


   cond.Z.en.GSM  <- ((GSM.results$Z == Z) & (GSM.results$Init.Energy == Init.Energy))
   LET.vector.GSM <- unique(GSM.results$LET[cond.Z.en.GSM])
   RBE.vector.GSM <- numeric()
   
   for( LET in LET.vector.GSM ){
      ###################### data preparation ###############################
     cat("******* LET = ",LET,"**********\n")
     E.MeV.u <- AT.E.MeV.u( LET , particle.no, material.no )[1] # energy from LET
     cond.LET.GSM  <- ((GSM.results$Z == Z) & (GSM.results$Init.Energy == Init.Energy) & (GSM.results$LET == LET))

     # read dose values from model simulation
     D.Gy <- GSM.results$D.Gy[cond.LET.GSM]
        
     # read dose values from model simulation (ion radiation)
     Survival.GSM.ion <- GSM.results$Survival.ion[cond.LET.GSM]
     
     # fit LQ model to simulated ion response
     GSM.ion.params <- AT.fit.linquad( D.Gy, Survival.GSM.ion, rep(1,length(D.Gy)), 0., 0. )
     cat( "GSM ion: alpha = ", GSM.ion.params[1], " beta = ",  GSM.ion.params[2],"\n")
     
     # calculate simulated (Katz) RBE value 
     RBE.GSM <- AT.RBE( GSM.ion.params[1], GSM.ion.params[2], gamma.LQ.params[1], gamma.LQ.params[2], survival = 10)
     RBE.vector.GSM <- c( RBE.vector.GSM , RBE.GSM )

   } # end LET loop

  ###################### plotting ###############################     
   # plot RBE
   lines( x = LET.vector.GSM/10., y = RBE.vector.GSM, pch = 4, col = "pink", lw = 2)     


   cond.Z.en.exp  <- ((experimental.data$Z == Z) & (experimental.data$Init.Energy == Init.Energy))
   LET.vector.exp <- unique(experimental.data$LET[cond.Z.en.exp])

   for( LET in LET.vector.exp ){
      ###################### data preparation ###############################
     cat("******* LET = ",LET,"**********\n")
     E.MeV.u <- AT.E.MeV.u( LET , particle.no, material.no )[1] # energy from LET
     cond.LET.exp  <- ((experimental.data$Z == Z) & (experimental.data$Init.Energy == Init.Energy) & (experimental.data$LET == LET))     

     # read data from experimental dataset
     D.Gy.exp <- experimental.data$D.Gy[cond.LET.exp]
     Survival.exp.ion <- experimental.data$Survival[cond.LET.exp]
     Survival.min.exp.ion <- experimental.data$Sur.min[cond.LET.exp]
     Survival.max.exp.ion <- experimental.data$Sur.max[cond.LET.exp]
     Survival.err.exp.ion <- experimental.data$Sur.err[cond.LET.exp]

     # fit LQ model to experimental ion response
     exp.ion.params <- AT.fit.linquad( D.Gy.exp, Survival.exp.ion, Survival.err.exp.ion , 0., 0. )
     cat( "Experiment ion: alpha = ", exp.ion.params[1], " beta = ",  exp.ion.params[2],"\n")
     
     # calculate experimental RBE value 
     RBE.exp <- AT.RBE( exp.ion.params[1], exp.ion.params[2], gamma.LQ.params[1], gamma.LQ.params[2], survival = 10)

     ###################### plotting ###############################     
     # plot RBE
     points( x = LET/10., y = RBE.exp, pch = 1, col = "green")     
   } # end LET loop
  } # end if
} # end dataset loop

dev.off()
