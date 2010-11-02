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

############ main script #########################

# read data into data frames
df.Survival = AT.survival.data.read()

# find alpha and beta parameters from gamma survival curve
cond             <- (df.Survival$radiation == "gamma")
gamma.data       <- data.frame( D.Gy = df.Survival$D.Gy[ cond ], Survival = df.Survival$Survival[ cond ], Sur.min = df.Survival$Sur.min[ cond ] , Sur.max = df.Survival$Sur.max[ cond ] , type = "gamma")   
gamma.LQ.params  <- AT.fit.linquad( gamma.data$D.Gy, gamma.data$Survival, gamma.data$Sur.min, gamma.data$Sur.max, 0., 0. )

# model setup - physics
material.no      <- 1                       # liquid water
RDD.model        <- 3                       # [ 3 - Geiss RDD ]
RDD.parameters   <- 1e-8                    # for Geiss RDD this means a0
ER.model         <- 2                       # [2 - Waligorski ER]
ER.parameters    <- 1
gamma.model      <- 5                       # LQ model
gamma.parameters <- c(gamma.LQ.params,20.0) # Gamma parameters + Dt = 20 Gy

# model setup - algorithms
N.runs           <-  10                     # number of runs
n.X              <- 1e2                     # number of pixels on the grid side (grid size = n.X * n.X)
grid.size.m      <- 1e-6                    # linear pixel size (length of side)


dose <- seq(1,6, by=0.1)
n <- length(dose)

simulation.results <- data.frame( D.Gy = 0, Survival.ion = 0, Survival.gamma = 0, Dose.Factor = 0, Z = 0, Init.Energy = 0, LET = 0, model = "GSM")

Z_init.energy <- unique(paste( df.Survival$Z , df.Survival$Init.Energy, sep = "-"))
# loop over all datasets
for( z_en in Z_init.energy){
  cond.z_en <- (paste( df.Survival$Z , df.Survival$Init.Energy, sep = "-") == z_en)
  Z <- unique(df.Survival$Z[cond.z_en])
  Init.Energy <- unique(df.Survival$Init.Energy[cond.z_en])
  LET.vector <- unique(df.Survival$LET[cond.z_en])
  cat("**************************************************\n")
  cat("Z = ",Z,"\n")
  cat("Init.Energy = ",Init.Energy,"\n")
  particle.no <- AT.particle.no.from.Z(Z) # ion number
#  if ((particle.no > 0) && (Z == 6) && (Init.Energy == 290)) {
  if ((particle.no > 0) & (Z != 10)) {
   # loop over all possible LET values
   for( LET in LET.vector ){
     cat("******* LET = ",LET,"**********\n")
  			#cond.LET  <- (df.Survival$Z == Z & df.Survival$Init.Energy == Init.Energy & df.Survival$LET == LET)
  			#ion.data <- data.frame( D.Gy = df.Survival$D.Gy[ cond.LET ], Survival = df.Survival$Survival[ cond.LET ], Sur.min = df.Survival$Sur.min[ cond.LET ] , Sur.max = df.Survival$Sur.max[ cond.LET ] , type = "ion")
  			E.MeV.u <- AT.E.MeV.u( LET , particle.no, material.no )[1] # energy from LET
  			
     # loop over all possible doses
     for( d in dose ){
       fluence.cm2 = -d
       res <- AT.GSM( E.MeV.u,
         particle.no,
         fluence.cm2,
         material.no,
         RDD.model,
         RDD.parameters,
         ER.model,
         ER.parameters,
         gamma.model,
         gamma.parameters,
         method = "grid",
         N.runs,
         N2 = 10,          
         fluence.factor = 1.0,
         write.output = F,
         n.X,
         lethal.events.mode = T,
         grid.size.m)

       current.results <- data.frame( D.Gy = d, Survival.ion = 100*res[3], Survival.gamma = 100*res[4], Dose.Factor = res[2]/d, Z = Z, Init.Energy = Init.Energy, LET = LET, model = "GSM")
       simulation.results <- rbind( simulation.results, current.results)
         
  			} # end dose loop		
   } # end LET loop
  } # end if
} # end dataset (Z-init.energy) loop

write.table(simulation.results,"GSM_results.dat")