################################################################################################
# R script for recreating plots from Tsuruoka article
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

# load wrapping script
source("../../WrappingScripts/AmTrack_S.ssc")

# load data access script
source("data_operations.R")

# necessary libraries for plotting (R-specific)
library("lattice")
library("gplots")

############ useful functions ####################

AT.plot.surv.data <- function( SurvTable ){
   p1 <- xyplot( log10(Survival/100.0) ~ D.Gy | Z, groups = radiation, distribute.type = TRUE, data = SurvTable , type = c("p","l"))
   return (p1)
}

############ main script #########################

# read data into data frames
df.Survival = AT.survival.data.read()

cond  <- (df.Survival$radiation == "gamma")
gamma.data <- data.frame( D.Gy = df.Survival$D.Gy[ cond ], Survival = df.Survival$Survival[ cond ], Sur.min = df.Survival$Sur.min[ cond ] , Sur.max = df.Survival$Sur.max[ cond ] , type = "gamma")   

############# reorganize data ######################

# panel A data, carbon ions, 290MeV
carbon.290.let = c( 13, 19, 38, 54, 64, 73, 76 )
ion.data.carbon.290 <- list()
for( let in carbon.290.let ){
			cond  <- (df.Survival$Z == 6 & df.Survival$Init.Energy == 290 & df.Survival$LET == 10*let)
   ion.data.carbon.290[[let]] <- data.frame( D.Gy = df.Survival$D.Gy[ cond ], Survival = df.Survival$Survival[ cond ], Sur.min = df.Survival$Sur.min[ cond ] , Sur.max = df.Survival$Sur.max[ cond ] , type = "ion")
}

# panel B data, carbon ions, 135MeV
carbon.135.let = c( 38, 55, 84, 91, 94, 98 )
ion.data.carbon.135 <- list()
for( let in carbon.135.let ){
			cond  <- (df.Survival$Z == 6 & df.Survival$Init.Energy == 135 & df.Survival$LET == 10*let)
   ion.data.carbon.135[[let]] <- data.frame( D.Gy = df.Survival$D.Gy[ cond ], Survival = df.Survival$Survival[ cond ], Sur.min = df.Survival$Sur.min[ cond ] , Sur.max = df.Survival$Sur.max[ cond ] , type = "ion")
}

# panel C data, neon ions, 400MeV
neon.400.let = c( 45, 59, 77, 105, 132, 158, 177 )
ion.data.neon.400 <- list()
for( let in neon.400.let ){
			cond  <- (df.Survival$Z == 10 & df.Survival$Init.Energy == 400 & df.Survival$LET == 10*let)
   ion.data.neon.400[[let]] <- data.frame( D.Gy = df.Survival$D.Gy[ cond ], Survival = df.Survival$Survival[ cond ], Sur.min = df.Survival$Sur.min[ cond ] , Sur.max = df.Survival$Sur.max[ cond ] , type = "ion")
}

# panel D data, neon ions, 230MeV
neon.230.let = c( 30,44,58,77,105,127,156,184 )
ion.data.neon.230 <- list()
for( let in neon.230.let ){
			cond  <- (df.Survival$Z == 10 & df.Survival$Init.Energy == 230 & df.Survival$LET == 10*let)
   ion.data.neon.230[[let]] <- data.frame( D.Gy = df.Survival$D.Gy[ cond ], Survival = df.Survival$Survival[ cond ], Sur.min = df.Survival$Sur.min[ cond ] , Sur.max = df.Survival$Sur.max[ cond ] , type = "ion")
}

# panel E data, silicon ions, 490MeV
silicon.490.let = c( 55, 59, 69, 113, 145, 173, 214 )
ion.data.silicon.490 <- list()
for( let in silicon.490.let ){
			cond  <- (df.Survival$Z == 14 & df.Survival$Init.Energy == 490 & df.Survival$LET == 10*let)
   ion.data.silicon.490[[let]] <- data.frame( D.Gy = df.Survival$D.Gy[ cond ], Survival = df.Survival$Survival[ cond ], Sur.min = df.Survival$Sur.min[ cond ] , Sur.max = df.Survival$Sur.max[ cond ] , type = "ion")
}

# panel F data, iron ions, 500MeV
iron.500.let = c( 200, 260, 300, 400 )
ion.data.iron.500 <- list()
for( let in iron.500.let ){
			cond  <- (df.Survival$Z == 26 & df.Survival$Init.Energy == 500 & df.Survival$LET == 10*let)
   ion.data.iron.500[[let]] <- data.frame( D.Gy = df.Survival$D.Gy[ cond ], Survival = df.Survival$Survival[ cond ], Sur.min = df.Survival$Sur.min[ cond ] , Sur.max = df.Survival$Sur.max[ cond ] , type = "ion")
}

# plot data

############# produce plots and save to file ######################

pch.table <- c(22,15,1,16,2,17,6,18,23)

pdf("ExperimentalData_SurvivalCurves.pdf")

#plot1

# panel A plot, carbon ions, 290MeV
i <- 1
for( let in carbon.290.let ){
	df <- ion.data.carbon.290[[let]]
	cat("let = ", let , " pch = ", pch.table[i+1] ,"\n")
	plotCI( x = df$D.Gy , y = df$Survival/100., ui = df$Sur.max/100., li = df$Sur.min/100., gap = 0 , xlab = "" , ylab = "" , log = "y" , axes = FALSE, pch = ".", xlim = c(0,6), ylim = c(0.001, 1.) )
	points( x = df$D.Gy , y = df$Survival/100. , type = "b", pch = pch.table[i+1], bg = "black")
	par(new=TRUE)
	i <- (i+1)
}
df <- gamma.data
plotCI( x = df$D.Gy , y = df$Survival/100., ui = df$Sur.max/100., li = df$Sur.min/100., gap = 0 , xlab = "Dose (Gy)" , ylab = "Surviving fraction" , log = "y" , axes = FALSE, pch = 22, xlim = c(0,6), ylim = c(0.001, 1.))
legend( "topright", legend = c("200 kV X rays", paste( "LET = " , carbon.290.let, "keV/um" )), pch = pch.table)
lines( x = df$D.Gy , y = df$Survival/100. )
axis( 1, at = c(0,1,2,3,4,5,6))
axis( 2, at = c(0.001,0.01,0.1,1))
title( main = "Panel A, carbon ions, 290MeV/amu" )

# panel B plot, carbon ions, 135MeV
i <- 1
for( let in carbon.135.let ){
	df <- ion.data.carbon.135[[let]]
	cat("let = ", let , " pch = ", pch.table[i+1] ,"\n")
	plotCI( x = df$D.Gy , y = df$Survival/100., ui = df$Sur.max/100., li = df$Sur.min/100., gap = 0 , xlab = "" , ylab = "" , log = "y" , axes = FALSE, pch = ".", xlim = c(0,6), ylim = c(0.001, 1.) )
	points( x = df$D.Gy , y = df$Survival/100. , type = "b", pch = pch.table[i+1], bg = "black")
	par(new=TRUE)
	i <- (i+1)
}
df <- gamma.data
plotCI( x = df$D.Gy , y = df$Survival/100., ui = df$Sur.max/100., li = df$Sur.min/100., gap = 0 , xlab = "Dose (Gy)" , ylab = "Surviving fraction" , log = "y" , axes = FALSE, pch = 22, xlim = c(0,6), ylim = c(0.001, 1.) )
legend( "topright", legend = c("200 kV X rays", paste( "LET = " , carbon.135.let, "keV/um" )), pch = pch.table)
lines( x = df$D.Gy , y = df$Survival/100. )
axis( 1, at = c(0,1,2,3,4,5,6))
axis( 2, at = c(0.001,0.01,0.1,1))
title( main = "Panel B, carbon ions, 135MeV/amu" )

# panel C plot, neon ions, 400MeV
i <- 1
for( let in neon.400.let ){
	df <- ion.data.neon.400[[let]]
	cat("let = ", let , " pch = ", pch.table[i+1] ,"\n")
	plotCI( x = df$D.Gy , y = df$Survival/100., ui = df$Sur.max/100., li = df$Sur.min/100., gap = 0 , xlab = "" , ylab = "" , log = "y" , axes = FALSE, pch = ".", xlim = c(0,6), ylim = c(0.001, 1.) )
	points( x = df$D.Gy , y = df$Survival/100. , type = "b", pch = pch.table[i+1], bg = "black")
	par(new=TRUE)
	i <- (i+1)
}
df <- gamma.data
plotCI( x = df$D.Gy , y = df$Survival/100., ui = df$Sur.max/100., li = df$Sur.min/100., gap = 0 , xlab = "Dose (Gy)" , ylab = "Surviving fraction" , log = "y" , axes = FALSE, pch = 22, xlim = c(0,6), ylim = c(0.001, 1.) )
legend( "topright", legend = c("200 kV X rays", paste( "LET = " , neon.400.let, "keV/um" )), pch = pch.table)
lines( x = df$D.Gy , y = df$Survival/100. )
axis( 1, at = c(0,1,2,3,4,5,6))
axis( 2, at = c(0.001,0.01,0.1,1))
title( main = "Panel C, neon ions, 400MeV/amu" )

# panel D plot, neon ions, 230MeV
i <- 1
for( let in neon.230.let ){
	df <- ion.data.neon.230[[let]]
	cat("let = ", let , " pch = ", pch.table[i+1] ,"\n")
	plotCI( x = df$D.Gy , y = df$Survival/100., ui = df$Sur.max/100., li = df$Sur.min/100., gap = 0 , xlab = "" , ylab = "" , log = "y" , axes = FALSE, pch = ".", xlim = c(0,6), ylim = c(0.001, 1.) )
	points( x = df$D.Gy , y = df$Survival/100. , type = "b", pch = pch.table[i+1], bg = "black")
	par(new=TRUE)
	i <- (i+1)
}
df <- gamma.data
plotCI( x = df$D.Gy , y = df$Survival/100., ui = df$Sur.max/100., li = df$Sur.min/100., gap = 0 , xlab = "Dose (Gy)" , ylab = "Surviving fraction" , log = "y" , axes = FALSE, pch = 22, xlim = c(0,6), ylim = c(0.001, 1.) )
legend( "topright", legend = c("200 kV X rays", paste( "LET = " , neon.230.let, "keV/um" )), pch = pch.table)
lines( x = df$D.Gy , y = df$Survival/100. )
axis( 1, at = c(0,1,2,3,4,5,6))
axis( 2, at = c(0.001,0.01,0.1,1))
title( main = "Panel D, neon ions, 230MeV/amu" )

# panel E plot, silicon ions, 490MeV
i <- 1
for( let in silicon.490.let ){
	df <- ion.data.silicon.490[[let]]
	cat("let = ", let , " pch = ", pch.table[i+1] ,"\n")
	plotCI( x = df$D.Gy , y = df$Survival/100., ui = df$Sur.max/100., li = df$Sur.min/100., gap = 0 , xlab = "" , ylab = "" , log = "y" , axes = FALSE, pch = ".", xlim = c(0,6), ylim = c(0.001, 1.) )
	points( x = df$D.Gy , y = df$Survival/100. , type = "b", pch = pch.table[i+1], bg = "black")
	par(new=TRUE)
	i <- (i+1)
}
df <- gamma.data
plotCI( x = df$D.Gy , y = df$Survival/100., ui = df$Sur.max/100., li = df$Sur.min/100., gap = 0 , xlab = "Dose (Gy)" , ylab = "Surviving fraction" , log = "y" , axes = FALSE, pch = 22, xlim = c(0,6), ylim = c(0.001, 1.) )
legend( "topright", legend = c("200 kV X rays", paste( "LET = " , silicon.490.let, "keV/um" )), pch = pch.table)
lines( x = df$D.Gy , y = df$Survival/100. )
axis( 1, at = c(0,1,2,3,4,5,6))
axis( 2, at = c(0.001,0.01,0.1,1))
title( main = "Panel E, silicon ions, 290MeV/amu" )

# panel F plot, iron ions, 500MeV
i <- 1
for( let in iron.500.let ){
	df <- ion.data.iron.500[[let]]
	cat("let = ", let , " pch = ", pch.table[i+1] ,"\n")
	plotCI( x = df$D.Gy , y = df$Survival/100., ui = df$Sur.max/100., li = df$Sur.min/100., gap = 0 , xlab = "" , ylab = "" , log = "y" , axes = FALSE, pch = ".", xlim = c(0,6), ylim = c(0.001, 1.) )
	points( x = df$D.Gy , y = df$Survival/100. , type = "b", pch = pch.table[i+1], bg = "black")
	par(new=TRUE)
	i <- (i+1)
}
df <- gamma.data
plotCI( x = df$D.Gy , y = df$Survival/100., ui = df$Sur.max/100., li = df$Sur.min/100., gap = 0 , xlab = "Dose (Gy)" , ylab = "Surviving fraction" , log = "y" , axes = FALSE, pch = 22, xlim = c(0,6), ylim = c(0.001, 1.) )
legend( "topright", legend = c("200 kV X rays", paste( "LET = " , iron.500.let, "keV/um" )), pch = pch.table)
lines( x = df$D.Gy , y = df$Survival/100. )
axis( 1, at = c(0,1,2,3,4,5,6))
axis( 2, at = c(0.001,0.01,0.1,1))
title( main = "Panel F, iron ions, 500MeV/amu" )

dev.off()
