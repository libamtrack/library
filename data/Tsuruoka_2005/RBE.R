# clear workspace
rm( list = ls() )

# load libAmTrack library
dyn.load("../../Release/libAmTrack.dll")

# load wrapping scripts
source("../../WrappingScripts/SGP.ssc")

library("gplots")
library("lattice")

############ read data from files #########################

SGP.survival.data.read	<-	function(	){

# read gamma data
	SurvTable <- read.table("Tsuruoka_gamma.dat")
 SurvTable$radiation = "gamma"
	SurvTable$Z = 0
	SurvTable$beta = 0
	SurvTable$zzbb = 0
	SurvTable$LET = 0
	SurvTable$Init.Energy = 0
	
# read ion data
	tmp.SurvTable <- read.table("Tsuruoka_ion.dat")
	tmp.SurvTable$radiation = "ion"
	SurvTable <- rbind( SurvTable , tmp.SurvTable )

 return (SurvTable)
}

############ constants #########################

SGP.survival.data.names <- function ( Z ){
 if( Z == 1 ){
 	return ("proton")
 } else if( Z == 6 ){
 	return ("carbon")
 } else if( Z == 10 ){
 	return ("neon")
 } else if( Z == 14 ){
 	return ("silicon")
 } else if( Z == 26 ){
 	return ("iron")
 } else {
 	return (Z)
 }

}

# read data into data frames
df.Survival = SGP.survival.data.read()

cond  <- (df.Survival$radiation == "gamma")
gamma.data <- data.frame( D.Gy = df.Survival$D.Gy[ cond ], Survival = df.Survival$Survival[ cond ], Sur.min = df.Survival$Sur.min[ cond ] , Sur.max = df.Survival$Sur.max[ cond ] , type = "gamma")   
gamma.params <- SGP.fit.linquad( gamma.data$D.Gy, gamma.data$Survival, gamma.data$Sur.min, gamma.data$Sur.max, 0., 0. )

cond <- paste( df.Survival$Z , df.Survival$Init.Energy, df.Survival$LET, sep="-")
df.Survival$cond <- cond

RBE.cond <- paste( df.Survival$Z[df.Survival$Z > 0] , df.Survival$Init.Energy[df.Survival$Z > 0], sep="-")
df.Summary <- data.frame( Z = df.Survival$Z[df.Survival$Z > 0], LET = df.Survival$LET[df.Survival$Z > 0], Init.Energy = df.Survival$Init.Energy[df.Survival$Z > 0], RBE = 0) 
df.Summary$cond <- RBE.cond

df.Summary <- unique(df.Summary)

RBE.vector <- numeric(0)

pdf("LQfit.pdf")

for( item in unique(df.Survival$cond) ){
 if( item != "0-0-0" ){
  cond <- ( df.Survival$cond == item )
  ion.params <- c(0.,0.)
  try(ion.params <- SGP.fit.linquad( df.Survival$D.Gy[ cond ], df.Survival$Survival[ cond ], df.Survival$Sur.min[ cond ] , df.Survival$Sur.max[ cond ], 0., 0. ))
	 df <- gamma.data
	 plotCI( x = df$D.Gy , y = df$Survival/100., ui = df$Sur.max/100., li = df$Sur.min/100., gap = 0 , xlab = "Dose (Gy)" , ylab = "Surviving fraction" , log = "y" , axes = FALSE, pch = 22, xlim = c(0,6), ylim = c(0.001, 1.))
  x <- seq( 0, 6, by=0.01)
  y <- exp( -gamma.params[1]*x - gamma.params[2]*x*x)
  lines( x = x , y = y )
	 df <- data.frame( D.Gy = df.Survival$D.Gy[ cond ], Survival = df.Survival$Survival[ cond ], Sur.min = df.Survival$Sur.min[ cond ] , Sur.max = df.Survival$Sur.max[ cond ])
	 par(new=TRUE)
	 plotCI( x = df$D.Gy , y = df$Survival/100., ui = df$Sur.max/100., li = df$Sur.min/100., gap = 0 , xlab = "Dose (Gy)" , ylab = "Surviving fraction" , log = "y" , axes = FALSE, pch = 22, xlim = c(0,6), ylim = c(0.001, 1.))
  legend( "topright", legend = c("200 kV X rays", paste( "LET = " , unique(df.Survival$LET[ cond ]) , "keV/um" ) ))
  x <- seq( 0, 6, by=0.01)
  y <- exp( -ion.params[1]*x - ion.params[2]*x*x)
  lines( x = x , y = y )
  abline( h = c(0.1) )
  axis( 1, at = c(0,1,2,3,4,5,6))
  axis( 2, at = c(0.001,0.01,0.1,1))
  RBE <- SGP.RBE( ion.params[1],  ion.params[2], gamma.params[1], gamma.params[2] )
  RBE.vector <- c( RBE.vector , RBE)
  title( main = paste("Z = ", unique(df.Survival$Z[ cond ]), " init energy = ", unique(df.Survival$Init.Energy[ cond ]) , "MeV/amu\n ion.alpha = ", round(ion.params[1],3), "ion.beta = ", round(ion.params[2],3), "\n RBE = ", round(RBE,3)) )
}	
} 
dev.off()

df.Summary$RBE <- RBE.vector 

pdf("RBE.pdf")

#plot( df.Summary$RBE ~ df.Summary$LET, groups = df.Summary$Z, xlab = "LET [keV/um]", ylab = "RBE" )
p1 <- xyplot( RBE ~ LET | paste("Z =",Z,lapply( df.Summary$Z, SGP.survival.data.names)), data = df.Summary , groups = cond, auto.key = list(), type = "p" , scales = c( list( x = list(log = 10)) ))
p1

dev.off()