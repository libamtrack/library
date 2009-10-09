# clear workspace
rm( list = ls() )

# load libAmTrack library
dyn.load("../../Release/libAmTrack.dll")

# load wrapping scripts
source("../../WrappingScripts/AmTrack_S.ssc")

# necessary libraries for plotting
library("lattice")
library("gplots")

############ read data from files #########################

AT.survival.data.read <- function( ){

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

############ main script #########################

# read data into data frames
df.Survival = AT.survival.data.read()

cond  <- (df.Survival$radiation == "gamma")
gamma.data <- data.frame( D.Gy = df.Survival$D.Gy[ cond ], Survival = df.Survival$Survival[ cond ], Sur.min = df.Survival$Sur.min[ cond ] , Sur.max = df.Survival$Sur.max[ cond ] , type = "gamma")   
gamma.LQ.params <- AT.fit.linquad( gamma.data$D.Gy, gamma.data$Survival, gamma.data$Sur.min, gamma.data$Sur.max, 0., 0. )

# plotting
pdf("Grid.pdf", width=15, height=15, pointsize=24)

# panel A data, carbon ions, 290MeV
# carbon.290.let = c( 13, 19, 38, 54, 64, 73, 76 )
carbon.290.let = c( 13 )
ion.data.carbon.290 <- list()
ion.model.data.carbon.290 <- list()
ion.model.params.carbon.290 <- list()
ion.gamma.data.carbon.290 <- list()

Katz.chi2 <- numeric()
LEM.chi2 <- numeric()

v.RBE.exp <- numeric()
v.RBE.Katz <- numeric()
v.RBE.LEM <- numeric()

for( let in carbon.290.let ){
  cat("****************** LET = ",let," [keV/um] ***************\n")

  exp.dose <- df.Survival$D.Gy[ cond ]
  exp.surv <- df.Survival$Survival[ cond ]/100.
  exp.err <- df.Survival$Sur.Err[ cond ]/100.

  ion.model.df.carbon.290 <- read.table(paste("grid_carbon290_let",let,".txt", sep=""))
  dose <- ion.model.df.carbon.290$D.Gy
  surv <- ion.model.df.carbon.290$Survival
  dose.factor <- ion.model.df.carbon.290$dose.factor

# plot of gamma data
  df <- gamma.data
  plotCI( x = df$D.Gy , y = df$Survival/100., ui = df$Sur.max/100., li = df$Sur.min/100., gap = 0 , xlab = "Dose (Gy)" , ylab = "Surviving fraction" , log = "y" , axis = FALSE, pch = 22, xlim = c(0,6), ylim = c(0.001, 1.), main=paste("Carbon ions (LET=",let," keV/um) and photons, cell survival") , lw = 4, cex.lab = 1.2)
  
 # plot of LQ fit to gamma data
  x <- seq( 0, 6, by=0.01)
  y <- exp( -gamma.LQ.params[1]*x - gamma.LQ.params[2]*x*x)
#  lines( x = x , y = y )
    
  # plot of LQ fit to simulated ion data by grid summation
  ion.params <- c(0.,0.)
  ion.params <- AT.fit.linquad( dose, surv, surv , surv, 0., 0. )
  x <- seq( 0, 6, by=0.01)
  y <- exp( -ion.params[1]*x - ion.params[2]*x*x)
  lines( x = x , y = y , col = "red", lw = 4)
  mod.L.surv <- exp( -ion.params[1]*exp.dose - ion.params[2]*exp.dose*exp.dose)
  RBE.LEM <- AT.RBE( ion.params[1],  ion.params[2], gamma.LQ.params[1], gamma.LQ.params[2] )
  cat("   LEM.ion.params",ion.params,"\n")
  cat("RBE LEM",RBE.LEM,"\n")
  v.RBE.LEM <- c( v.RBE.LEM, RBE.LEM )

  # plot of experimental ion data 
  cond  <- (df.Survival$Z == 6 & df.Survival$Init.Energy == 290 & df.Survival$LET == 10*let)
  ion.data.carbon.290 <- data.frame( D.Gy = df.Survival$D.Gy[ cond ], Survival = df.Survival$Survival[ cond ], Sur.min = df.Survival$Sur.min[ cond ] , Sur.max = df.Survival$Sur.max[ cond ] , type = "ion")
  par(new=TRUE)
  plotCI( x = ion.data.carbon.290$D.Gy , y = ion.data.carbon.290$Survival/100., ui = ion.data.carbon.290$Sur.max/100., li = ion.data.carbon.290$Sur.min/100., gap = 0 , xlab = "Dose (Gy)" , ylab = "Surviving fraction" , log = "y" , axis = FALSE, pch = 2, xlim = c(0,6), ylim = c(0.001, 1.), lw = 4, cex.lab = 1.2)
  
  # plot of LQ fit to experimental ion data
  ion.params <- c(0.,0.)
  try(ion.params <- AT.fit.linquad( df.Survival$D.Gy[ cond ], df.Survival$Survival[ cond ], df.Survival$Sur.min[ cond ] , df.Survival$Sur.max[ cond ], 0., 0. ))
  x <- seq( 0, 6, by=0.01)
  y <- exp( -ion.params[1]*x - ion.params[2]*x*x)
#  lines( x = x , y = y )
  abline( h = c(0.1) ,lw = 2)
  axis( 1, at = c(0,1,2,3,4,5,6), lw = 4)
  axis( 2, at = c(0.001,0.01,0.1,1), lw = 4)

  RBE.exp <- AT.RBE( ion.params[1],  ion.params[2], gamma.LQ.params[1], gamma.LQ.params[2] )
  cat("   ion.params",ion.params,"\n")
  cat("RBE exp",RBE.exp,"\n")
  v.RBE.exp <- c( v.RBE.exp, RBE.exp )

  cat("Reading file... ", paste("model_let",let,".dat", sep=""))
  katz.model.df <- read.table(paste("model_let",let,".dat", sep=""))
#  points( x = katz.model.df$V1 , y = katz.model.df$V2/100. , pch="*")
  ion.params <- c(0.,0.)
  try(ion.params <- AT.fit.linquad( katz.model.df$V1, katz.model.df$V2, katz.model.df$V2 , katz.model.df$V2, 0., 0. ))
  x <- seq( 0, 6, by=0.01)
  y <- exp( -ion.params[1]*x - ion.params[2]*x*x)
  lines( x = x , y = y , col = "green", lw = 4)
  mod.K.surv <- exp( -ion.params[1]*exp.dose - ion.params[2]*exp.dose*exp.dose)
  RBE.Katz <- AT.RBE( ion.params[1],  ion.params[2], gamma.LQ.params[1], gamma.LQ.params[2] )
  cat("   Katz.ion.params",ion.params,"\n")
  cat("RBE Katz",RBE.Katz,"\n")
  v.RBE.Katz <- c( v.RBE.Katz, RBE.Katz )

 # cat("LET = ",let,"\n")
 # cat("A",exp.err,"\n")
 # cat("B",exp.surv,"\n")
 # cat("C",mod.surv,"\n")
 # cat("length",length(exp.dose),"\n")
  cat("Katz LET = ",let,"chi2", sum((exp.surv-mod.K.surv)*(exp.surv-mod.K.surv)/(exp.err*exp.err))/length(exp.dose),"\n")
  cat("LEM LET = ",let,"chi2", sum((exp.surv-mod.L.surv)*(exp.surv-mod.L.surv)/(exp.err*exp.err))/length(exp.dose),"\n")

  kc <- sum((exp.surv-mod.K.surv)*(exp.surv-mod.K.surv)/(exp.err*exp.err))/length(exp.dose)
  Katz.chi2 <- c(Katz.chi2, kc )
  lc <- sum((exp.surv-mod.L.surv)*(exp.surv-mod.L.surv)/(exp.err*exp.err))/length(exp.dose)
  LEM.chi2 <- c(LEM.chi2, lc )

  legend( "topright", legend = c("exp. data: photons", "exp. data: carbon ions", paste( "Grid summation model" ),paste( "Katz model" ) ), pch = c(22,2,NA,NA), lty = c(NA,NA,1,1), lw = 4, col=c(1,1,"red","green"))
  #par(new=TRUE)
  #plot(x = NA, main=paste("Carbon beam (LET=",let," keV/um) and X rays, cell survival"), log = "y" , axes = FALSE, pch = 3, xlim = c(0,6), ylim = c(0.001, 1.))
}

plot(Katz.chi2 ~ carbon.290.let, xlim = c(0,100), ylim = c(0,500), pch = 1, xlab = "LET", ylab = "chi2/n")
points(LEM.chi2 ~ carbon.290.let, pch = 2)
#legend( "topright", legend = c("Katz","LEM"), pch = c(1,2))

dev.off()

pdf("RBE-comparison.pdf")

letu <- c(0,0,1,3,5,7,8)
letl <- c(0,0,1,3,6,11,13)
#p1 <- plotCI( x = carbon.290.let, y = v.RBE.exp, uiw = letu, liw = letl, err = "x", xlim = c(0,80), ylim = c(0,4), xlab = "LET [keV/um]", ylab = "RBE" , main = "RBE: models vs experiment", lw = 4)
p1 <- plotCI( x = carbon.290.let, y = v.RBE.exp, uiw = letu, liw = letl, err = "x", xlim = c(0,80), ylim = c(0,4), xlab = substitute(paste("LET [keV/",mu,"m]")), ylab = "RBE" , lw = 4, cex.lab = 1.5)
par(new=TRUE)
rbeu <- c(0.01,0.05,0.01,0.08,0.02,0.16,0.3)
rbel <- c(0.01,0.05,0.01,0.08,0.02,0.16,0.3)
plotCI( x = carbon.290.let, y = v.RBE.exp, uiw = rbeu, liw = rbel, err = "y", xlim = c(0,80), ylim = c(0,4), xlab = NA, ylab = NA, axis = FALSE, lw = 4, add = TRUE)
lines( v.RBE.LEM ~ carbon.290.let, col="red", lw = 4)
lines( v.RBE.Katz ~ carbon.290.let, col="green", lw = 4)
legend( "topright", legend = c("experimental data","Katz model","Grid summ. model"), pch=1, col=c("black","green","red"), lw = 4)
p1

dev.off()
