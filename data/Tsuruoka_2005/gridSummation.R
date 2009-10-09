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

# panel A data, carbon ions, 290MeV
#carbon.290.let = c( 13, 19, 38, 54, 64, 73, 76 )
carbon.290.let = c( 13 )
ion.data.carbon.290 <- list()
ion.model.data.carbon.290 <- list()
ion.model.params.carbon.290 <- list()
ion.gamma.data.carbon.290 <- list()

dose <- seq(0,6, by=1)
n <- length(dose)

for( let in carbon.290.let ){
  cond  <- (df.Survival$Z == 6 & df.Survival$Init.Energy == 290 & df.Survival$LET == 10*let)
  ion.data.carbon.290[[let]] <- data.frame( D.Gy = df.Survival$D.Gy[ cond ], Survival = df.Survival$Survival[ cond ], Sur.min = df.Survival$Sur.min[ cond ] , Sur.max = df.Survival$Sur.max[ cond ] , type = "ion")

  LET.MeV.cm2.g = 10*let
  particle.no <- 18 # carbon ion
  fluence.cm2 <- -1.0
  material.no <- 1 # water
  E.MeV.u <- AT.E.MeV.u( LET.MeV.cm2.g , particle.no, material.no )[1]
  cat("E MeV", E.MeV.u, "\n")
  RDD.model <- 3
  RDD.parameters <- 1e-8
  ER.model <- 2
  ER.parameters <- 1
  gamma.model <- 5
  gamma.parameters <- c(gamma.LQ.params,20.0)

  HCP = numeric()
  gamma = numeric()
  dose.factor = numeric()
   
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
      N.runs = 50,
      N2 = 10,          
      fluence.factor = 1.0,
      write.output = F,
      n.X = 100,
      lethal.events.mode = T,
      grid.size.m = 1e-7)
 
    HCP = c(HCP, 100*res[3])
    gamma = c(gamma,100*res[4])
    dose.factor = c(dose.factor, res[2]/d)

    ion.model.data.carbon.290[[let]] <- HCP
    ion.model.params.carbon.290[[let]] <- AT.fit.linquad( dose, HCP, HCP , HCP, 0., 0. )
        
    ion.gamma.data.carbon.290[[let]] <- gamma
   }
   
   ion.model.df.carbon.290 <- data.frame( D.Gy = dose, Survival = HCP , dose.factor = dose.factor)
   write.table(ion.model.df.carbon.290,paste("grid_carbon290_let",let,".txt", sep=""))
 
}
