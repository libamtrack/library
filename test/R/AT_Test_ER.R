# clear workspace
rm( list = ls() )

# load libAmTrack library
dyn.load("../../Release/libAmTrack.dll")

# load wrapping scripts
source("../../wrapper/R/AmTrack.R")

# necessary library for plotting
library("lattice")

# energy range definitions:
#E.MeV.u <- c( 1, 10, 60, 100, 250, 1000)

expn <- seq (-3, 4, by=1e-1)
E.MeV.u <- 10^expn

#E.MeV.u <- seq( 1, 200, by=1)

# models definition
#er.models <- c(2,3,4,5,7)
er.models.names <- c("simple test ER model",  "Butts & Katz' ER model (linear)",  "Waligorski's ER model (power-law wmax)",  "Geiss' ER model (power-law E)", "Scholz' ER model (power-law E)", "Edmund' ER model (power-law wmax)","Tabata  ER model")

er.models <- c(1,2,3,4,5,6,7)

# other parameters
n <- length(E.MeV.u) 
material.number <- rep(1, 1) # Water

# data frame setup
df <- expand.grid( E.MeV.u = E.MeV.u , er.models = er.models )

df$er.models.name	<- as.character(er.models.names[df$er.models])
df$range.m	<- numeric(nrow(df))
df$wmax <- 0

# calculations
j <- 0
for( i in er.models ){
 j <- j+1
 ii					<-	df$er.models == i
# cat("i = ",i,"\n")
 df$er.models.name[ii] <- er.models.names[j]
# beta <- AT.beta.from.E( df$E.MeV.u[ii] )
# cat("E_MeV = ",df$E.MeV.u[ii],"\n")
  wmax_MeV <- AT.max.E.transfer.MeV( E.MeV.u = df$E.MeV.u[ii] )
# cat("wmax_MeV = ",wmax_MeV,"\n")
 df$wmax[ii] <- wmax_MeV
	df$range.m[ii] <- 	AT.max.electron.range(	E.MeV.u = df$E.MeV.u[ii], material.number, i)
}

# plots...
#logplot <- xyplot( 1e2*range.m ~ E.MeV.u, groups = er.models.name, ref = TRUE, data=df, pch = ".", lty = 1, type = "l", xlab = "Energy [MeV]", ylab = "Range [cm]", auto.key = list(title = "Range of delta electrons in liquid water",points = FALSE, lines = TRUE), scales = list(log = 10))
#linplot <- xyplot( 1e2*range.m ~ E.MeV.u, groups = er.models.name, ref = TRUE, data=df, pch = ".", lty = 1, type = "l", xlab = "Energy [MeV]", ylab = "Range [cm]", auto.key = list(title = "Range of delta electrons in liquid water",points = FALSE, lines = TRUE))

logplot <- xyplot( 1e2*range.m ~ wmax, groups = er.models.name, ref = TRUE, data=df, pch = ".", lty = 1, type = "l", xlab = "wmax [MeV]", ylab = "Range [cm]", auto.key = list(title = "Range of delta electrons in liquid water",points = FALSE, lines = TRUE), scales = list(log = 10))
linplot <- xyplot( 1e2*range.m ~ wmax, groups = er.models.name, ref = TRUE, data=df, pch = ".", lty = 1, type = "l", xlab = "wmax [MeV]", ylab = "Range [cm]", auto.key = list(title = "Range of delta electrons in liquid water",points = FALSE, lines = TRUE))


pdf("ER.pdf")

logplot
linplot

dev.off()


