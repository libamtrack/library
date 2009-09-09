rm( list = ls() )
dyn.load("../Release/libAmTrack.dll")
source("../WrappingScripts/SGP.ssc")
library("lattice")

# energy range definitions:
#E.MeV.u <- c( 1, 10, 60, 100, 250, 1000)

expn <- seq (0, 3, by=1e-1)
E.MeV.u <- 10^expn

#E.MeV.u <- seq( 1, 200, by=1)

# models definition
er.models <- c(1,2,3,4)
er.models.names <- c("KatzScholz","Waligorski et al. 1986","Geiss, 1998","Scholz, 2001")

# other parameters
n <- length(E.MeV.u) 
particle.number <- rep(1, n)
material.number <- rep(1, 1) # Water

# data frame setup
df <- expand.grid( E.MeV.u = E.MeV.u , er.models = er.models )

df$er.models.name	<- as.character(er.models.names[df$er.models])
df$range.m	<- numeric(nrow(df))

# calculations
for( i in er.models ){
 ii					<-	df$er.models == i
	df$range.m[ii] <- 	SGP.max.electron.range(	E.MeV.u = df$E.MeV.u[ii], particle.number, material.number, i)
}

df

# plots...
logplot <- xyplot( 1e2*range.m ~ E.MeV.u, groups = er.models.name, ref = TRUE, data=df, pch = ".", lty = 1, type = "l", xlab = "Energy [MeV]", ylab = "Range [cm]", auto.key = list(title = "Protons in liquid water, range of delta electrons",points = FALSE, lines = TRUE), scales = list(log = 10))
linplot <- xyplot( 1e2*range.m ~ E.MeV.u, groups = er.models.name, ref = TRUE, data=df, pch = ".", lty = 1, type = "l", xlab = "Energy [MeV]", ylab = "Range [cm]", auto.key = list(title = "Protons in liquid water, range of delta electrons",points = FALSE, lines = TRUE))

pdf("ER.pdf")

logplot
linplot

dev.off()