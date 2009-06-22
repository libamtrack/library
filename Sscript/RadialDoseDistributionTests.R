rm( list = ls() )
dyn.load("/home/grzanka/workspace/DKFZ/Steffen/SGParticleLibrary/Release/libSGPL.so")
source("SGP.ssc")
library("lattice")

pdf("/home/grzanka/workspace/DKFZ/Steffen/SGParticleLibrary/SGPL-R/RDD.pdf")

####################### DOSE(DISTANCE) function test #########################################

# radius definitions:
r.expn <- seq (-12, -7, by=0.1)
r.m <- 10^r.expn 
#r.m <- c(1e-7,5e-8)

# energy range definitions:
#E.MeV.u <- c( 1, 10, 60, 100, 250, 1000)
#E.MeV.u <- c( 60, 250)
# other parameters
#nn <- length(E.MeV.u) 
#particle.index <- rep(1, nn)
#material.name <- "Water, Liquid"

# electron range models definition
#ER.model <- c(1,2,3,4)
ER.model.names <- c("KatzScholz","Waligorski, 1986","Geiss, 1998","Scholz, 2001")
ER.model <- c(2)

# RDD models definition
#RDD.model <- c(1,2,3)
RDD.model.names <- c("Test","KatzPoint","Geiss","Site")
RDD.model <- c(0,1,2,3)

# RDD parameters
RDD.parameters <- list(c(100,1,0,1e-10),c(100, 1, 0, 1e-12), c(100,1,0,1e-10),c(100,1,0,1e-10))

# data frame setup
df1 <- expand.grid( r.m = r.m, ER.model = ER.model, RDD.model = RDD.model )

df1$ER.model.name	<- as.character(ER.model.names[df1$ER.model])
df1$RDD.model.name	<- as.character(RDD.model.names[df1$RDD.model+1])

df1$D.Gy	<- numeric(nrow(df1))

for( i in ER.model) {
 ii 			<- df1$ER.model == i
	for( j in RDD.model) {
 	jj 			<- ((df1$RDD.model == j) & (df1$ER.model == i)) 
		df1$D.Gy[jj] <- SGP.RDD.D.Gy(r.m = df1$r.m[jj], ER.model = i, RDD.model = j, RDD.parameters=RDD.parameters[[j+1]] ) 	
		}
}

df1

# plots...

p1 <- xyplot( D.Gy ~ r.m | ER.model.name , groups = RDD.model.name, ref = TRUE, data=df1, pch = ".", lty = 1, type = "l", xlab = "Distance [m]", ylab = "Dose [Gy]", auto.key = list(title = "Protons in liquid water, RDD",points = FALSE, lines = TRUE), scales = list(log = 10))
#p2 <- xyplot( D.Gy ~ r.m , groups = RDD.model.name, ref = TRUE, data=df1, pch = ".", lty = 1, type = "l", xlab = "Distance [m]", ylab = "Dose [Gy]", auto.key = list(title = "Protons in liquid water, RDD",points = FALSE, lines = TRUE), scales = list(log = 10))

p1
#p2

####################### DISTANCE(DOSE) function test #########################################

# energy range definitions:
#D.Gy <- c( 1, 10, 60, 100, 250, 1000)
#D.Gy <- c( 7.115367e-02, 7.116197e+02)
D.expn <- seq (-2, 7, by=0.1)
D.Gy <- 10^D.expn
#D.Gy <- seq( 1, 200, by=1)

ER.model <- c(2)
RDD.model <- c(1,2,3)

df2 <- expand.grid( D.Gy = D.Gy , ER.model = ER.model, RDD.model = RDD.model )

df2$ER.model.name	<- as.character(ER.model.names[df2$ER.model])
df2$RDD.model.name	<- as.character(RDD.model.names[df2$RDD.model+1])
df2$r.m	<- numeric(nrow(df2))

for( i in ER.model) {
 ii 			<- df2$ER.model == i
for( j in RDD.model) {
 	jj 			<- ((df2$RDD.model == j) & (df2$ER.model == i)) 
		df2$r.m[jj] <- SGP.RDD.r.m(D.Gy = df2$D.Gy[jj], ER.model = i, RDD.model = j, RDD.parameters=RDD.parameters[[j]] ) 	
		}
}

df2

# plots...

p3 <- xyplot( r.m ~ D.Gy | ER.model.name , groups = RDD.model.name, ref = TRUE, data=df2, pch = ".", lty = 1, type = "l", xlab = "Dose [Gy]", ylab = "Distance [m]", auto.key = list(title = "Protons in liquid water, inverse RDD",points = FALSE, lines = TRUE), scales = list(log = 10))
#p4 <- xyplot( r.m ~ D.Gy, groups = RDD.model.name, ref = TRUE, data=df2, pch = ".", lty = 1, type = "l", xlab = "Dose [Gy]", ylab = "Distance [m]", auto.key = list(title = "Protons in liquid water, inverse RDD",points = FALSE, lines = TRUE), scales = list(log = 10))

p3
#p4

dev.off()