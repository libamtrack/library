# clear workspace
rm( list = ls() )

# load libAmTrack library
dyn.load("../Release/libAmTrack.dll")

# load wrapping scripts
source("../WrappingScripts/SGP.ssc")

# necessary library for plotting
library("lattice")

####################### DOSE(DISTANCE) function test #########################################

# radius definitions:
r.m <- 10^seq (-14, -3.5, length.out=1000)
#r.m <- seq( 1e-9,1.5*1e-8, length.out=10)
#r.m <- c(1e-6, 1e-7, 1e-8, 1e-9, 1e-12)

# energy range definitions:
#E.MeV.u <- c( 1, 10, 60, 100, 250, 1000)
E.MeV.u <- c( 100)
# other parameters
nn <- length(E.MeV.u) 
particle.no <- rep(1, 2*nn)
material.no <- c(1)

# electron range models definition
#ER.model <- c(1,2,3,4)
ER.model.names <- c("KatzScholz","Waligorski, 1986","Geiss, 1998","Scholz, 2001")
ER.model <- c(2)

# RDD models definition
RDD.model <- c(1,2,3,4,5)
RDD.model.names <- c("Test","KatzPoint","Geiss","Site","Extended Target")
#RDD.model <- c(0,1,2,3,4)

nn <- length(E.MeV.u)
nnn <-  length(RDD.model)
particle.no <- rep(1, nnn*nn)
material.no <- c(1)

# RDD parameters
RDD.parameters <- list(c(1e-10),c(1e-10,1e-10), c(1e-10),c(1e-8,1e-10),c(1e-10,1e-8,1e-10))

# data frame setup
df1 <- expand.grid( r.m = r.m, ER.model = ER.model, RDD.model = RDD.model )

df1$ER.model.name	<- as.character(ER.model.names[df1$ER.model])
df1$RDD.model.name	<- as.character(RDD.model.names[df1$RDD.model])

df1$D.Gy	<- numeric(nrow(df1))

E.MeV.u = rep( 100, nrow(df1))
particle.no = rep( 1, nrow(df1))

ER.parameters <- c(0)

for( i in ER.model) {
 ii 			<- df1$ER.model == i
	for( j in RDD.model) {
 	jj 			<- ((df1$RDD.model == j) & (df1$ER.model == i)) 
			df1$D.Gy[jj] <- SGP.RDD.D.Gy(r.m = df1$r.m[jj], E.MeV.u, particle.no, material.no, ER.model = i, ER.parameters = ER.parameters, RDD.model = j, RDD.parameters=RDD.parameters[[j]] ) 	
		}
}

#df1

# plots...

p1 <- xyplot( D.Gy ~ r.m | ER.model.name , groups = RDD.model.name, ref = TRUE, data=df1, pch = ".", lty = 1, type = "l", xlab = "Distance [m]", ylab = "Dose [Gy]", auto.key = list(title = "Protons in liquid water, RDD",points = FALSE, lines = TRUE), scales = list(log = 10))
p2 <- xyplot( D.Gy ~ r.m | ER.model.name , groups = RDD.model.name, ref = TRUE, data=df1, pch = ".", lty = 1, type = "l", xlab = "Distance [m]", ylab = "Dose [Gy]", auto.key = list(title = "Protons in liquid water, RDD",points = FALSE, lines = TRUE))

p1
p2

####################### DISTANCE(DOSE) function test #########################################

# energy range definitions:
#D.Gy <- c( 1e-8,1e-6,1e-4)
#D.Gy <- c( 7.115367e-02, 7.116197e+02)
#D.expn <- seq (-2, 7, by=0.1)
D.Gy <- 10^seq(-10,-1,length.out=100)

ER.model <- c(2)
RDD.model <- c(2,4,5)

df2 <- expand.grid( D.Gy = D.Gy , ER.model = ER.model, RDD.model = RDD.model )

df2$ER.model.name	<- as.character(ER.model.names[df2$ER.model])
df2$RDD.model.name	<- as.character(RDD.model.names[df2$RDD.model])
df2$r.m	<- numeric(nrow(df2))

E.MeV.u = rep( 100, nrow(df2))
particle.no = rep( 1, nrow(df2))

#df2

#RDD.parameters

for( i in ER.model) {
 ii 			<- df2$ER.model == i
for( j in RDD.model) {
 	jj 			<- ((df2$RDD.model == j) & (df2$ER.model == i)) 
		df2$r.m[jj] <- SGP.RDD.r.m(D.Gy = df2$D.Gy[jj], E.MeV.u, particle.no, material.no, ER.model = i, ER.parameters = ER.parameters, RDD.model = j, RDD.parameters=RDD.parameters[[j]] ) 	
		}
}

#df2

# plots...

p3 <- xyplot( r.m ~ D.Gy | ER.model.name , groups = RDD.model.name, ref = TRUE, data=df2, pch = ".", lty = 1, type = "l", xlab = "Dose [Gy]", ylab = "Distance [m]", auto.key = list(title = "Protons in liquid water, inverse RDD",points = FALSE, lines = TRUE), scales = list(log = 10))
p4 <- xyplot( r.m ~ D.Gy, groups = RDD.model.name, ref = TRUE, data=df2, pch = ".", lty = 1, type = "l", xlab = "Dose [Gy]", ylab = "Distance [m]", auto.key = list(title = "Protons in liquid water, inverse RDD",points = FALSE, lines = TRUE), scales = list(log = 10))

# saving plots to the file

pdf("RDD.pdf")

p1
p2
p3
p4

dev.off()
