rm( list = ls() )
dyn.load("/home/grzanka/workspace/DKFZ/Steffen/SGParticleLibrary/Release/libSGPL.so")
source("SGP.ssc")

# test 1
#model <- "LEM"
#RRD.parameters <- c(60,1,1e-8)
#material.name <- "Water, Liquid"
#parameter <- 1e-8
#expn <- seq (-10, -6, by=0.01)
#r.m <- c(1e-6,1e-7,1e-8,1e-9,1e-10)
#r.m <- 10^expn 
#r.m <- c(1e-7,5e-8)
#res <- SGP.RDD.D.Gy(r.m,model,RRD.parameters,material.name,parameter)
#res
#plot( r.m, res )
#plot( expn, res )

# test 2 - test of RDD
# 0 <- RDD_Test        , parameters: (a0)
# 1 <- RDD_KatzPoint   , parameters: (?)
# 2 <- RDD_Geiss       , parameters: (E_MeV_u, particle_index, material_name, a0)
#RDD.model <- 0
#RRD.parameters <- c(1e-8)
#RDD.model <- 2
#RRD.parameters <- c(60,1,0,1e-8)


# sample distances
#expn <- seq (-10, -6, by=1e0)
#r.m <- 10^expn 

# calculations...
#res <- SGP.RDD.D.Gy.L(r.m,RDD.model,RRD.parameters)

# plotting
#plot( expn, res )

# test 3 - test of SC
#############
E.MeV.u = c(60)
particle.index = c(1)
fluence.cm2 = c(1e8)
material.name <- "Water, Liquid"
rdd.parameter <- c(1e-8)
N2 <- 10
cat("step 1\n")
f.list <- SGP.SC.get.f1(	E.MeV.u, particle.index, fluence.cm2, material.name, rdd.parameter, N2)
write.output <- T
shrink.tails <- F
shrink.tails.under <- 1e-30
adjust.N2 <- F
fluence.factor <- 1.0
cat("step 2\n")
res <- SGP.SC.do.SC(	f.list,	fluence.factor,	write.output, shrink.tails,	shrink.tails.under,	adjust.N2 )
	
