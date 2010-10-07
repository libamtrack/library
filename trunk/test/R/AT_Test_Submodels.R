################################################################################################
# S/R test script for implemented submodels: RDDs, ERs, gamma response
################################################################################################
# Copyright 2006, 2010 The libamtrack team
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
################################################################################################

# clear workspace
rm( list = ls() )

# load libAmTrack library
try(dyn.load("../../lib/libamtrack.dll"))
try(dyn.load("../../lib/libamtrack.so"))
try(dyn.load("../../lib/libamtrack.dylib"))

# load wrapping scripts
source("../../wrapper/R/AmTrack.R")

# necessary library for plotting
library("lattice")

debug				<-	F

################################################################################################
# RDDs
################################################################################################
E.MeV.u			<-	c(1, 10, 100)
particle.no		<-	c(1001, 2004, 6012)					# p, He-4, C-12

material.no		<-	c(1, 2)

RDD.model <- c(1,2,3,4,5,6,7)
RDD.model.names <- c("Simple step test function",  "Katz' point target", "Geiss'", "Site", "Cucinotta", "KatzExtTarget", "CucinottaExtTarget")
RDD.parameters <- list(c(1),c(1e-10,1e-10), c(1e-10),c(1e-8,1e-10),c(5e-11,1e-10),c(1e-10,1e-8,1e-10),c(5e-11,1e-8,1e-10))

ER.model			<-	2 								# (Waligorski)

df.RDD				<-	expand.grid(	r.m				=	10^seq(-12, 0, by = 0.02),
										E.MeV.u		=	E.MeV.u,
										particle.no	=	particle.no,
										material.no	=	material.no,
										RDD.model	= 	RDD.model)

df.RDD$RDD.model.name	<-	as.character(RDD.model.names[df.RDD$RDD.model])

df.RDD$D.Gy				<-	numeric(nrow(df.RDD))

# Conditioning variable
df.RDD$EPMM				<-	paste(df.RDD$E.MeV.u, df.RDD$particle.no, df.RDD$material.no, df.RDD$RDD.model, sep = "-")

for (cur.EPMM in unique(df.RDD$EPMM)){
	#cur.EPMM<-unique(df$EPMM)[1]
	ii					<-	cur.EPMM == df.RDD$EPMM
	df.RDD$D.Gy[ii]	<-	AT.D.RDD.Gy(		r.m					=	df.RDD$r.m[ii],
												E.MeV.u			=	df.RDD$E.MeV.u[ii],
												particle.no		=	df.RDD$particle.no[ii],
												material.no		=	df.RDD$material.no[ii],
												ER.model			=	ER.model,
												RDD.model			=	unique(df.RDD$RDD.model[ii]),
												RDD.parameters	=	RDD.parameters[[unique(df.RDD$RDD.model[ii])]])			
}

p1 <- xyplot(		data = df.RDD,
			log10(D.Gy) ~ log10(r.m)|paste(E.MeV.u, "MeV/u")*paste("particle no.", sprintf("%02d",particle.no))*paste("material no.", material.no),
			groups		=	RDD.model.name,
			type		=	'l',
			as.table	=	T)

p2 <- xyplot(		data = df.RDD,
			log10(D.Gy) ~ log10(r.m)|RDD.model.name*paste("particle no.", sprintf("%02d",particle.no)),
			groups		=	paste(E.MeV.u, "MeV/u"),
			type		=	'l',
			as.table	=	T)
							
################################################################################################
# ERs
################################################################################################

E.MeV.u 					<- 10^seq (0, 3, by = 0.1)

ER.models 					<- c(1, 2, 3, 4, 5)
ER.model.names 			<- c("Test", "Butts & Katz", "Waligorski", "Geiss", "Scholz")

n 							<- length(E.MeV.u) 
material.no	 			<- 1 				# Water

df.ER						<- expand.grid( 	E.MeV.u 			= E.MeV.u , 
												ER.model 			= ER.models,
												material.no		= material.no)

df.ER$ER.model.name		<- as.character(ER.model.names[df.ER$ER.model])
df.ER$range.m				<- numeric(nrow(df.ER))

# Conditioning variable
df.ER$PMM					<-	paste(df.ER$material.no, df.ER$ER.model, sep = "-")

for(cur.PMM in unique(df.ER$PMM) ){
 	#cur.PMM<-unique(df.ER$PMM)[1]
	ii					<-	cur.PMM == df.ER$PMM
	df.ER$range.m[ii] 		<- 	AT.max.electron.range(	E.MeV.u 		= df.ER$E.MeV.u[ii], 
																material.no	= unique(df.ER$material.no[ii]), 
																ER.model		= unique(df.ER$ER.model[ii]))
}

p3 <- xyplot(		data = df.ER,
			log10(range.m) ~ log10(E.MeV.u)|paste("material no.", material.no),
			groups		=	ER.model.name,
			type		=	'l',
			as.table	=	T)

################################################################################################
# gamma responses
################################################################################################

D.Gy	 					<- 10^seq (-2, 3, by = 0.1)

GR.models 					<- c(1, 2, 3, 4, 5)
GR.model.names 			<- c(		"Test", 
										"General hit/target", 
										"Radioluminescence", 
										"Exp.-saturation", 
										"Linear-quadratic")
R									<-	1
Smax								<-	1
k1									<-	Smax * (R / 100)
k2									<-	Smax * (1 - R / 100)
Jens.gamma.parameters.peak.A	<-	c(	k1 = k1, D01 = 0.36, c1 = 1, m1 = 1,
											k2 = k2, D02 = 3.06, c2 = 1, m2 = 2,
											0)

GR.parameters				<- list(	c(1, 0.1), 
										Jens.gamma.parameters.peak.A,
										c(1, 10, 10),
										c(1, 10),
										c(1, 1))
 

df.GR						<- expand.grid( 	D.Gy	 			= D.Gy, 
												GR.model 			= GR.models)

df.GR$GR.model.name		<- as.character(GR.model.names[df.GR$GR.model])
df.GR$S					<- numeric(nrow(df.GR))

for(cur.M in unique(df.GR$GR.model)){
 	#cur.M <- unique(df.GR$GR.model)[1]
	ii					<-	cur.M == df.GR$GR.model
	df.GR$S[ii] 		<- 	AT.gamma.response(	d.Gy	 			= df.GR$D.Gy[ii], 
													gamma.model		= unique(df.GR$GR.model[ii]),
													gamma.parameters	= GR.parameters[[unique(df.GR$GR.model[ii])]],
													lethal.events.mode = F)
}

p4 <- xyplot(		data = df.GR,
			log10(S) ~ log10(D.Gy),
			groups		=	GR.model.name,
			type		=	'l',
			ylim		=  c(-3, 3),
			xlim		=  c(-2, 3),
			as.table	=	T)
			
			
pdf("Submodels.pdf")

p1
p2
p3
p4

dev.off()
