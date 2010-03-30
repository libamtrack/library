################o G################################################################################
# S/R wrapping function interfacing libamtrack library
################################################################################################
# This script replaces (together with the libamtrack.dll) the old S-Plus or mixed S-Plus/C-versions
# v1.0 - v2.x of the AmTrack particle library (or LGC, TIM, SGParticle...)
#
# Some function are not yet implemented in the
# library as C-code but here in S/R
#
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
#
# TODO shall we keep changelog in the script file, or rather rely on commit log messages ?
# changelog does not really tells much to the user
#
# 2006       : Started library under S
# 2009-Apr-27: Started this version of wrapping script, sgre
# 2009-Jun-12: Some minor typo correction before adding to the new repository system, sgre
# 2009-Jun-15: Wrapping function for general RDD methods, debugging printouts, lgrz
# 2009-Jun-16: New wrapping function for r(D) inverse RDDs, sgre
# 2009-Jun-18: New wrapping function for delta electron range models, lgrz
# 2009-Jun-23: Revised wrapping functions for new RDD and SC functions, removed
#              material.name as character cas cause trouble, replaces by material.no (int)
# 2009-Jun-26: Added Edmund transport (CSDA)
# 2010-Feb-08: Minor changes
# 2010-Feb-20: Transport functions removed as they were removed already from the library, lgrz
# 2010-Mar-23: Removed most functions except for those defined in AT_Wrapper_R.h to match the
#              new (but non-SWIG) access from R to the library. As the old version of AmTrack.R
#              is kept in the repository (pre-rev 304) function can be copied into the new 
#              script whenever needed
################################################################################################

################################################################################################
# FUNCTION LIST
#
# AT.gamma.response                          looks OK, to be tested
# AT.LET.MeV.cm2.g                           OK
# AT.max.E.transfer.MeV                      OK
# AT.max.electron.range                      OK
# AT.RDD.D.Gy                                OK
# AT.RDD.D.ext.Gy                            OK
# AT.SC.get.f1.array.size					 OK
# AT.SC.get.f1								 OK
# AT.SC.get.f.array.size                     OK
# 
################################################################################################

debug 		<- F
AT.version 	<- "libamtrack S/R wrapping script - 2010/03/23"

###########
AT.D.RDD.Gy					<-	function(	r.m,
											E.MeV.u,
											particle.no,
											material.no,
											ER.model,
											ER.parameters,
											RDD.model,
											RDD.parameters){
	n					<-	length(r.m)
	D.Gy				<-	numeric(n)
  		
    res					<-	.C(	"AT_D_RDD_Gy_R",	n						=	as.integer(n),
													r.m						=	as.single(r.m),
													E.MeV.u					=	as.single(E.MeV.u),
													particle.no				=	as.integer(particle.no),
													material.no				=	as.integer(material.no),
													RDD.model				=	as.integer(RDD.model),
													RDD.parameters			=	as.single(RDD.parameters),
													ER.model				=	as.integer(ER.model),
													ER.parameters			=	as.single(ER.parameters),
													D.Gy					=	as.single(D.Gy))		
			
	 return(res$D.Gy)						
}

###########
AT.D.RDD.extended.target.Gy	<-	function(	r.m,
											E.MeV.u,
											particle.no,
											material.no,
											ER.model,
											ER.parameters,
											RDD.model,
											RDD.parameters){
	n					<-	length(r.m)
	D.Gy				<-	numeric(n)
  		
    res					<-	.C(	"AT_D_RDD_extended_target_Gy_R",	n						=	as.integer(n),
																	r.m						=	as.single(r.m),
																	E.MeV.u					=	as.single(E.MeV.u),
																	particle.no				=	as.integer(particle.no),
																	material.no				=	as.integer(material.no),
																	RDD.model				=	as.integer(RDD.model),
																	RDD.parameters			=	as.single(RDD.parameters),
																	ER.model				=	as.integer(ER.model),
																	ER.parameters			=	as.single(ER.parameters),
																	D.Gy					=	as.single(D.Gy))		
			
	 return(res$D.Gy)						
}

##################
AT.run.SPIFF	<-	function(	E.MeV.u,
								particle.no,
								fluence.cm2,
								material.no,
								RDD.model,
								RDD.parameters,
								ER.model,
								ER.parameters,
								gamma.model,
								gamma.parameters,
								N2						= 20,
								fluence.factor 			= 1.0,
								write.output 			= F,
								shrink.tails 			= T,
								shrink.tails.under		= 1e-30,
								adjust.N2 				= T,
								lethal.events.mode 		= F){
	
		results			<-	numeric(10)
		N2.tmp			<-	numeric(1)
		N2.tmp			<-	N2
		res				<-	.C(	"AT_run_SPIFF_R",	n					= 	as.integer(length(E.MeV.u)),
													E.MeV.u				=	as.single(E.MeV.u),
													particle.no			=	as.integer(particle.no),
													fluence.cm2			=	as.single(fluence.cm2),
													material.no			=	as.integer(material.no),
													RDD.model			=	as.integer(RDD.model),
													RDD.parameters		=	as.single(RDD.parameters),
													ER.model			=	as.integer(ER.model),
													ER.parameters		=	as.single(ER.parameters),
													gamma.model			=	as.integer(gamma.model),
													gamma.parameters	=	as.single(gamma.parameters),
													N2					=	as.integer(N2.tmp),
													fluence.factor		=	as.single(fluence.factor),
													write.output		=	as.integer(write.output),
													shrink.tails		=	as.integer(shrink.tails),
													shrink.tails.under	=	as.single(shrink.tails.under),
													adjust.N2			=	as.integer(adjust.N2),
													lethal.events.mode 	= 	as.integer(lethal.events.mode),
													results				=	as.single(results))
		results			<-	res$results
	
	return(results)
}
											


#################
AT.gamma.response	<-	function(	d.Gy,
										gamma.model,
										gamma.parameters){
											
	n						<-	length(d.Gy)
	S						<-	numeric(n)
	res						<-	.C(	"AT_gamma_response_R",		n					= as.integer(n),
																d.Gy				= as.single(d.Gy),
																gamma.model			= as.integer(gamma.model),
																gamma.parameters	= as.single(gamma.parameters),
																S					= as.single(S))
	return(res$S)
}

#################
AT.LET.MeV.cm2.g		<-	function(	E.MeV.u,
										particle.no,
										material.no){
	n					<-	length(E.MeV.u)
	LET.MeV.cm2.g		<-	numeric(n)
	res					<-	.C(	"AT_LET_MeV_cm2_g_R",		n					=	as.integer(n),
															E.MeV.u				=	as.single(E.MeV.u),
															particle.no			=	as.integer(particle.no),
															material.no			=	as.integer(material.no),
															LET.MeV.cm2.g		=	as.single(LET.MeV.cm2.g))
	return(res$LET.MeV.cm2.g)						
}

#####################
AT.max.E.transfer.MeV		<-	function(	E.MeV.u ){
	n							<-	length(E.MeV.u)
	max.E.transfer.MeV		<-	numeric(n)
	res							<-	.C(	"AT_max_E_transfer_MeV_R",		n					=	as.integer(n),
																		E.MeV.u				=	as.single(E.MeV.u),
																		max.E.transfer.MeV	=	as.single(max.E.transfer.MeV))
	return(res$max.E.transfer.MeV)						
}

#####################
AT.max.electron.range					<-	function(	E.MeV.u,
												material.no,
												ER.model){
	n					<-	length(E.MeV.u)
	range.m				<-	numeric(n)
	res				<-	.C("AT_max_electron_ranges_m_R", 	n						= as.integer(n),
															E.MeV.u					= as.single(E.MeV.u),
															material.no			= as.integer(material.no),
															ER.model				= as.integer(ER.model),
															range.m					= as.single(range.m))
	return(res$range.m)
}

##########
AT.r.RDD.m					<-	function(	D.Gy,
											E.MeV.u,
											particle.no,
											material.no,
											ER.model,
											ER.parameters,
											RDD.model,
											RDD.parameters){
	n					<-	length(D.Gy)
	r.m					<-	numeric(n)
  
				
		res					<-	.C(	"AT_r_RDD_m_R",	n						=	as.integer(n),
														D.Gy					=	as.single(D.Gy),
														E.MeV.u				=	as.single(E.MeV.u),
														particle.no			=	as.integer(particle.no),
														material.no			=	as.integer(material.no),
														RDD.model				=	as.integer(RDD.model),
														RDD.parameters		=	as.single(RDD.parameters),
														ER.model				=	as.integer(ER.model),
														ER.parameters			=	as.single(ER.parameters),
														r.m						=	as.single(r.m))		
			
	 return(res$r.m)						
}

#######################
AT.SC.get.f1.array.size		<-	function(		E.MeV.u,
												particle.no,
												fluence.cm2,
												material.no,
												RDD.model,
												RDD.parameters,
												ER.model,
												ER.parameters,
												N2){
											
	n				<-	length(E.MeV.u)
	f1.parameters	<-	numeric(9 * n)
	n.bins.f1 <- 0
	
	res				<-	.C("AT_SC_get_f1_array_size_R", 	n					= as.integer(n),
															E.MeV.u				= as.single(E.MeV.u),
															particle.no			= as.integer(particle.no),
															material.no			= as.integer(material.no),
															RDD.model			= as.integer(RDD.model),
															RDD.parameters		= as.single(RDD.parameters),
															ER.model			= as.integer(ER.model),
															ER.parameters		= as.single(ER.parameters),
															N2					= as.integer(N2),
															n.bins.f1			= as.integer(n.bins.f1),
															f1.parameters		= as.single(f1.parameters))
		
	return(list(n.bins.f1 = res$n.bins.f1, f1.parameters = res$f1.parameters))
}

############
AT.SC.get.f1		<-	function(		E.MeV.u,
										particle.no,
										fluence.cm2,
										material.no,
										RDD.model,
										RDD.parameters,
										ER.model,
										ER.parameters,
										N2,
										n.bins.f1,
										f1.parameters){

	n								<-	length(E.MeV.u)
	
	f.parameters					<-	numeric(7)
	
	norm.fluence					<-	numeric(n)
	dose.contribution.Gy			<-	numeric(n)
	
	f1.d.Gy							<-	numeric(n.bins.f1)
	f1.dd.Gy						<-	numeric(n.bins.f1)
	f1								<-	numeric(n.bins.f1)	

	res				<-	.C("AT_SC_get_f1_R", 					n						= as.integer(n),
																E.MeV.u					= as.single(E.MeV.u),
																particle.no				= as.integer(particle.no),
																fluence.cm2				= as.single(fluence.cm2),
																material.no				= as.integer(material.no),
																RDD.model				= as.integer(RDD.model),
																RDD.parameters			= as.single(RDD.parameters),
																ER.model				= as.integer(ER.model),
																ER.parameters			= as.single(ER.parameters),
																N2						= as.integer(N2),
																n.bins.f1				= as.integer(n.bins.f1),
																f1.parameters			= as.single(f1.parameters),
																norm.fluence			= as.single(norm.fluence),
																dose.contribution.Gy	= as.single(dose.contribution.Gy),
																f.parameters			= as.single(f.parameters),
																f1.d.Gy					= as.single(f1.d.Gy),
																f1.dd.Gy				= as.single(f1.dd.Gy),
																f1						= as.single(f1))
	
	results	<-	list(	fluence.cm2					= 	res$fluence.cm2,
						norm.fluence				= 	res$norm.fluence,
						dose.contribution.Gy		= 	res$dose.contribution.Gy,
						f.parameters				=	res$f.parameters,											
						f1							=	data.frame(	f1.d.Gy						=	res$f1.d.Gy,
																	f1.dd.Gy					=	res$f1.dd.Gy,
																	f1							=	res$f1))
	return(results)
}


######################
AT.SC.get.f.array.size		<-	function(		u,
												fluence.factor,
												N2,
												n.bins.f1,
												f1.d.Gy,
												f1.dd.Gy,
												f1){
											
	n.bins.f		<-	numeric(1)
	u.start			<-	numeric(1)
	n.convolutions	<-	numeric(1)
	
	res				<-	.C("AT_SC_get_f_array_size_R", 		u					= as.single(u),
															fluence.factor		= as.single(fluence.factor),
															N2					= as.integer(N2),
															n.bins.f1			= as.integer(n.bins.f1),
															f1.d.Gy				= as.single(f1.d.Gy),
															f1.dd.Gy			= as.single(f1.dd.Gy),
															f1					= as.single(f1),
															n.bins.f			= as.integer(n.bins.f),
															u.start				= as.single(u.start),
															n.convolutions		= as.integer(n.convolutions))
		
	return(list(n.bins.f = res$n.bins.f, u.start = res$u.start, n.convolutions = res$n.convolutions))
}

######################
AT.SC.get.f.start		<-	function(		u,
											N2,
											n.bins.f1,
											f1.d.Gy,
											f1.dd.Gy,
											f1,
											n.bins.f){
											
	f.d.Gy							<-	numeric(n.bins.f)
	f.dd.Gy							<-	numeric(n.bins.f)
	f								<-	numeric(n.bins.f)	
	
	res						<-	.C("AT_SC_get_f_start_R", 		u						= as.single(u),
																n.bins.f1				= as.integer(n.bins.f1),
																N2						= as.integer(N2),
																f1.d.Gy					= as.single(f1.d.Gy),
																f1.dd.Gy				= as.single(f1.dd.Gy),
																f1						= as.single(f1),
																n.bins.f				= as.integer(n.bins.f),
																f.d.Gy					= as.single(f.d.Gy),
																f.dd.Gy					= as.single(f.dd.Gy),
																f						= as.single(f))

	return	(			data.frame(	f.start.d.Gy					=	res$f.d.Gy,
								f.start.dd.Gy					=	res$f.dd.Gy,
								f.start						=	res$f))
}

############################
AT.SC.SuccessiveConvolutions		<-	function(	u,
													n.bins.f,
													N2,
													n.bins.f.used,
													f.d.Gy,
													f.dd.Gy,
													f,
													write.output,
													shrink.tails,
													shrink.tails.under,
													adjust.N2){

	f0					<-	numeric(1)
	fdd					<-	numeric(n.bins.f)
	dfdd				<-	numeric(n.bins.f)
	d					<-	numeric(1)
	
	res						<-	.C("AT_SuccessiveConvolutions_R", 		u						= as.single(u),
																		n.bins.f				= as.integer(n.bins.f),
																		N2						= as.integer(N2),
																		n.bins.f.used			= as.integer(n.bins.f.used),
																		f.d.Gy					= as.single(f.d.Gy),
																		f.dd.Gy					= as.single(f.dd.Gy),
																		f						= as.single(f),
																		f0						= as.single(f0),
																		fdd						= as.single(fdd),
																		dfdd					= as.single(dfdd),
																		d						= as.single(d),
																		write.output			= as.integer(write.output),
																		shrink.tails			= as.integer(shrink.tails),
																		shrink.tails.under		= as.single(shrink.tails.under),
																		adjust.N2				= as.integer(adjust.N2))

	return(				list(	N2				=	res$N2,
								n.bins.f.used	=	res$n.bins.f.used,
								f0				=	res$f0,
								d				=	res$d,
								f				=	data.frame(	f.d.Gy						=	res$f.d.Gy,
																f.dd.Gy						=	res$f.dd.Gy,
																f							=	res$f,
																fdd							=	res$fdd,
																dfdd						=	res$dfdd)))																	
}
	
######################################################################################################################################

print(AT.version)
print(" ")
print("List of functions:")
print(sort(ls()[grep("AT.", ls())]))
print("end.")