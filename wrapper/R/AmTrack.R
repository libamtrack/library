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
# 2010-Jun-01: Added IGK methode, and variable names for all methods' result arrays
################################################################################################

debug 		<- F
AT.version 	<- "libamtrack S/R wrapping script - 2010/07/21"

##########################
AT.convert.beam.parameters	<-	function(	N 				= 0,
											FWHM.mm			= 0,
											fluence.cm2		= 0,
											sigma.cm		= 0){
	n					<-	max(length(N), length(FWHM.mm), length(fluence.cm2), length(sigma.cm))
	if(N == 0){
		N					<-	numeric(n)
	}
	if(FWHM.mm == 0){
		FWHM.mm				<-	numeric(n)
	}
	if(fluence.cm2 == 0){
		fluence.cm2			<-	numeric(n)
	}
	if(sigma.cm == 0){
		sigma.cm			<-	numeric(n)
	}
	
	res					<-	.C(	"AT_convert_beam_parameters_R",						n					=	as.integer(n),
																					fluence.cm2			=	as.single(fluence.cm2),
																					sigma.cm			=	as.single(sigma.cm),
																					N					=	as.single(N),
																					FWHM.mm				=	as.single(FWHM.mm))
	return(res)						
}

#################
AT.CSDA.range.g.cm2		<-	function(	E.MeV.u,
										particle.no,
										material.no){
	n					<-	length(E.MeV.u)
	CSDA.range.g.cm2	<-	numeric(n)
	res					<-	.C(	"AT_CSDA_range_g_cm2_R",		n					=	as.integer(n),
															    E.MeV.u				=	as.single(E.MeV.u),
															    particle.no			=	as.integer(particle.no),
															    material.no			=	as.integer(material.no),
															    CSDA.range.g.cm2	=	as.single(CSDA.range.g.cm2))
	return(res$CSDA.range.g.cm2)						
}

#############
AT.D.Gy					<-	function(	E.MeV.u,
										particle.no,
										fluence.cm2,
										material.no){
	n					<-	length(E.MeV.u)
	D.Gy			<-	numeric(n)
  		
    res					<-	.C(	"AT_D_Gy_R",	n						=	as.integer(n),
												E.MeV.u					=	as.single(E.MeV.u),
												particle.no				=	as.integer(particle.no),
												fluence.cm2				=	as.single(fluence.cm2),
												material.no				=	as.integer(material.no),
												D.Gy					=	as.single(D.Gy))		
			
	 return(res$D.Gy)						
}

###########
AT.D.RDD.Gy					<-	function(	r.m,
											E.MeV.u,
											particle.no,
											material.no,
											ER.model,
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
													D.Gy					=	as.single(D.Gy))		
			
	 return(res$D.Gy)						
}

###########
AT.RDD.f1.parameters.mixed.field					<-	function(	E.MeV.u,
											particle.no,
											material.no,
											ER.model,
											RDD.model,
											RDD.parameters){
	AT.SC.F1.PARAMETERS.SINGLE.LENGTH	<-	8
	n						<-	length(E.MeV.u)
	f1.parameters				<-	numeric(n * AT.SC.F1.PARAMETERS.SINGLE.LENGTH)
  		
    res					<-	.C(	"AT_RDD_f1_parameters_mixed_field_R",	n						=	as.integer(n),
													E.MeV.u					=	as.single(E.MeV.u),
													particle.no				=	as.integer(particle.no),
													material.no				=	as.integer(material.no),
													RDD.model				=	as.integer(RDD.model),
													RDD.parameters			=	as.single(RDD.parameters),
													ER.model				=	as.integer(ER.model),
													f1.parameters			=	as.single(f1.parameters))		
			
	df					<-	data.frame(	LET.MeV.cm2.g			= numeric(0),
									r.min.m				= numeric(0),
									r.max.m				= numeric(0),
									d.min.Gy				= numeric(0),
									d.max.Gy				= numeric(0),
									normalization.constant		= numeric(0),
									single.impact.fluence.cm2	= numeric(0),
									single.impact.dose.Gy		= numeric(0))

	for(i in 1:n){
		# i <- 2
		df[i,]	<-	res$f1.parameters[(((i - 1) * AT.SC.F1.PARAMETERS.SINGLE.LENGTH)+1):(i * AT.SC.F1.PARAMETERS.SINGLE.LENGTH)]
	}
	
	return(df)						
}

#############
AT.fluence.cm2			<-	function(	E.MeV.u,
										particle.no,
										D.Gy,
										material.no){
	n					<-	length(E.MeV.u)
	fluence.cm2			<-	numeric(n)
  		
    res					<-	.C(	"AT_fluence_cm2_R",	n						=	as.integer(n),
												E.MeV.u					=	as.single(E.MeV.u),
												particle.no				=	as.integer(particle.no),
												D.Gy					=	as.single(D.Gy),
												material.no				=	as.integer(material.no),
												fluence.cm2				=	as.single(fluence.cm2))		
			
	 return(res$fluence.cm2)						
}

##################
AT.run.GSM.method	<-	function(	E.MeV.u,
									particle.no,
									fluence.cm2.or.dose.Gy,
									material.no,
									RDD.model,
									RDD.parameters,
									ER.model,
									gamma.model,
									gamma.parameters,
									N.runs,
									write.output,
									nX,
									voxel.size.m,
									lethal.events.mode){
		results			<-	numeric(10)
		res				<-	.C(	"AT_run_GSM_method_R",	n					= 	as.integer(length(E.MeV.u)),
														E.MeV.u				=	as.single(E.MeV.u),
														particle.no			=	as.integer(particle.no),
														fluence.cm2.or.dose.Gy			=	as.single(fluence.cm2.or.dose.Gy),
														material.no			=	as.integer(material.no),
														RDD.model			=	as.integer(RDD.model),
														RDD.parameters		=	as.single(RDD.parameters),
														ER.model			=	as.integer(ER.model),
														gamma.model			=	as.integer(gamma.model),
														gamma.parameters	=	as.single(gamma.parameters),
														N.runs				=	as.integer(N.runs),
														write.output		=	as.integer(write.output),
														nX					=	as.integer(nX),
														voxel.size.m		=	as.single(voxel.size.m),
														lethal.events.mode 	= 	as.integer(lethal.events.mode),
														results				=	as.single(results))
		results			<-	res$results
		names(results)		<-	c(	"efficiency", 
								"D.check.Gy", 
								"response.HCP", 
								"response.gamma", 
								"n.particles",	
								"err.efficiency", 
								"err.D.check.Gy", 
								"err.response.HCP", 
								"err.response.gamma", 
								"err.n.particles")
	return(results)
}
											
##################
AT.run.SPIFF.method	<-	function(	E.MeV.u,
									particle.no,
									fluence.cm2.or.dose.Gy,
									material.no,
									RDD.model,
									RDD.parameters,
									ER.model,
									gamma.model,
									gamma.parameters,
									N2,
									fluence.factor,
									write.output,
									shrink.tails,
									shrink.tails.under,
									adjust.N2,
									lethal.events.mode){
	
		results			<-	numeric(10)
		N2.tmp			<-	numeric(1)
		N2.tmp			<-	N2
		res				<-	.C(	"AT_run_SPIFF_method_R",	n					= 	as.integer(length(E.MeV.u)),
															E.MeV.u				=	as.single(E.MeV.u),
															particle.no			=	as.integer(particle.no),
															fluence.cm2.or.dose.Gy			=	as.single(fluence.cm2.or.dose.Gy),
															material.no			=	as.integer(material.no),
															RDD.model			=	as.integer(RDD.model),
															RDD.parameters		=	as.single(RDD.parameters),
															ER.model			=	as.integer(ER.model),
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
		names(results)		<-	c(	"efficiency", 
								"D.check.Gy", 
								"response.HCP", 
								"response.gamma", 
								"---",	
								"u", 
								"u.start", 
								"n.convolutions", 
								"lower.Jensen.bound", 
								"upper.Jensen.bound")
	
	return(results)
}

											
##################
AT.run.IGK.method	<-	function(	E.MeV.u,
									particle.no,
									fluence.cm2.or.dose.Gy,
									material.no,
									RDD.model,
									RDD.parameters,
									ER.model,
									gamma.model,
									gamma.parameters,
									saturation.cross.section.factor,
									write.output){
	
		results			<-	numeric(10)
		res				<-	.C(	"AT_run_IGK_method_R",		n					= 	as.integer(length(E.MeV.u)),
															E.MeV.u				=	as.single(E.MeV.u),
															particle.no			=	as.integer(particle.no),
															fluence.cm2.or.dose.Gy			=	as.single(fluence.cm2.or.dose.Gy),
															material.no			=	as.integer(material.no),
															RDD.model			=	as.integer(RDD.model),
															RDD.parameters		=	as.single(RDD.parameters),
															ER.model			=	as.integer(ER.model),
															gamma.model			=	as.integer(gamma.model),
															gamma.parameters	=	as.single(gamma.parameters),
															saturation.cross.section.factor 	=	as.single(saturation.cross.section.factor),
															write.output		=	as.integer(write.output),
															results				=	as.single(results))
		results			<-	res$results
		names(results)		<-	c(	"efficiency", 
								"---", 
								"response.HCP", 
								"response.gamma", 
								"---",	
								"ion.cross.section.cm2", 
								"gamma.dose.Gy", 
								"pi.Ion", 
								"pi.Gamma", 
								"---")
	
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

################################
AT.Katz.inactivation.probability		<-	function(	r.m,
													E.MeV.u,
													particle.no,
													material.no,
													rdd.model,
													rdd.parameters,
													er.model,
													gamma.parameters){

	n					<-	length(r.m)
	inactivation.probability					<-	numeric(n)

	res						<-	.C("AT_KatzModel_inactivation_probability_R", 		n						= as.integer(n),
																		r.m					               			= as.single(r.m),
																		E.MeV.u								         		= as.single(E.MeV.u),
																		particle.no			 				       = as.integer(particle.no),
																		material.no								      	= as.integer(material.no),
																		rdd.model				        					= as.integer(rdd.model),
																		rdd.parameters						      = as.single(rdd.parameters),
																		er.model	            					= as.integer(er.model),
																		gamma.parameters					    	= as.single(gamma.parameters),
																		inactivation.probability	 = as.single(inactivation.probability))

	 return(res$inactivation.probability)
}
	
#####################################
AT.Katz.inactivation.cross.section.m2		<-	function(		E.MeV.u,
													particle.no,
													material.no,
													rdd.model,
													rdd.parameters,
													er.model,
													gamma.parameters){

	n					<-	length(E.MeV.u)
	inactivation.cross.section.m2					<-	numeric(n)

	res						<-	.C("AT_KatzModel_inactivation_cross_section_m2_R", 		n						= as.integer(n),
																		E.MeV.u								         		= as.single(E.MeV.u),
																		particle.no			 				       = as.integer(particle.no),
																		material.no								      	= as.integer(material.no),
																		rdd.model				        					= as.integer(rdd.model),
																		rdd.parameters						      = as.single(rdd.parameters),
																		er.model	            					= as.integer(er.model),
																		gamma.parameters					    	= as.single(gamma.parameters),
																		inactivation.cross.section.m2	 = as.single(inactivation.cross.section.m2))

	 return(res$inactivation.cross.section.m2)
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

#################
AT.LET.keV.um		<-	function(	E.MeV.u,
										particle.no,
										material.no){
	n					<-	length(E.MeV.u)
	LET.keV.um			<-	numeric(n)
	res					<-	.C(	"AT_LET_keV_um_R",			n					=	as.integer(n),
															E.MeV.u				=	as.single(E.MeV.u),
															particle.no			=	as.integer(particle.no),
															material.no			=	as.integer(material.no),
															LET.keV.um			=	as.single(LET.keV.um))
	return(res$LET.keV.um)						
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

#################################
AT.particle.name.from.particle.no		<-	function(		particle.no){

	n					<-	length(particle.no)
	particle.name		<-	character(n)
	
	for (i in 1:n){
		cur.particle.name	<-	character(1)
		res					<-	.C("AT_particle_name_from_particle_no_R", 		particle.no				= as.integer(particle.no[i]),
																	particle.name			= as.character(cur.particle.name))
		particle.name[i]	<-	res$particle.name
	}		
	return(particle.name)
}
	
#################################	
AT.particle.no.from.particle.name		<-	function(		particle.name ){

	n					<-	length(particle.name)
	particle.no		<-	numeric(n)
	
	for (i in 1:n){
		cur.particle.no	<-	numeric(1)
		res					<-	.C("AT_particle_no_from_particle_name_R", 		particle.name				= as.character(particle.name[i]),  
																	 particle.no			= as.integer(cur.particle.no))
		particle.no[i]	<-	res$particle.no
	}		
	return(particle.no)
}

#################################
AT.material.name.from.material.no		<-	function(		material.no){

	n					<-	length(material.no)
	material.name		<-	character(n)
	
	for (i in 1:n){
		cur.material.name	<-	character(1)
		res					<-	.C("AT_material_name_from_number_R", 		material.no				= as.integer(material.no[i]),
																	material.name			= as.character(cur.material.name))
		material.name[i]	<-	res$material.name
	}		
	return(material.name)
}
	
#################################	
AT.material.no.from.material.name		<-	function(		material.name ){

	n					<-	length(material.name)
	material.no		<-	numeric(n)
	
	for (i in 1:n){
		cur.material.no	<-	numeric(1)
		res					<-	.C("AT_material_number_from_name_R", 		material.name				= as.character(material.name[i]),  
																	 material.no			= as.integer(cur.material.no))
		material.no[i]	<-	res$material.no
	}		
	return(material.no)
}

#################################
AT.particle.no.from.Z.and.A		<-	function(	Z, A){

	n					<-	length(Z)
	particle.no			<-	numeric(n)
	
	res					<-	.C("AT_particle_no_from_Z_and_A_R", 		n  						= as.integer(n),
																		Z  						= as.integer(Z),
																		A  						= as.integer(A),
																		particle.no				= as.integer(particle.no))

	return(res$particle.no)
}
	
##########
AT.r.RDD.m					<-	function(	D.Gy,
											E.MeV.u,
											particle.no,
											material.no,
											ER.model,
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
														r.m						=	as.single(r.m))		
			
	 return(res$r.m)						
}

#######################
AT.SC.get.f1.array.size		<-	function(		E.MeV.u,
												particle.no,
												fluence.cm2.or.dose.Gy,  # TODO this parameter is not needed !
												material.no,
												RDD.model,
												RDD.parameters,
												ER.model,
												N2){
											
	n				<-	length(E.MeV.u)
	f1.parameters	<-	numeric(8 * n)
	n.bins.f1 <- 0
	
	res				<-	.C("AT_SC_get_f1_array_size_R", 	n					= as.integer(n),
															E.MeV.u				= as.single(E.MeV.u),
															particle.no			= as.integer(particle.no),
															material.no			= as.integer(material.no),
															RDD.model			= as.integer(RDD.model),
															RDD.parameters		= as.single(RDD.parameters),
															ER.model			= as.integer(ER.model),
															N2					= as.integer(N2),
															n.bins.f1			= as.integer(n.bins.f1),
															f1.parameters		= as.single(f1.parameters))
		
	return(list(n.bins.f1 = res$n.bins.f1, f1.parameters = res$f1.parameters))
}

############
AT.SC.get.f1		<-	function(		E.MeV.u,
										particle.no,
										fluence.cm2.or.dose.Gy,
										material.no,
										RDD.model,
										RDD.parameters,
										ER.model,
										N2,
										n.bins.f1,
										f1.parameters){

	n								<-	length(E.MeV.u)
		
	f1.d.Gy							<-	numeric(n.bins.f1)
	f1.dd.Gy						<-	numeric(n.bins.f1)
	f1								<-	numeric(n.bins.f1)	

	res				<-	.C("AT_SC_get_f1_R", 					n						= as.integer(n),
																E.MeV.u					= as.single(E.MeV.u),
																particle.no				= as.integer(particle.no),
																fluence.cm2.or.dose.Gy	= as.single(fluence.cm2.or.dose.Gy),
																material.no				= as.integer(material.no),
																RDD.model				= as.integer(RDD.model),
																RDD.parameters			= as.single(RDD.parameters),
																ER.model				= as.integer(ER.model),
																N2						= as.integer(N2),
																n.bins.f1				= as.integer(n.bins.f1),
																f1.parameters			= as.single(f1.parameters),
																f1.d.Gy					= as.single(f1.d.Gy),
																f1.dd.Gy				= as.single(f1.dd.Gy),
																f1						= as.single(f1))
	
	results	<-	list(	fluence.cm2					= 	res$fluence.cm2,
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
AT.SC.get.f.start		<-	function(		N2,
											n.bins.f1,
											f1.d.Gy,
											f1.dd.Gy,
											f1,
											n.bins.f){
											
	f.d.Gy							<-	numeric(n.bins.f)
	f.dd.Gy							<-	numeric(n.bins.f)
	f								<-	numeric(n.bins.f)	
	
	res						<-	.C("AT_SC_get_f_start_R", n.bins.f1				= as.integer(n.bins.f1),
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

########################
AT.SC.get.gamma.response		<-	function(	d.Gy,
												dd.Gy,
												f,
												f0,
												gamma.model,
												gamma.parameters,
												lethal.events.mode){

	number.of.bins		<-	length(d.Gy)
	S					<-	numeric(number.of.bins)
	S.HCP				<-	numeric(1)
	S.gamma				<-	numeric(1)
	efficiency			<-	numeric(1)
	
	res						<-	.C("AT_SC_get_gamma_response_R", 		number.of.bins			= as.integer(number.of.bins),
																		d.Gy  					= as.single(d.Gy),
																		dd.Gy					= as.single(dd.Gy),
																		f						= as.single(f),
																		f0						= as.single(f0),
																		gamma.model				= as.integer(gamma.model),
																		gamma.parameters		= as.single(gamma.parameters),
																		lethal.events.mode		= as.integer(lethal.events.mode),
																		S						= as.single(S),
																		S.HCP					= as.single(S.HCP),
																		S.gamma					= as.single(S.gamma),
																		efficiency				= as.single(efficiency))
																		
																		
	return(				list(	S.HCP			=	res$S.HCP,
								S.gamma			=	res$S.gamma,
								efficiency		=	res$efficiency,
								S				=	data.frame(	d.Gy						=	d.Gy,
																dd.Gy						=	dd.Gy,
																f							=	f,
																S							=	res$S)))																	
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

#############
AT.fluence.weighted.stopping.power.ratio.R					<-	function(	E.MeV.u,
											particle.no,
											fluence.cm2,
											material.no,
											reference.material.no){
	n										<-	length(E.MeV.u)
	fluence.weighted.stopping.power.ratio	<-	numeric(1)
  		
    res					<-	.C(	"AT_fluence_weighted_stopping_power_ratio_R",	n						=	as.integer(n),
																				E.MeV.u					=	as.single(E.MeV.u),
																				particle.no				=	as.integer(particle.no),
																				fluence.cm2				=	as.single(fluence.cm2),
																				material.no				=	as.integer(material.no),
																				reference.material.no	=	as.integer(reference.material.no),
																				result					=	as.single(fluence.weighted.stopping.power.ratio))		
			
	 return(res$result)						
}

#############
AT.total.D.Gy					<-	function(	E.MeV.u,
											particle.no,
											fluence.cm2,
											material.no){
	n					<-	length(E.MeV.u)
	total.D.Gy			<-	numeric(1)
  		
    res					<-	.C(	"AT_total_D_Gy_R",	n						=	as.integer(n),
													E.MeV.u					=	as.single(E.MeV.u),
													particle.no				=	as.integer(particle.no),
													fluence.cm2				=	as.single(fluence.cm2),
													material.no				=	as.integer(material.no),
													total.D.Gy				=	as.single(total.D.Gy))		
			
	 return(res$total.D.Gy)						
}

#####################
AT.get.materials.data		<-	function( material.no){

	number.of.materials			<-	length(material.no)
	density.g.cm3				<-	numeric(number.of.materials)
    electron.density.m3			<-	numeric(number.of.materials)
    I.eV						<-	numeric(number.of.materials)
    alpha.g.cm2.MeV				<-	numeric(number.of.materials)
    p.MeV						<-	numeric(number.of.materials)
    m.g.cm2						<-	numeric(number.of.materials)
    average.A					<-	numeric(number.of.materials)
    average.Z					<-	numeric(number.of.materials)
	
	res							<-	.C("AT_get_materials_data_R", 	number.of.materials		= as.integer(number.of.materials),
															material.no				= as.integer(material.no),		
															density.g.cm3			=	as.single(density.g.cm3),
																	electron.density.m3		=	as.single(electron.density.m3),
																	I.eV					=	as.single(I.eV),
																	alpha.g.cm2.MeV			=	as.single(alpha.g.cm2.MeV),
																	p.MeV					=	as.single(p.MeV),
																	m.g.cm2					=	as.single(m.g.cm2),
																	average.A				=	as.single(average.A),
																	average.Z				=	as.single(average.Z))

	return(	data.frame(	material.no				=	material.no,
						density.g.cm3			=	res$density.g.cm3,
						electron.density.m3		=	res$electron.density.m3,
						I.eV					=	res$I.eV,
						alpha.g.cm2.MeV			=	res$alpha.g.cm2.MeV,
						p.MeV					=	res$p.MeV,
						m.g.cm2					=	res$m.g.cm2,
						average.A				=	res$average.A,
						average.Z				=	res$average.Z))
}

#####################
AT.A.from.particle.no		<-	function(	particle.no){
	
	n					<-	length(particle.no)
	A					<-	integer(n)
  		
    res					<-	.C(	"AT_A_from_particle_no_R",	n						=	as.integer(n),
													particle.no				=	as.integer(particle.no),
													A						=	as.integer(A))		
			
	 return(res$A)						
}

#####################
AT.Z.from.particle.no		<-	function(	particle.no){
	
	n					<-	length(particle.no)
	Z					<-	integer(n)
  		
    res					<-	.C(	"AT_Z_from_particle_no_R",	n						=	as.integer(n),
													particle.no				=	as.integer(particle.no),
													Z						=	as.integer(Z))		
			
	 return(res$Z)						
}

###########################
AT.fluence.weighted.E.MeV.u		<-	function(	E.MeV.u,
												fluence.cm2)
{
	n							<-	length(E.MeV.u)
	fluence.weighted.E.MeV.u	<-	numeric(1)
    res					<-	.C(	"AT_fluence_weighted_E_MeV_u_R",	n							=	as.integer(n),
																	E.MeV.u						=	as.single(E.MeV.u),
																	fluence.cm2					=	as.single(fluence.cm2),
																	fluence.weighted.E.MeV.u	=	as.single(fluence.weighted.E.MeV.u))		
	
	return(res$fluence.weighted.E.MeV.u)
}
								
########################
AT.dose.weighted.E.MeV.u		<-	function(	E.MeV.u,
												particle.no,
												fluence.cm2,
												material.no)
{
	n							<-	length(E.MeV.u)
	dose.weighted.E.MeV.u		<-	numeric(1)
    res							<-	.C(	"AT_dose_weighted_E_MeV_u_R",	n							=	as.integer(n),
																		E.MeV.u						=	as.single(E.MeV.u),
																		particle.no					=	as.integer(particle.no),
																		fluence.cm2					=	as.single(fluence.cm2),
																		material.no					=	as.integer(material.no),
																		dose.weighted.E.MeV.u		=	as.single(dose.weighted.E.MeV.u))		
	return(res$dose.weighted.E.MeV.u)
}
			
#################################
AT.fluence.weighted.LET.MeV.cm2.g	<-	function(	E.MeV.u,
													particle.no,
													fluence.cm2,
													material.no)
{
	n								<-	length(E.MeV.u)
	fluence.weighted.LET.MeV.cm2.g	<-	numeric(1)
    res								<-	.C(	"AT_fluence_weighted_LET_MeV_cm2_g_R",	n								=	as.integer(n),
																					E.MeV.u							=	as.single(E.MeV.u),
																					particle.no						=	as.integer(particle.no),
																					fluence.cm2						=	as.single(fluence.cm2),
																					material.no						=	as.integer(material.no),
																					fluence.weighted.LET.MeV.cm2.g	=	as.single(fluence.weighted.LET.MeV.cm2.g))		
	return(res$fluence.weighted.LET.MeV.cm2.g)
}
			
##############################
AT.dose.weighted.LET.MeV.cm2.g		<-	function(	E.MeV.u,
													particle.no,
													fluence.cm2,
													material.no)
{
	n							<-	length(E.MeV.u)
	dose.weighted.LET.MeV.cm2.g	<-	numeric(1)
    res							<-	.C(	"AT_dose_weighted_LET_MeV_cm2_g_R",	n							=	as.integer(n),
																			E.MeV.u						=	as.single(E.MeV.u),
																			particle.no					=	as.integer(particle.no),
																			fluence.cm2					=	as.single(fluence.cm2),
																			material.no					=	as.integer(material.no),
																			dose.weighted.LET.MeV.cm2.g	=	as.single(dose.weighted.LET.MeV.cm2.g))		
	
	return(res$dose.weighted.LET.MeV.cm2.g)
}

##############################
AT.total.u      <-  function(   E.MeV.u,
                                                    particle.no,
                                                    fluence.cm2,
                                                    material.no,
                                                    er.model)
{
    n                           <-  length(E.MeV.u)
    u                           <-  numeric(1)
    res                         <-  .C( "AT_total_u_R", n                           =   as.integer(n),
                                                                            E.MeV.u                     =   as.single(E.MeV.u),
                                                                            particle.no                 =   as.integer(particle.no),
                                                                            fluence.cm2                 =   as.single(fluence.cm2),
                                                                            material.no                 =   as.integer(material.no),
                                                                            er.model                    =   as.integer(er.model),
                                                                            u                           =   as.single(u))
                                                                            

    return(res$u)
}

##############################
AT.u      <-  function(   E.MeV.u,
                                                    particle.no,
                                                    fluence.cm2,
                                                    material.no,
                                                    er.model)
{
    n                           <-  length(E.MeV.u)
    u_cur                          <-  numeric(1)
    u	                          <-  numeric(n)
    for (i in 1:n){
		res                         <-  .C( "AT_total_u_R", n                           =   as.integer(1),
                                                                           E.MeV.u                     =   as.single(E.MeV.u[i]),
                                                                           particle.no                 =   as.integer(particle.no[i]),
                                                                           fluence.cm2                 =   as.single(fluence.cm2[i]),
                                                                            material.no                 =   as.integer(material.no),
                                                                            er.model                    =   as.integer(er.model),
                                                                            u                           =   as.single(u_cur))
                                                                            
		u[i]	<-	res$u
	}
    
	return(u)
}

########################################
AT.Bethe.Mass.Stopping.Power.MeV.cm2.g     <-  function(   E.MeV.u,
                                                    particle.no,
                                                    material.no,
                                                    E.restricted.keV)
{
    n                           						<-  length(E.MeV.u)
    Bethe.Mass.Stopping.Power.MeV.cm2.g                 <-  numeric(n)
	res                         <-  .C( "AT_Bethe_Mass_Stopping_Power_MeV_cm2_g_R", n                           		=   as.integer(n),
																					E.MeV.u                     		=   as.single(E.MeV.u),
																					particle.no                 		=   as.integer(particle.no),
																					material.no                 		=   as.integer(material.no),
																					E.restricted.keV            		=   as.single(E.restricted.keV),
																					Bethe.Mass.Stopping.Power.MeV.cm2.g =   as.single(Bethe.Mass.Stopping.Power.MeV.cm2.g))
                                                                            
	return(res$Bethe.Mass.Stopping.Power.MeV.cm2.g)
}

###############################
AT.GSM.calculate.dose.histogram	<-	function(	E.MeV.u,
									particle.no,
									fluence.cm2.or.dose.Gy,
									material.no,
									RDD.model,
									RDD.parameters,
									ER.model,
									nX,
									pixel.size.m,
									N.runs,
									dose.bin.centers.Gy){
	
		n				<-	length(E.MeV.u)
		n.bins			<-	length(dose.bin.centers.Gy)
		
		if(fluence.cm2.or.dose.Gy[1] < 0){
			fluence.cm2.or.dose.Gy	<-	AT.fluence.cm2(	E.MeV.u,
														particle.no,
														-1.0 * fluence.cm2.or.dose.Gy,
														material.no)	
		}
		
		dose.frequency.Gy		<-	numeric(n.bins)
		zero.dose.fraction	<-	numeric(1)

		res				<-	.C(	"AT_GSM_calculate_dose_histogram_R",	n						= 	as.integer(n),
																		E.MeV.u					=	as.single(E.MeV.u),
																		fluence.cm2.or.dose.Gy	=	as.single(fluence.cm2.or.dose.Gy),
																		particle.no				=	as.integer(particle.no),
																		material.no				=	as.integer(material.no),
																		RDD.model				=	as.integer(RDD.model),
																		RDD.parameters			=	as.single(RDD.parameters),
																		ER.model				=	as.integer(ER.model),
																		nX						=	as.integer(nX),
																		pixel.size.m			=	as.single(pixel.size.m),
																		N.runs				=	as.integer(N.runs),
																		number.of.bins			=	as.integer(n.bins),
																		dose.bin.centers.Gy		=	as.single(dose.bin.centers.Gy),
																		zero.dose.fraction		=	as.single(zero.dose.fraction),
																		dose.frequency.Gy		=	as.single(dose.frequency.Gy))


	results			<-	data.frame(		dose.bin.centers.Gy		= dose.bin.centers.Gy,
										dose.frequency.Gy		= res$dose.frequency.Gy,
										zero.dose.fraction		= rep(res$zero.dose.fraction, n.bins))

	
	return(results)
}

#########################################
AT.GSM.calculate.multiple.dose.histograms	<-	function(	E.MeV.u,
															particle.no,
															fluence.cm2.or.dose.Gy,
															material.no,
															RDD.model,
															RDD.parameters,
															ER.model,
															nX,
															pixel.size.m,
															N.runs,
															N.repetitions,
															dose.bin.centers.Gy){
	
		n				<-	length(E.MeV.u)
		n.bins			<-	length(dose.bin.centers.Gy)
		
		if(fluence.cm2.or.dose.Gy[1] < 0){
			fluence.cm2.or.dose.Gy	<-	AT.fluence.cm2(	E.MeV.u,
														particle.no,
														-1.0 * fluence.cm2.or.dose.Gy,
														material.no)	
		}
		
		dose.bin.widths.Gy		<-	numeric(n.bins)
		mean.dose.frequency.Gy	<-	numeric(n.bins)
		sd.dose.frequency.Gy	<-	numeric(n.bins)
		mean.zero.dose.fraction	<-	numeric(1)
		sd.zero.dose.fraction	<-	numeric(1)
		mean.d.check.Gy			<-	numeric(1)
		sd.d.check.Gy			<-	numeric(1)

		res				<-	.C(	"AT_GSM_calculate_multiple_dose_histograms_R",	n						= 	as.integer(n),
																				E.MeV.u					=	as.single(E.MeV.u),
																				fluence.cm2.or.dose.Gy	=	as.single(fluence.cm2.or.dose.Gy),
																				particle.no				=	as.integer(particle.no),
																				material.no				=	as.integer(material.no),
																				RDD.model				=	as.integer(RDD.model),
																				RDD.parameters			=	as.single(RDD.parameters),
																				ER.model				=	as.integer(ER.model),
																				nX						=	as.integer(nX),
																				pixel.size.m			=	as.single(pixel.size.m),
																				N.runs					=	as.integer(N.runs),
																				N.repetitions			=	as.integer(N.repetitions),
																				number.of.bins			=	as.integer(n.bins),
																				dose.bin.centers.Gy		=	as.single(dose.bin.centers.Gy),
																				dose.bin.widths.Gy		=	as.single(dose.bin.widths.Gy),
																				mean.d.check.Gy			=	as.single(mean.d.check.Gy),
																				sd.d.check.Gy			=	as.single(sd.d.check.Gy),
																				mean.zero.dose.fraction	=	as.single(mean.zero.dose.fraction),
																				sd.zero.dose.fraction	=	as.single(sd.zero.dose.fraction),
																				mean.dose.frequency.Gy	=	as.single(mean.dose.frequency.Gy),
																				sd.dose.frequency.Gy	=	as.single(sd.dose.frequency.Gy))


	results			<-	data.frame(		dose.bin.centers.Gy		= dose.bin.centers.Gy,
										dose.bin.width.Gy		= res$dose.bin.widths.Gy,
										mean.dose.frequency.Gy	= res$mean.dose.frequency.Gy,
										sd.dose.frequency.Gy	= res$sd.dose.frequency.Gy,
										mean.zero.dose.fraction	= rep(res$mean.zero.dose.fraction, n.bins),
										sd.zero.dose.fraction	= rep(res$sd.zero.dose.fraction, n.bins),
										mean.d.check.Gy			= rep(res$mean.d.check.Gy, n.bins),
										sd.d.check.Gy			= rep(res$sd.d.check.Gy, n.bins))
	
	return(results)
}


######################################################################################################################################

print(AT.version)
print(" ")
print("List of functions:")
print(sort(ls()[grep("AT.", ls())]))
print("end.")