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
################################################################################################

################################################################################################
# FUNCTION LIST
#
# :::GENERAL FUNCTIONS:::
# AT.beta.from.E                             wrapper needed
# AT.effective.charge.from.beta              wrapper needed
# AT.effective.charge.from.particle.no       wrapper needed
# AT.gamma.response                          looks OK, to be tested
# AT.max.E.transfer.MeV                      OK
# AT.particle.properties                     function missing				   
# AT.max.electron.range                      OK
# AT.D.Gy                                    wrapper needed
# AT.convert.beam.parameters                 wrapper needed
# AT.CSDA.range.g.cm2                        wrapper needed
# AT.Z.from.particle.no                      wrapper needed
# AT.A.from.particle.no                      wrapper needed
# AT.read.spectrum                           C function not implemented
# AT.scaled.energy.from.particle.no          wrapper needed
#
# :::MATERIAL FUNCTIONS:::
# AT.get.material.data                       wrapper needed
# AT.density.g.cm3                           wrapper needed
# AT.LET.MeV.cm2.g                           wrapper needed, should LET be in "material functions" ?
# AT.LET.keV.um                              wrapper needed, should LET be in "material functions" ?
# AT.electron.density.m3                     wrapper needed
# AT.E.MeV.u                                 wrapper needed, function should be renamed, should LET be in "material functions" ?
#
# ::: KATZ MODEL TEST FUNCTIONS:::
# will be implemented soon
#
# :::RADIAL DOSE FUNCTIONS:::
# AT.RDD.D.Gy                                OK
# AT.RDD.D.ext.Gy                            OK
# AT.RDD.r.m                                 OK
# AT.RDD.f1.parameters                       wrapper needed
#
# :::Successive convolution FUNCTIONS:::
# AT.SC.get.f1                               wrapper needed, function needs to be refactored
# AT.SC.get.gamma.response                   wrapper written, not used, why not use here AT.gamma.response ? 
# AT.SC.do.SC                                wrapper needed, function needs to be refactored
#
# :::Compute efficiency FUNCTIONS:::
# AT.GSM                                     wrapper needed
# AT.IGK                                     wrapper needed
# AT.SPIFF                                   wrapper needed
# AT.RBE                                     wrapper needed, function needs to be refactored
#
# AT.efficiency                              function needs to be refactored
# AT.SPIFF.short                             function needs to be refactored
# AT.fit.linquad.chi2                        needs to be implemented in C and interfaced
# AT.fit.linquad.chi2.grad                   needs to be implemented in C and interfaced
# AT.fit.linquad                             needs to be implemented in C and interfaced
# AT.Edmund.Transport                        needs to be implemented in C and interfaced
# E.delta                                    needs to be implemented in C and interfaced, is it needed ?
# 
################################################################################################

debug <- F

print("libamtrack S/R wrapping script - 2010/02/21")

##################
AT.beta.from.E	<-	function(	E.MeV.u ){
	n		<-	length(E.MeV.u)
	beta	<-	numeric(n)
	res		<-	.C(	"AT_beta_from_E",		n					=	as.integer(n),
												E.MeV.u			=	as.single(E.MeV.u),
												beta				=	as.single(beta))
	return(res$beta)						
}


##############################
# TODO function AT_density_g_cm3_from_material_no takes only single variable
AT.density.g.cm3		<-	function(	material.no){
	density.g.cm3			<-	numeric(1)
	res						<-	.C(	"AT_density_g_cm3_from_material_no",			material.no		=	as.character(material.no),
																		density.g.cm3		=	as.single(density.g.cm3))
	return(res$density.g.cm3)						
}

########################################
AT.effective.charge.from.particle.no	<-	function(	energy 				= E.MeV.u,
																particle.no		= particle.no){
	n					<-	length(energy)
	effective.charge	<-	numeric(n)
	res					<-	.C(	"AT_effective_charge_from_particle_no",	n					=	as.integer(n),
																					E.MeV.u			=	as.single(energy),
																					particle.no	=	as.integer(particle.no),
																					effective.charge	=	as.single(effective.charge))
	return(res$effective.charge)						
}

##############################
AT.effective.charge.from.beta	<-	function(	beta					= beta,
													Z						= Z){
	n					<-	length(beta)
	effective.charge	<-	numeric(n)
	res					<-	.C(	"AT_effective_charge_from_beta",	n					=	as.integer(n),
																		beta				=	as.single(beta),
																		Z					=	as.integer(Z),
																		effective.charge	=	as.single(effective.charge))
	return(res$effective.charge)						
}

##############################
# TODO function AT_electron_density_m3_from_material_no takes only single variable
AT.electron.density.m3	<-	function(	material.no){
	el.dens.m3			<-	numeric(1)
	res					<-	.C(	"AT_electron_density_m3_from_material_no",			material.no		=	as.integer(material.no),
																		el.dens.m3			=	as.single(el.dens.m3))
	return(res$el.dens.m3)						
}

##################
AT.gamma.response	<-	function(	d.Gy,
										gamma.model,
										gamma.parameters){
											
	n						<-	length(d.Gy)
	S						<-	numeric(n)
	res					<-	.C(	"AT_gamma_response_R",		n					= as.integer(n),
																d.Gy				= as.single(d.Gy),
																gamma.model		= as.integer(gamma.model),
																gamma.parameters	= as.single(gamma.parameters),
																S					= as.single(S))
	return(res$S)
}

##############################
AT.get.material.data	<-	function(	material.no){
	density.g.cm3		<-	numeric(1)
	el.dens.m3			<-	numeric(1)
	I.eV				<-	numeric(1)
	alpha.g.cm2.MeV	<-	numeric(1)
	p.MeV				<-	numeric(1)
	m.g.cm2			<-	numeric(1)
	n					<-	1
	
# average A and Z needs to be added
	
	res					<-	.C(	"AT_get_materials_data",				n					=	as.integer(n),
																		material.no		=	as.integer(material.no),
																		density.g.cm3		=	as.single(density.g.cm3),
																		el.dens.m3			=	as.single(el.dens.m3),
																		I.eV				=	as.single(I.eV),
																		alpha.g.cm2.MeV	=	as.single(alpha.g.cm2.MeV),
																		p.MeV				=	as.single(p.MeV),
																		m.g.cm2			=	as.single(m.g.cm2))
	return(res)						
}

#################
AT.D.Gy		<-	function(	E.MeV.u,
							particle.no,
							fluence.cm2,
							material.no){
	n					<-	length(E.MeV.u)
	D.Gy				<-	numeric(n)
	res					<-	.C(	"AT_D_Gy",		n					=	as.integer(n),
												E.MeV.u				=	as.single(E.MeV.u),
												particle.no			=	as.integer(particle.no),
												fluence.cm2			=	as.single(fluence.cm2),
												material.no			=	as.integer(material.no),
												D.Gy				=	as.single(D.Gy))
	return(res$D.Gy)						
}

#################
AT.LET.MeV.cm2.g		<-	function(	E.MeV.u,
										particle.no,
										material.no){
	n					<-	length(E.MeV.u)
	LET.MeV.cm2.g		<-	numeric(n)
	res					<-	.C(	"AT_LET_MeV_cm2_g",		n						=	as.integer(n),
															E.MeV.u				=	as.single(E.MeV.u),
															particle.no			=	as.integer(particle.no),
															material.no			=	as.integer(material.no),
															LET.MeV.cm2.g			=	as.single(LET.MeV.cm2.g))
	return(res$LET.MeV.cm2.g)						
}

##############
AT.LET.keV.um	<-	function(	E.MeV.u,
									particle.no,
									material.no){
	n					<-	length(E.MeV.u)
	LET.keV.um			<-	numeric(n)
	res					<-	.C(	"AT_LET_keV_um",			n						=	as.integer(n),
															E.MeV.u				=	as.single(E.MeV.u),
															particle.no			=	as.integer(particle.no),
															material.no			=	as.integer(material.no),
															LET.keV.um				=	as.single(LET.keV.um))
	return(res$LET.keV.um)						
}

###########
AT.E.MeV.u			<-	function(	LET.MeV.cm2.g,
										particle.no,
										material.no){
	n					<-	length(LET.MeV.cm2.g)
	E.MeV.u			<-	numeric(n)
	res					<-	.C(	"AT_E_MeV_from_LET",	n						=	as.integer(n),
															LET.MeV.cm2.g			=	as.single(LET.MeV.cm2.g),
															particle.no			=	as.integer(particle.no),
															material.no			=	as.integer(material.no),
															E.MeV.u				=	as.single(E.MeV.u))
	return(res$E.MeV.u)						
}


####################
AT.CSDA.range.g.cm2	<-	function(		E.MeV.u,
											particle.no,
											material.no){
	n					<-	length(E.MeV.u)
	CSDA.range.g.cm2	<-	numeric(n)
	res					<-	.C(	"AT_CSDA_range_g_cm2",		n						=	as.integer(n),
																E.MeV.u				=	as.single(E.MeV.u),
																particle.no			=	as.integer(particle.no),
																material.no			=	as.integer(material.no),
																CSDA.range.g.cm2		=	as.single(CSDA.range.g.cm2))
	return(res$CSDA.range.g.cm2)						
}


######################
AT.max.E.transfer.MeV		<-	function(	E.MeV.u ){
	n							<-	length(E.MeV.u)
	max.E.transfer.MeV		<-	numeric(n)
	res							<-	.C(	"AT_max_E_transfer_MeV_R",	n						=	as.integer(n),
																		E.MeV.u				=	as.single(E.MeV.u),
																		max.E.transfer.MeV	=	as.single(max.E.transfer.MeV))
	return(res$max.E.transfer.MeV)						
}


######################
AT.Z.from.particle.no	<-	function( particle.no){
	n						<-	length(particle.no)
	Z						<-	integer(n)

	res						<-	.C(	"AT_Z_from_particle_no",	n				= as.integer(n),
																particle.no		= as.integer(particle.no),
																Z				= as.integer(Z))
	return(res$Z)	
}


######################
AT.A.from.particle.no	<-	function( particle.no){
	n						<-	length(particle.no)
	A						<-	integer(n)

	res						<-	.C(	"AT_A_from_particle_no",	n				= as.integer(n),
																particle.no		= as.integer(particle.no),
																A				= as.integer(A))
	return(res$A)	
}

#################
AT.read.spectrum		<-	function(	file.name){
	nLines					<-	numeric(1)
	res						<-	.C(	"AT_browseSpectrumS",	file.name			= as.character(file.name),
																nLines				= as.integer(nLines))
	
	nLines					<-	res$nLines
	
	E.MeV.u				<-	numeric(nLines)
	particle.no		<-	integer(nLines)
	fluence.cm2			<-	numeric(nLines)
	slab.no				<-	integer(nLines)
	
	res						<-	.C(	"AT_readSpectrumS",		file.name			= as.character(file.name),
																nLines				= as.integer(nLines),
																E.MeV.u			= as.single(E.MeV.u),
																particle.no	= as.integer(particle.no),
																fluence.cm2		= as.single(fluence.cm2),
																slab.no			= as.integer(slab.no))
	return(				data.frame(	E.MeV.u			=	res$E.MeV.u,
											particle.no	=	res$particle.no,
											fluence.cm2		=	res$fluence.cm2,
											slab.no			=	res$slab.no))
}

#####################################
AT.scaled.energy.from.particle.no	<-	function(	energy 				= E.MeV.u,
															particle.no		= particle.no){
	n					<-	length(energy)
	scaled.energy		<-	numeric(n)
	res					<-	.C(	"AT_scaled_energy",								n					=	as.integer(n),
																					E.MeV.u			=	as.single(energy),
																					particle.no	=	as.integer(particle.no),
																					scaled.energy		=	as.single(scaled.energy))
	return(res$scaled.energy)						
}

#############################
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
	
	res					<-	.C(	"AT_convert_beam_parameters",						n					=	as.integer(n),
																					fluence.cm2			=	as.single(fluence.cm2),
																					sigma.cm			=	as.single(sigma.cm),
																					N					=	as.single(N),
																					FWHM.mm				=	as.single(FWHM.mm))
	return(res)						
}


################################################################################################
# Max delta electron range
################################################################################################

############
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

################################################################################################
# Katz model test functions
################################################################################################

# Katz model test functions:	

# TODO to be implemented

################################################################################################
# RDD functions
################################################################################################

############
AT.RDD.D.Gy					<-	function(	r.m,
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
														E.MeV.u				=	as.single(E.MeV.u),
														particle.no			=	as.integer(particle.no),
														material.no			=	as.integer(material.no),
														RDD.model				=	as.integer(RDD.model),
														RDD.parameters		=	as.single(RDD.parameters),
														ER.model				=	as.integer(ER.model),
														ER.parameters			=	as.single(ER.parameters),
														D.Gy					=	as.single(D.Gy))		
			
	 return(res$D.Gy)						
}

###########
AT.RDD.r.m					<-	function(	D.Gy,
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


#####################
AT.RDD.f1.parameters	<-	function(	E.MeV.u,
											particle.no,
											material.no,
											RDD.model,
											RDD.parameters,
											ER.model,
											ER.parameters){

	n								<-	length(E.MeV.u)
	f1.parameters					<-	data.frame(	LET.MeV.cm2.g						= numeric(n),
														r.min.m							= numeric(n),
														r.max.m							= numeric(n),
														d.min.Gy							= numeric(n),
														d.max.Gy							= numeric(n),
														k.Gy								= numeric(n),
														single.impact.fluence.cm2		= numeric(n),
														single.impact.dose.Gy			= numeric(n),
														dEdx.MeV.cm2.g					= numeric(n))
	f1.tmp							<-	numeric(9)
														
	for (i in 1:n){
		#i<-1
		res				<-	.C("AT_RDD_f1_parameters", 	E.MeV.u						= as.single(E.MeV.u[i]),
																particle.no					= as.integer(particle.no[i]),
																material.no					= as.integer(material.no),
																RDD.model						= as.integer(RDD.model),
																RDD.parameters				= as.single(RDD.parameters),
																ER.model						= as.integer(ER.model),
																ER.parameters					= as.single(ER.parameters),
																f1.tmp							= as.single(f1.tmp))
		f1.parameters[i,]		<-	res$f1.tmp
	}
	
	return(f1.parameters)
}
										
################################################################################################
# Successive convolutions functions
################################################################################################

############
AT.SC.do.SC	<-	function(	f.list,
								fluence.factor = 1.0,
								write.output = F,
								shrink.tails = T,
								shrink.tails.under = 1e-30,
								adjust.N2 = T){
	n.bins.f				<-	integer(1)
	u.start				<-	numeric(1)
	n.convolutions		<-	integer(1)

	if(debug == T) cat("1 n.convolutions = ",n.convolutions,"\n")
	if(debug == T) cat("1 n.bins.f1 = ",nrow(f.list$f1),"\n")
	if(debug == T) cat("1 f1.dd.Gy = ",f.list$f1$f1.dd.Gy,"\n")
	if(debug == T) cat("1 f1 = ",f.list$f1$f1,"\n")
		
	res						<-	.C("AT_SC_get_f_array_size", 	u						= as.single(f.list$f.parameters[1]),
																	fluence.factor		= as.single(fluence.factor),
																	N2						= as.integer(unique(f.list$f1$N2)),
																	n.bins.f1				= as.integer(nrow(f.list$f1)),
																	f1.d.Gy				= as.single(f.list$f1$f1.d.Gy),
																	f1.dd.Gy				= as.single(f.list$f1$f1.dd.Gy),
																	f1						= as.single(f.list$f1$f1),
																	n.bins.f				= as.integer(n.bins.f),
																	u.start				= as.single(u.start),
																	n.convolutions		= as.integer(n.convolutions))
	n.bins.f				<-	res$n.bins.f
	u.start				<-	res$u.start
	n.convolutions		<-	res$n.convolutions
	
	if(debug == T) cat("2 n.convolutions = ",n.convolutions,"\n")
	
	f.list$f				<-	data.frame(	f.d.Gy			= numeric(n.bins.f),
												f.dd.Gy		= numeric(n.bins.f),
												f				= numeric(n.bins.f),
												n.bins.f		= rep(n.bins.f, n.bins.f),
												N2				= rep(unique(f.list$f1$N2), n.bins.f),
												u.start		= rep(u.start, n.bins.f))

	if(debug == T) cat("3\n")
	if(debug == T) cat("3 n.bins.f1 = ",nrow(f.list$f1),"\n")
	if(debug == T) cat("3 unique(f.list$f1$u) = ",unique(unique(f.list$f.parameters[1])),"\n")
	if(debug == T) cat("3 unique(f.list$f1$N2) = ",unique(f.list$f1$N2),"\n")
												
	res						<-	.C("AT_SC_get_f_start", 		u						= as.single(unique(f.list$f.parameters[1])),
																	n.bins.f1				= as.integer(nrow(f.list$f1)),
																	N2						= as.integer(unique(f.list$f1$N2)),
																	f1.d.Gy				= as.single(f.list$f1$f1.d.Gy),
																	f1.dd.Gy				= as.single(f.list$f1$f1.dd.Gy),
																	f1						= as.single(f.list$f1$f1),
																	n.bins.f				= as.integer(n.bins.f),
																	f.d.Gy					= as.single(f.list$f$f.d.Gy),
																	f.dd.Gy				= as.single(f.list$f$f.dd.Gy),
																	f						= as.single(f.list$f$f))
	
	f.list$f$f.d.Gy				<- res$f.d.Gy
	f.list$f$f.dd.Gy				<- res$f.dd.Gy
	f.list$f$f						<- res$f
	
	f.list$f.start				<-	f.list$f
	
	f.list$f$fdd					<-	numeric(nrow(f.list$f))
	f.list$f$dfdd					<-	numeric(nrow(f.list$f))
	
	f0						<- numeric(1)
	d						<- numeric(1)
	
	if(debug == T) cat("5\n")

	res						<-	.C("AT_SuccessiveConvolutions", 	u						= as.single(unique(f.list$f.parameters[1])),
																		n.bins.f				= as.integer(n.bins.f),
																		N2						= as.integer(unique(f.list$f1$N2)),
																		n.bins.f.used			= as.integer(nrow(f.list$f1)),
																		f.d.Gy					= as.single(f.list$f$f.d.Gy),
																		f.dd.Gy				= as.single(f.list$f$f.dd.Gy),
																		f						= as.single(f.list$f$f),
																		f0						= as.single(f0),
																		fdd						= as.single(f.list$f$fdd),
																		dfdd					= as.single(f.list$f$dfdd),
																		d						= as.single(d),
																		write.output			= as.integer(write.output),
																		shrink.tails			= as.integer(shrink.tails),
																		shrink.tails.under	= as.single(shrink.tails.under),
																		adjust.N2				= as.integer(adjust.N2))
																		
	if(debug == T) cat("6\n")
																		
	f.list$f$f.d.Gy				<-	res$f.d.Gy
	f.list$f$f.dd.Gy				<-	res$f.dd.Gy
	f.list$f$f						<-	res$f

	if(debug == T) cat("6a\n")

	f.list$f$fdd					<-	res$fdd
	f.list$f$dfdd					<-	res$dfdd
	f.list$f$f0					<-	rep(res$f0, n.bins.f)
	f.list$f$d						<-	rep(res$d, n.bins.f)
		if(debug == T) cat("6b\n")

	f.list$f$u.start				<-	rep(u.start, n.bins.f)
	f.list$f$u						<-	rep(res$u, n.bins.f)
	f.list$f$n.convolutions		<-	rep(n.convolutions, n.bins.f)
	f.list$f$n.bins.f.used		<-	rep(res$n.bins.f.used, n.bins.f)
	f.list$f$write.output		<-	rep(write.output, n.bins.f)
	f.list$f$shrink.tails		<-	rep(shrink.tails, n.bins.f)
	
		if(debug == T) cat("6c\n")
	
	f.list$f$shrink.tails.under	<-	rep(shrink.tails.under, n.bins.f)
	f.list$f$adjust.N2			<-	rep(adjust.N2, n.bins.f)
	f.list$f$N2.adjusted			<-	rep(res$N2, n.bins.f)
	tmp <- res$n.bins.f.used
		if(debug == T) cat("6d",tmp,"\n")
	
	f.list$f						<-	f.list$f[1:res$n.bins.f.used,]

	if(debug == T) cat("7\n")

	return(f.list)
}


#############
AT.SC.get.f1		<-	function(		E.MeV.u,
										particle.no,
										fluence.cm2,
										material.no,
										RDD.model,
										RDD.parameters,
										ER.model,
										ER.parameters,
										N2){
											
	#####################################################
	# FIRST: get array size for single impact function f1
	#####################################################

	n				<-	length(E.MeV.u)
	f1.parameters	<-	numeric(9 * n)
	n.bins.f1 <- 0
	
	res				<-	.C("AT_SC_get_f1_array_size", 	n						= as.integer(n),
																E.MeV.u				= as.single(E.MeV.u),
																particle.no			= as.integer(particle.no),
																material.no			= as.integer(material.no),
																RDD.model				= as.integer(RDD.model),
																RDD.parameters			= as.single(RDD.parameters),
																ER.model				= as.integer(ER.model),
																ER.parameters			= as.single(ER.parameters),
																N2						= as.integer(N2),
																n.bins.f1				= as.integer(n.bins.f1),
																f1.parameters			= as.single(f1.parameters))
		
	n.bins.f1		<- 	res$n.bins.f1
	f1.parameters	<-	res$f1.parameters
	
	#####################################################
	# THEN: get f1 and f1 parameters
	#####################################################

	f.parameters					<-	numeric(7)
	
	norm.fluence					<-	numeric(n)
	dose.contribution.Gy			<-	numeric(n)
	
	f1.d.Gy						<-	numeric(n.bins.f1)
	f1.dd.Gy						<-	numeric(n.bins.f1)
	f1								<-	numeric(n.bins.f1)	

	res				<-	.C("AT_SC_get_f1", 					n						= as.integer(n),
																E.MeV.u				= as.single(E.MeV.u),
																particle.no			= as.integer(particle.no),
																fluence.cm2			= as.single(fluence.cm2),
																material.no			= as.integer(material.no),
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
																f1.d.Gy				= as.single(f1.d.Gy),
																f1.dd.Gy				= as.single(f1.dd.Gy),
																f1						= as.single(f1))
	
	results	<-	list(	general				=	data.frame(	E.MeV.u						=	res$E.MeV.u,
																		particle.no					= 	res$particle.no,
																		fluence.cm2					= 	res$fluence.cm2,
																		fluence.cm2.set				= 	fluence.cm2,
																		norm.fluence					= 	res$norm.fluence,
																		dose.contribution.Gy			= 	res$dose.contribution.Gy,
																		material.no					=	rep(material.no, n)),
																		f1.parameters			=	f1.parameters,											
							           f.parameters			=	res$f.parameters,
           							f1						=	data.frame(	f1.d.Gy						=	res$f1.d.Gy,
																		f1.dd.Gy						=	res$f1.dd.Gy,
																		f1								=	res$f1,
																		material.no					=	rep(material.no, n.bins.f1),
																		n.bins.f1						=	rep(res$n.bins.f1, n.bins.f1),
																		N2								=	rep(N2, n.bins.f1)))
	return(results)
}


#########################
AT.SC.get.gamma.response		<-	function(		f.list,
														gamma.model,
														gamma.parameters){
	n						<-	nrow(f.list$f)
	S						<-	numeric(n)
	S.HCP					<-	numeric(1)
	S.gamma				<-	numeric(1)
	efficiency				<-	numeric(1)

	n.gamma.parameters 	<-	length(gamma.parameters)
	if(gamma.model == 1){
		n.gamma.parameters 	<-	length(gamma.parameters) / 4
	}

	res					<-	.C(	"AT_get_gamma_response",		n						= as.integer(n),
																	f.d.Gy					= as.single(f.list$f$f.d.Gy),
																	f.dd.Gy				= as.single(f.list$f$f.dd.Gy),
																	f						= as.single(f.list$f$f),
																	f0						= as.single(unique(f.list$f$f0)),
																	gamma.model			= as.integer(gamma.model),
																	gamma.parameters		= as.single(gamma.parameters),
																	S						= as.single(S),
																	S.HCP					= as.single(S.HCP),
																	S.gamma				= as.single(S.gamma),
																	efficiency				= as.single(efficiency))

	f.list$f$n.gamma.parameters			<-	rep(n.gamma.parameters, n)
	f.list$f$gamma.model					<-	rep(gamma.model, n)
	
	# create one-row data-frame from gamma parameters
	
	# This is not working correctly in R
#	df.gamma				<-	t(as.data.frame(gamma.parameters))
#	rownames(df.gamma)	<-	""
#	colnames(df.gamma)	<-	paste("gamma.parameters.", 1:ncol(df.gamma), sep = "")
#	df.gamma				<-	df.gamma[rep(1, n),]
#	f.list$f								<-	cbind.data.frame(f.list$f, df.gamma)

	f.list$f$S				<-	res$S
	f.list$f$S.HCP		<-	rep(res$S.HCP, n)
	f.list$f$S.gamma		<-	rep(res$S.gamma, n)
	f.list$f$efficiency	<-	rep(res$efficiency, n)
	
	return(f.list)
}


##############
AT.efficiency		<-	function(		E.MeV.u,
											particle.no,
											fluence.cm2,
											material.no,
											RDD.model,
											RDD.parameters,
											ER.model,
											ER.parameters,
											gamma.model,
											gamma.parameters,
											method					= "SC",
											N2						= 20,
											fluence.factor 		= 1.0,
											write.output 			= F,
											shrink.tails 			= T,
											shrink.tails.under	= 1e-30,
											adjust.N2 				= T){
	
	debug 				<-	F
	f.list				<-	AT.SC.get.f1(					E.MeV.u,
																particle.no,
																fluence.cm2,
																material.no,
																RDD.model,
																RDD.parameters,
																ER.model,
																ER.parameters,
																N2)
	
	f.list				<-	AT.SC.do.SC( f.list)

	f.list				<-	AT.SC.get.gamma.response(		f.list,
																gamma.model,
																gamma.parameters)
														
	f.list$results	<-	data.frame(	efficiency  	=	unique(f.list$f$efficiency),
											d.check.Gy		=	unique(f.list$f$d))

	return(f.list)
}


AT.SPIFF.short	<-	function(	E.MeV.u,
										particle.no,
										fluence.cm2,
										material.no,
										RDD.model,
										RDD.parameters,
										ER.model,
										ER.parameters,
										gamma.model,
										gamma.parameters,
										method					= "SC",
										N2						= 20,
										fluence.factor 		= 1.0,
										write.output 			= F,
										shrink.tails 			= T,
										shrink.tails.under	= 1e-30,
										adjust.N2 				= T,
          lethal.events.mode = F,									
										transport				= F,
										E.min.MeV.u			= 1,
										p.tol					= 0.1,
										z.total.cm				= 0.1,
					 				N.max					= 1e3){
	
		results			<-	numeric(10)
		if (transport == T){
		transport		<-	AT.Edmund.Transport(	E.start.MeV.u		= E.MeV.u,
														E.min.MeV.u		= E.min.MeV.u,
														particle.no		= particle.no,
														material.no		= material.no,
														p.tol				= p.tol,
														z.total.cm			= z.total.cm,
					 						  			N.max				= N.max)
		print(transport)
		n.loop			<-	nrow(transport)
		for (i in 1:n.loop){
			slab.results		<-	numeric(10)
			res					<-	.C(	"AT_SPIFF",		n						= 	1,
																E.MeV.u				=	as.single(transport$E.MeV[i]),
																particle.no			=	as.integer(particle.no),
																fluence.cm2			=	as.single(fluence.cm2),
																material.no			=	as.integer(material.no),
																RDD.model				=	as.integer(RDD.model),
																RDD.parameters		=	as.single(RDD.parameters),
																ER.model				=	as.integer(ER.model),
																ER.parameters			=	as.single(ER.parameters),
																gamma.model			=	as.integer(gamma.model),
																gamma.parameters		=	as.single(gamma.parameters),
																N2						=	as.integer(N2),
																fluence.factor		=	as.single(transport$D.rel[i]),					# Bragg curve passed by fluence factor
																write.output			=	as.integer(write.output),
																shrink.tails			=	as.integer(shrink.tails),
																shrink.tails.under	=	as.single(shrink.tails.under),
																adjust.N2				=	as.integer(adjust.N2),
                lethal.events.mode = as.integer(lethal.events.mode),
																results				=	as.single(slab.results))
			print(cat("Done slab no. ", i, "\n"))
			print(cat("results:  ", res$results, "\n"))
			results[1]			<-	results[1] + res$results[1] * transport$w.Edmund[i]			# efficiency
			results[2]			<-	results[2] + res$results[2] * transport$w.Edmund[i]			# d.check
			results[3]			<-	results[3] + res$results[3] * transport$w.Edmund[i]			# S.HCP
			results[4]			<-	results[4] + res$results[4] * transport$w.Edmund[i]			# S.gamma
		}
	}else{
		res					<-	.C(	"AT_SPIFF",		n						= 	as.integer(length(E.MeV.u)),
															E.MeV.u				=	as.single(E.MeV.u),
															particle.no			=	as.integer(particle.no),
															fluence.cm2			=	as.single(fluence.cm2),
															material.no			=	as.integer(material.no),
															RDD.model				=	as.integer(RDD.model),
															RDD.parameters		=	as.single(RDD.parameters),
															ER.model				=	as.integer(ER.model),
															ER.parameters			=	as.single(ER.parameters),
															gamma.model			=	as.integer(gamma.model),
															gamma.parameters		=	as.single(gamma.parameters),
															N2						=	as.integer(N2),
															fluence.factor		=	as.single(fluence.factor),
															write.output			=	as.integer(write.output),
															shrink.tails			=	as.integer(shrink.tails),
															shrink.tails.under	=	as.single(shrink.tails.under),
															adjust.N2				=	as.integer(adjust.N2),
															lethal.events.mode = as.integer(lethal.events.mode),
															results				=	as.single(results))
		results			<-	res$results
	}
	return(results)
}
														
AT.GSM	<-	function(	E.MeV.u,
										particle.no,
										fluence.cm2,
										material.no,
										RDD.model,
										RDD.parameters,
										ER.model,
										ER.parameters,
										gamma.model,
										gamma.parameters,
										method					= "grid",
										N.runs					=	100,
										N2					=	10,
										fluence.factor		=	1.0,
										write.output			=	T,
										n.X						=	100,
										lethal.events.mode  = F,
										grid.size.m			=	1e-9){
	
		results			<-	numeric(10)
		res					<-	.C(	"AT_GSM",		n				 = 	as.integer(length(E.MeV.u)),
															E.MeV.u			    	    =	as.single(E.MeV.u),
															particle.no			     =	as.integer(particle.no),
															fluence.cm2			     =	as.single(fluence.cm2),
															material.no			     =	as.integer(material.no),
															RDD.model				      =	as.integer(RDD.model),
															RDD.parameters    	=	as.single(RDD.parameters),
															ER.model    				   =	as.integer(ER.model),
															ER.parameters		    =	as.single(ER.parameters),
															gamma.model			     =	as.integer(gamma.model),
															gamma.parameters	  =	as.single(gamma.parameters),
															N.runs			        		=	as.integer(N.runs),
															N2      		      			=	as.integer(N2),
															fluence.factor		   =	as.single(fluence.factor),
															write.output		    	=	as.integer(write.output),
															n.X			          			=	as.integer(n.X),
															grid.size.m			     =	as.single(grid.size.m),
															lethal.events.mode = as.integer(lethal.events.mode),
															results		        		=	as.single(results))
	results			<-	res$results
	return(results)
}

AT.IGK	<-	function(		E.MeV.u,
										particle.no,
										fluence.cm2,
										material.no,
										RDD.model,
										RDD.parameters,
										ER.model,
										ER.parameters,
										gamma.model,
										gamma.parameters,
										s0.factor){
	
	 	results			<-	numeric(10)
		res					<-	.C(	"AT_IGK",		n						= 	as.integer(length(E.MeV.u)),
													E.MeV.u				=	as.single(E.MeV.u),
													particle.no			=	as.integer(particle.no),
													fluence.cm2			=	as.single(fluence.cm2),
													material.no			=	as.integer(material.no),
													RDD.model				=	as.integer(RDD.model),
													RDD.parameters		=	as.single(RDD.parameters),
													ER.model				=	as.integer(ER.model),
													ER.parameters			=	as.single(ER.parameters),
													gamma.model			=	as.integer(gamma.model),
													gamma.parameters		=	as.single(gamma.parameters),
													s0.factor				=	as.single(s0.factor),
													results				=	as.single(results))
	results			<-	res$results
	return(results)
}


AT.fit.linquad.chi2 <- function( params, D.Gy, 
															Survival){
	alpha <- params[1]
	beta <- params[2]
 res <- 0
 tmp <- (log(Survival/100.0) + D.Gy*alpha + D.Gy*D.Gy*beta)
 res <- sum(tmp*tmp)
	return (res)										 
}

AT.fit.linquad.chi2.grad <- function( params, D.Gy, 
															Survival){
	alpha <- params[1]
	beta <- params[2]
 res <- c(0., 0.)
 tmp <- (log(Survival/100.0) + D.Gy*alpha + D.Gy*D.Gy*beta)
 res[1] <- sum(2.*tmp*D.Gy)
 res[2] <- sum(2.*tmp*D.Gy*D.Gy)
	return (res)										 
}

AT.fit.linquad <- function( D.Gy, 
															Survival,
															Survival.min,
															Survival.max,
															alpha.zero,
															beta.zero){
 params <- c(0.,0.)
 params <- optim( par = c(0.1,0.01), fn = AT.fit.linquad.chi2, gr = AT.fit.linquad.chi2.grad, method = "L-BFGS-B", lower = c(0.,0.), upper = c(1000.,1000.), D.Gy = D.Gy, Survival = Survival)
	return (params[[1]])										 
}


AT.RBE <- function( alpha.ion, 
															beta.ion,
															alpha.gamma,
															beta.gamma,
															survival = 10){
 RBE <- length(alpha.ion)
 if( beta.ion == 0 ){
   D.ion <- -log( survival/100.0 ) / alpha.ion
	} else {
		 delta.ion <- alpha.ion*alpha.ion - 4. * beta.ion * log( survival / 100.0 )
		 D.ion <- (-alpha.ion + sqrt( delta.ion ) ) / (2.*beta.ion)	
	}  
 if( beta.gamma == 0 ){
   D.gamma <- -log( survival/100.0 ) / alpha.gamma
	} else {
   delta.gamma <- alpha.gamma*alpha.gamma - 4. * beta.gamma * log( survival / 100.0 )
   D.gamma <- (-alpha.gamma + sqrt( delta.gamma ) ) / (2.*beta.gamma)  
 }
 RBE <- D.gamma / D.ion
	return (RBE)										 
}

AT.Edmund.Transport	<- 	function(			E.start.MeV.u,
												E.min.MeV.u		=	1e-2,
												particle.no		=	1,				# proton
												material.no		=	1,				# water
												p.tol				=	0.05,
												z.total.cm			=	0.1,
					 						  	N.max				=	1e3){
	
	# This function keeps track of the energy account
	# as the proton penetrats the crystal of thickness.
	# E.start		= entrance proton energy 		[MeV]
	# p.tol		= energy tolerance 				[%]
	# LET.tab		= table to look up LET values	[MeV*cm^2/g]
	# names should be "MeV","LET,""R.max"
	# z.total.cm		= thickness of crystal			[cm]	
	# N.max		= maximum number of iterations [dimless]
	# dens			= material density				[g/cm^3]
	# E.min		= minimum proton energy			[MeV]
	
	# Created 08.27.06 by JE as ### Ez.account.1 ###
	# Rewritten for use with SGP library, 10.06.09, Steffen
	
	
	# Start slap, initial values
	HCP.properties	<-	AT.particle.properties(	particle.no)
	E.start.MeV		<-	E.start.MeV.u * HCP.properties$A
	E.min.MeV			<-	E.min.MeV.u * HCP.properties$A
	mat.properties	<-	AT.get.material.data(	material.no)
	LET.MeV.cm			<- 	AT.LET.MeV.cm2.g(	E.MeV.u			= E.start.MeV.u,
													particle.no		= particle.no,
													material.no		= material.no) * mat.properties$density.g.cm3
													
	dE.MeV				<- E.start.MeV * p.tol
		
	dz.cm 				<- dE.MeV / LET.MeV.cm					# Approximation: LET constant in dE



	E.MeV				<- E.start.MeV
	z.cm				<- 0
	E.dep.MeV			<- dE.MeV
	i 					<- 1

	df					<-	data.frame(	i				=	i,
											E.MeV			=	E.MeV,
											LET.MeV.cm		=	LET.MeV.cm,
											z.cm			=	z.cm,
											dz.cm			=	dz.cm,
											dE.MeV			=	dE.MeV,
											E.dep.MeV		=	E.dep.MeV,
											D.rel			=	1,
											w.Edmund		=	dz.cm / z.total.cm)												# weight for signal from this slab by relative mass: 
											
	
	# Iteration
	
	repeat{
		
		# dz.cm >> z.total.cm (large energies or small crystals)
		if (dz.cm >= z.total.cm){
			dz.cm 		<- z.total.cm - z.cm
			dE.MeV		<- dz.cm * LET.MeV.cm
			E.dep.MeV	<- dE.MeV

			df			<-	data.frame(	i			=	i,
											E.MeV		=	E.MeV,
											LET.MeV.cm	=	LET.MeV.cm,
											z.cm		=	z.cm,
											dz.cm		=	dz.cm,
											dE.MeV		=	dE.MeV,
											E.dep.MeV	=	E.dep.MeV,
											D.rel		=	dE.MeV / dz.cm / (df$dE.MeV[1] / df$dz.cm[1]),
											w.Edmund	=	dz.cm / z.total.cm)
			
			print("Only one energy deposition (large energy/small crystal)")
			break			
		}
		
		# dE.MeV << E.min.MeV (large crystals or small energies)
		if (E.MeV - dE.MeV <= E.min.MeV){
			dE.MeV			<- E.MeV - E.min.MeV
			dz.cm 			<- dE.MeV / LET.MeV.cm
			E.dep.MeV		<- dE.MeV

			df.add		<-	data.frame(	i			=	i,
											E.MeV		=	E.MeV,
											LET.MeV.cm	=	LET.MeV.cm,
											z.cm		=	z.cm,
											dz.cm		=	dz.cm,
											dE.MeV		=	dE.MeV,
											E.dep.MeV	=	E.dep.MeV,
											D.rel		=	dE.MeV / dz.cm / (df$dE.MeV[1] / df$dz.cm[1]),
											w.Edmund	=	dz.cm / z.total.cm)
			df			<-	rbind.data.frame(df, df.add)

			print("Only one energy deposition (small energy/large crystal)")
			break			
		}

		# Maximum number of iterations
		if (i > N.max){
			print ("Max number of iterations reached")
			break
		}
		
		# Normal iteration
		E.MeV			<- E.MeV - dE.MeV 									# E for next slap
		E.MeV.u		<- E.MeV / HCP.properties$A
		z.cm 			<- z.cm + dz.cm	 									# start depth of next slap

		LET.MeV.cm		<- 	AT.LET.MeV.cm2.g(	E.MeV.u			= E.MeV.u,
													particle.no	 	= particle.no,
													material.no		= material.no) * mat.properties$density.g.cm3

		dE.MeV			<- E.MeV * p.tol				 						# dE next slap
		dz.cm			<- dE.MeV / LET.MeV.cm								# size next slap
		E.dep.MeV		<- E.dep.MeV + dE.MeV								# accumulated energy
		i 				<- i+1													# number of iterations
		print(paste("Slab", i))

		# Last slap
		if (z.cm + dz.cm >= z.total.cm){
			dz.cm 			<- z.total.cm - z.cm
			dE.MeV			<- dz.cm * LET.MeV.cm
			E.dep.MeV 		<- E.dep.MeV + dE.MeV

			df.add			<-	data.frame(	i			=	i,
												E.MeV		=	E.MeV,
												LET.MeV.cm	=	LET.MeV.cm,
												z.cm		=	z.cm,
												dz.cm		=	dz.cm,
												dE.MeV		=	dE.MeV,
												E.dep.MeV	=	E.dep.MeV,
												D.rel		=	dE.MeV / dz.cm / (df$dE.MeV[1] / df$dz.cm[1]),
												w.Edmund	=	dz.cm / z.total.cm)
			df				<-	rbind.data.frame(df, df.add)

			print("Crystal thickness reached")
			break 
		}
		
		# Protons thermalized
		if (E.MeV - dE.MeV <= E.min.MeV){
			dE.MeV 			<- E.MeV - E.min.MeV
			dz.cm				<- dE.MeV / LET.MeV.cm
			z.cm				<- z.cm + dz.cm
			E.dep.MeV			<- E.dep.MeV + dE.MeV

			# Update vectors

			df.add			<-	data.frame(	i			=	i,
												E.MeV		=	E.MeV,
												LET.MeV.cm	=	LET.MeV.cm,
												z.cm		=	z.cm,
												dz.cm		=	dz.cm,
												dE.MeV		=	dE.MeV,
												E.dep.MeV	=	E.dep.MeV,
												D.rel		=	dE.MeV / dz.cm / (df$dE.MeV[1] / df$dz.cm[1]),
												w.Edmund	=	dz.cm / z.total.cm)
			df				<-	rbind.data.frame(df, df.add)

			print("Minimum energy reached")
			break 

		}
		
		df.add			<-	data.frame(	i			=	i,
											E.MeV		=	E.MeV,
											LET.MeV.cm	=	LET.MeV.cm,
											z.cm		=	z.cm,
											dz.cm		=	dz.cm,
											dE.MeV		=	dE.MeV,
											E.dep.MeV	=	E.dep.MeV,
											D.rel		=	dE.MeV / dz.cm / (df$dE.MeV[1] / df$dz.cm[1]),
											w.Edmund	=	dz.cm / z.total.cm)
		df				<-	rbind.data.frame(df, df.add)

	}

	return(df)
	
}


####### Energy deposition by delta-rays

E.delta <- function(E,R.e=delta.max,a0=1e-8,point.integrand=point.dose.int,
						dens=3.97,err=1e-16,N.el=1.2e24){
	# This function calculates the integral energy deposition
	# from delta-rays by the the Katz formula from a0 to R.max in
	# an infinitisimal thin slap of the detector. 
	# The output is given in J/cm.
	# The function is meant as input for the Hansen model.
	
	# E			= proton energy				 	[MeV]
	# R.max	= maximum electorn range	 	[cm]
	# a0		= core radius (start point) 	[cm]
	# dens		= density of material			[g/cm^3]
	
	# Created 10.09.2006 by JE
	
	R.max <- R.e(E)
	
	# Total energy deposited in track
	
	E.track <- 1e-3*dens*2*pi*integrate(point.integrand,lower=a0,upper=R.max,rel.tol=err,
				 N.el=N.el,E.p=E,dens=dens)$integral # J/cm	

}

