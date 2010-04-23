run.SPIFF		<-	function(	E.MeV.u,
					particle.no,
					fluence.cm2,
					material.no				=	1,				# Liquid water
					RDD.model				=	3,				# Geiß' RDD
					RDD.parameters			=	5e-8,			# a0 = 50 nm
					ER.model				=	4,				# Geiß' ER
					gamma.model				=	2,				# General hit-target
					gamma.parameters		=	c(1,10,1,1,0),	# One single-hit-single-target (exp-sat) component, characteristic dose 10 Gy
					N2						=	20,				# 20 bins per factor 2 in histograms
					fluence.factor			=	1,				# use fluence as given
					write.output			=	F,				# no log file
					shrink.tails			=	T,				# cut insignificant tails
					shrink.tails.under		=	1e-30,			# cut them in case contribution to first moment is lower than
					adjust.N2				=	T,				# adjust bin width during convolution
					lethal.events.mode		=	F,				# use survival instead of activiation
					verbose					=	F){				# return distributions etc.

if(verbose){
	results.1	<-	AT.SC.get.f1.array.size(E.MeV.u = E.MeV.u,
								particle.no = particle.no,
								fluence.cm2 = fluence.cm2,
								material.no = material.no,
								RDD.model = RDD.model,
								RDD.parameters = RDD.parameters,
								ER.model = ER.model,
								N2 = N2)

	results.2	<-	AT.SC.get.f1(	E.MeV.u = E.MeV.u,
							particle.no = particle.no,
							fluence.cm2 = fluence.cm2,
							material.no = material.no,
							RDD.model = RDD.model,
							RDD.parameters = RDD.parameters,
							ER.model = ER.model,
							N2 = N2,
							n.bins.f1 = results.1$n.bins.f1,
							f1.parameters = results.1$f1.parameters)

	results.3	<-	AT.SC.get.f.array.size(	u = results.2$f.parameters[1],
								fluence.factor = fluence.factor,
								N2 = N2,
								n.bins.f1 = results.1$n.bins.f1,
								f1.d.Gy = results.2$f1$f1.d.Gy,
								f1.dd.Gy = results.2$f1$f1.dd.Gy,
								f1 = results.2$f1$f1)

	results.4	<-	AT.SC.get.f.start(	N2 = N2,
								n.bins.f1 = results.1$n.bins.f1,
								f1.d.Gy = results.2$f1$f1.d.Gy,
								f1.dd.Gy = results.2$f1$f1.dd.Gy,
								f1 = results.2$f1$f1,
								n.bins.f = results.3$n.bins.f)

	results.5	<-	AT.SC.SuccessiveConvolutions(	u = results.2$f.parameters[1],
									n.bins.f = results.3$n.bins.f,
									N2 = N2,
									n.bins.f.used = results.1$n.bins.f1,
									f.d.Gy = results.4$f.start.d.Gy,
									f.dd.Gy = results.4$f.start.dd.Gy,
									f = results.4$f.start,
									write.output = write.output,
									shrink.tails = shrink.tails,
									shrink.tails.under = shrink.tails.under,
									adjust.N2 = adjust.N2)

	results.6	<-	AT.SC.get.gamma.response(	d.Gy = results.5$f$f.d.Gy,
									dd.Gy = results.5$f$f.dd.Gy,
									f = results.5$f$f,
									f0 = results.5$f0,
									gamma.model = gamma.model,
									gamma.parameters = gamma.parameters,
									lethal.events.mode = lethal.events.mode)
	
	index			<-	1:(length(results.1$f1.parameters) / 8) * 8 - 8
	df.f1.parameters	<-	data.frame(	E.MeV.u				=	E.MeV.u,
							particle.name			=	AT.particle.name.from.particle.no(particle.no),
							LET.MeV.cm2.g			=	results.1$f1.parameters[index + 1],
							r.min.m				=	results.1$f1.parameters[index + 2],
							r.max.m				=	results.1$f1.parameters[index + 3],
							d.min.Gy				=	results.1$f1.parameters[index + 4],
							d.max.Gy				=	results.1$f1.parameters[index + 5],
							normalization.constant.Gy 	=	results.1$f1.parameters[index + 6],
							single.impact.fluence.cm2	=	results.1$f1.parameters[index + 7],
							single.impact.dose.Gy		=	results.1$f1.parameters[index + 8])

	df.f.parameters	<-	data.frame(	u					=	results.2$f.parameters[1],
							total.fluence.cm2			=	results.2$f.parameters[2],
							total.dose.Gy			=	results.2$f.parameters[3],
							ave.E.MeV				=	results.2$f.parameters[4],
							dw.E.MeV				=	results.2$f.parameters[5],
							ave.LET.MeV.cm2.g			=	results.2$f.parameters[6],
							dw.LET.MeV.cm2.g			=	results.2$f.parameters[7])

	return(	list(		E.MeV.u = E.MeV.u,
					particle.no = particle.no,
					fluence.cm2 = fluence.cm2,
					norm.fluence = results.2$norm.fluence,
					dose.contribution.Gy = results.2$dose.contribution.Gy,
					material.no = material.no,
					RDD.model = RDD.model,
					RDD.parameters = RDD.parameters,
					ER.model = ER.model,
					N2.set = N2,
					N2 = results.5$N2,
					f1.parameters = df.f1.parameters,
					f.parameters = df.f.parameters,
					n.bins.f1 = results.2$n.bins.f1,
					n.bins.f = results.3$n.bins.f,
					n.bins.f.used = results.5$n.bins.f.used,
					u.start = results.3$u.start,
					n.convolutions = results.3$n.convolutions,
					fluence.factor = fluence.factor,
					write.output = write.output,
					shrink.tails = shrink.tails,
					shrink.tails.under = shrink.tails.under,
					adjust.N2 = adjust.N2,
					d = results.5$d,
					f1 = results.2$f1,
					f.start = results.4,
					f0 = results.5$f0,
					f = results.5$f))
	}else{
		return(AT.run.SPIFF(	E.MeV.u 		= E.MeV.u,
							particle.no 	= particle.no,
							fluence.cm2 	= fluence.cm2,
							material.no 	= material.no,
							RDD.model 		= RDD.model,
							RDD.parameters 	= RDD.parameters,
							ER.model 		= ER.model,
							gamma.model		= gamma.model,
							gamma.parameters	= gamma.parameters,
							N2 			= N2,
							fluence.factor 	= fluence.factor,
							write.output 	= write.output,
							shrink.tails	= shrink.tails,
							shrink.tails.under = shrink.tails.under,
							adjust.N2 		= adjust.N2,
							lethal.events.mode	= lethal.events.mode)$results)
	}
}
