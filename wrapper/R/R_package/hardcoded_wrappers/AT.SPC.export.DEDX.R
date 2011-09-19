#############################################################
# Export of stopping powers from libamtrack for use in TRiP98
#
# S. Greilich, Jul 2011 / Sep 2011
#############################################################
AT.SPC.export.DEDX <- function(stopping.power.source.no, file.name.DEDX = NULL, particle.names = NULL, energy.MeV.u = NULL, write = TRUE)
{
# for documentation only: this code will convert a text file in ICRU49/73 format directly
#	df                <- read.table("Water.txt", header = FALSE, sep = "\t")
#	particle.names    <- c("H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar")
#	names(df)         <- c("E.MeV.u", particle.names)
#
#	df2               <- reshape( df,
#								  direction       = "long",
#								  idvar           = "E.MeV.u",
#								  timevar         = "particle.name",
#								  times           = particle.names,
#								  varying         = list(2:ncol(df)),
#								  v.names         = "LET.MeV.cm2.g")
#								  
#	df2$particle.name <- factor( df2$particle.name,
#								 levels     = particle.names,
#								 ordered    = TRUE)
#
#	row.names(df2)    <- 1:nrow(df2)

	if(is.null(file.name.DEDX)){
		file.name.DEDX    <- "libamtrack.dedx"
	}
	# If no particles given, use ICRU49/73 - use A = 2*Z for simplicity
	if(is.null(particle.names)){
		particle.names    <- c(	"H",  "He", "Li", "Be", "B", 
								"C",  "N",  "O",  "F",  "Ne", 
								"Na", "Mg", "Al", "Si", "P", 
								"S",  "Cl", "Ar")
		particle.nos      <- c(	"1H",   "4He",  "6Li",  "8Be",  "10B", 
								"12C",  "14N",  "16O",  "18F",  "20Ne", 
								"22Na", "24Mg", "26Al", "28Si", "30P", 
								"32S",  "34Cl", "36Ar")
	}
	# if no energy grid given, use ICRU49/73
	if(is.null(energy.MeV.u)){
		energy.MeV.u  	  <- c( 2.5e-02, 3.0e-02, 4.0e-02, 5.0e-02, 6.0e-02, 7.0e-02, 8.0e-02, 9.0e-02, 1.0e-01, 1.5e-01, 2.0e-01, 2.5e-01, 3.0e-01, 4.0e-01,
								5.0e-01, 6.0e-01, 7.0e-01, 8.0e-01, 9.0e-01, 1.0e+00, 1.5e+00, 2.0e+00, 2.5e+00, 3.0e+00, 4.0e+00, 5.0e+00, 6.0e+00, 7.0e+00,
								8.0e+00, 9.0e+00, 1.0e+01, 1.5e+01, 2.0e+01, 2.5e+01, 3.0e+01, 4.0e+01, 5.0e+01, 6.0e+01, 7.0e+01, 8.0e+01, 9.0e+01, 1.0e+02,
								1.5e+02, 2.0e+02, 2.5e+02, 3.0e+02, 4.0e+02, 5.0e+02, 6.0e+02, 7.0e+02, 8.0e+02, 9.0e+02, 1.0e+03)
	}
	
	# Get stopping power data
	df <- expand.grid( energy.MeV.u      = energy.MeV.u,
	                   particle.names    = particle.names,
					   particle.no       = 0,
					   stopping.power.MeV.cm2.g = 0)

	df$particle.no   <- particle.nos[which(df$particle.names, particle.names)]
	
	
	
	output, <-, "!filetype, , , , dEdx
	!fileversion, , , , 19980515
	!filedate, , , , "

	output, <-, paste(output,, date(),, "\n",, sep, =, "")
	output, <-, paste(output,, "!material, , , , H2O
	!density, , , , 1
	#############################################################
	#Export, of, SHIELD-HIT, Bethe, stopping, powers, for, use, in, TRiP98
	#
	#, S., Greilich,, Jul, 2011
	#############################################################
	\n",, sep, =, "")


	for(cur.particle, in, unique(df2$particle.name)){
		#, cur.particle, <-, unique(df2$particle.name)[1]
		ii, <-, df2$particle.name, ==, cur.particle
		output, <-, paste(output,, 
						"!projectile, , , , ",, cur.particle,, "\n",, 
						"#, E/(MeV/u), dE/dx(MeVcm**2/g)\n",
						"!dedx\n",, 
						sep, =, "")

		for(i, in, 1:nrow(df2[ii,])){
			output, <-, paste(output,
							df2[ii,]$E.MeV.u[i],, ", ",
							df2[ii,]$LET.MeV.cm2.g[i],, "\n",
							sep, =, "")
		}
	}

	if(write, ==, TRUE){
		write(output,, file, =, file.name)
	}
	return(df)
}