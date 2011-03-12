# R function for reading of spc files 
#
# started 2010/02/18, sgre
# as R script
#
# revised 2010/05/18, fk
# added: show.bin/show.bin.width feature, set type of mean.
# Bug fixed: tmp.data$Cum changed to tmp.data$HCum after tmp.data$H.
#
# revised 2010/07/28, sgre
# removed DDD handling
#
# revised 2010/08/04, sgre
# now choice of endian
#
# revised 2010/11/02, sgre
# added to R package
#
# revised 2010/11/10, sgre
# switched to matrix ops instead of dataframes, optimized for libamtrack, removed raw option
#
# revised 2010/11/11, sgre
# speed-up due to working on index tables, prototype for C version
#
# revised 2011/03/12, sgre
# after being superseeded by C version (r781) which
# still at troubles at r904 rolled back and combined
# with C version
# Now both R ('vanilla') and C (unstable) version can be used

AT.SPC.read <- function( file.name, 
                         flavour        = "vanilla",
						 endian         = c("big", "little")[1], 
                         mean           = c("geometric", "arithmetic")[2],
                         compress       = TRUE)
{
	if(flavour == "vanilla"){
    # R version
		to.read                 <- file(file.name, "rb")

		# Scan file first to get bin sizes
		seek(to.read, where = 0, origin = 'end')  
		end.pos                 <- seek(to.read, where = NA)  
		seek(to.read, where = 0, origin = 'start')
		tags                    <- 1
		max.tags                <- 20000
		mtag                    <- matrix( data = 0, ncol = 4, nrow = max.tags)
		repeat{
		   mtag[tags, 1]  <- seek(to.read, where = NA)  
		   mtag[tags, 2]  <- readBin(to.read, integer(), endian = endian)
		   if(mtag[tags, 2] > 20){
			  cat("Strange data read. Probably wrong endianess given.\n")
			  close(to.read)
			  stop()
		   }
		   mtag[tags, 3]  <- readBin(to.read, integer(), endian = endian)
		   if( mtag[tags, 2] %in% c(9,12,16,18)){
			   mtag[tags, 4]  <- readBin(    to.read, integer(), size = 8, signed = FALSE, n = 1, endian = endian)
		   }else if(mtag[tags, 2] == 10){
			   mtag[tags, 4]  <- readBin(    to.read, double(), size = 8, n = 1, endian = endian)
		   }else if(mtag[tags, 2] == 13){
			   seek(to.read, where = 16, origin = 'current')
			   mtag[tags, 4]  <- 1000 * readBin( to.read, integer(), size = 4, signed = TRUE, n = 1, endian = endian) + readBin( to.read, integer(), size = 4, signed = TRUE, n = 1, endian = endian)
		   } else{
			   seek(to.read, where = mtag[tags, 3], origin = 'current')
		   }
		   tags <- tags + 1
		   if(tags >= max.tags){
			  cat("Exceeded maximum number of tags. Please choose higher number.\n")
			  close(to.read)
			  stop()
		   }
		   if(seek(to.read, where = NA) >= end.pos){
			   break()
		   }
		}

		mtag               <- mtag[mtag[,1] != 0,]
		n.depth.steps      <- mtag[mtag[,2] == 9, 4]
		depth.g.cm2        <- mtag[mtag[,2] == 10, 4]
		n.particle.species <- mtag[mtag[,2] == 12, 4]
		ref.species        <- mtag[mtag[,2] == 18, 4]
		particle.no        <- mtag[mtag[,2] == 13, 4]
		bins               <- mtag[mtag[,2] == 16, 4]

		seek(to.read, where = mtag[mtag[,2] == 5,1] + 8, origin = 'start')
		projectile         <- rawToChar(readBin(	to.read, raw(), n = mtag[mtag[,2] == 5,3], endian = endian))	
		seek(to.read, where = mtag[mtag[,2] == 4,1] + 8, origin = 'start')
		target.material    <- rawToChar(readBin(	to.read, raw(), n = mtag[mtag[,2] == 4,3], endian = endian))	
		seek(to.read, where = mtag[mtag[,2] == 6,1] + 8, origin = 'start')
		beam.energy.MeV.u  <- readBin( to.read, double(), size = 8, n = floor(mtag[mtag[,2] == 6,3]/8), endian = endian)
		seek(to.read, where = mtag[mtag[,2] == 7,1] + 8, origin = 'start')
		peak.position.g.cm2<- readBin( to.read, double(), size = 8, n = floor(mtag[mtag[,2] == 7,3]/8), endian = endian)
		
		mm                 <- matrix(data = 0, ncol = 7, nrow = sum(bins))
		mm[,1]             <- rep.int(rep.int(1:n.depth.steps, n.particle.species), bins)
		mm[,2]             <- rep.int(rep.int(depth.g.cm2, n.particle.species), bins)
		mm[,3]             <- rep.int(sequence(n.particle.species), bins)
		mm[,4]             <- rep.int(particle.no, bins)

		fluence.tags       <- mtag[mtag[,2] == 19,]
		idx                <- 1
		for(i in 1:nrow(fluence.tags)){
			# i <- 1
			seek(to.read, where = fluence.tags[i,1] + 8, origin = 'start')
			size                   <- floor(fluence.tags[i,3]/8)
			mm[idx:(idx+size-1),7] <- readBin( to.read, double(), size = 8, n = size, endian = endian)
			idx                    <- idx + size
		}

		E.grid.tags        <- cbind( mtag[mtag[,2] %in% c(17,18),], 
									 rep.int(1:n.depth.steps, n.particle.species),   # depth step
									 sequence(n.particle.species),                   # species
									 cumsum(bins)-bins[1]+1,                         # index in mm
									 bins)                                           # size
		idx                <- 1
		for(i in 1:nrow(E.grid.tags)){
			# i <- 2
			if(E.grid.tags[i,2] == 17){
				seek(to.read, where = E.grid.tags[i,1] + 8, origin = 'start')
				size                   <- floor(E.grid.tags[i,3]/8)
				E.bins.MeV.u           <- readBin( to.read, double(), size = 8, n = size, endian = endian)
				E.low.MeV.u            <- E.bins.MeV.u[-length(E.bins.MeV.u)]
				E.high.MeV.u           <- E.bins.MeV.u[-1]
				mm[idx:(idx+size-2),6] <- E.high.MeV.u - E.low.MeV.u
				if(mean == "geometric"){  
					mm[idx:(idx+size-2),5]  <- sqrt(E.low.MeV.u * E.high.MeV.u)
				}else{
					mm[idx:(idx+size-2),5]  <- (E.low.MeV.u + E.high.MeV.u)/2
				}            
				idx                    <- idx + size - 1
			}else{
				ii                      <- E.grid.tags[,5] == E.grid.tags[i,5]
				ref.idx                 <- which((E.grid.tags[i,4]+1) == E.grid.tags[ii,6])
				size                    <- E.grid.tags[ii,8][ref.idx]
				from                    <- E.grid.tags[ii,7][ref.idx]
				mm[idx:(idx+size-1),5]  <- mm[from:(from+size-1),5]
				mm[idx:(idx+size-1),6]  <- mm[from:(from+size-1),6]
				idx                     <- idx + size
			}
		}
		close(to.read)

		# Remove zero fluence bins if requested
		if(compress == TRUE){
		  ii         <- mm[,7] != 0
		  mm         <- mm[ii,]
		}
		mm           <- mm[,-3]
		mm[,6]       <- mm[,6] * mm[,5]         # convert fluence / binwidth -> fluence
		df           <- as.data.frame(mm)
		names(df)    <- c("depth.step", "depth.g.cm2", "particle.no", "E.MeV.u", "DE.MeV.u", "fluence.cm2")

		cat(paste("Read ", n.depth.steps, " depth steps for projectile ", projectile, " on ", target.material, " with ", beam.energy.MeV.u, " MeV/u and peak at ", peak.position.g.cm2, " g/cm2.\n", sep = ""))
		
	}else{
	# C version
		 spc.size            <- numeric(1)
		 res                 <- .C( "AT_SPC_get_size_from_filename_R",
									file.name          = as.character(file.name),
									spc.size           = as.integer(spc.size),
									PACKAGE            = "libamtrack")

		n                    <- res$spc.size
		depth.step           <- integer(n)
		depth.g.cm2          <- numeric(n)
		E.MeV.u              <- numeric(n)
		DE.MeV.u             <- numeric(n)
		particle.no          <- integer(n)
		fluence.cm2          <- numeric(n)
		n.bins.read          <- integer(1)
		
		res                  <- .C( "AT_SPC_read_data_from_filename_R",
									file.name          = as.character(file.name),
									n                  = as.integer(n),
									depth.step         = depth.step,
									depth.g.cm2        = depth.g.cm2,
									E.MeV.u            = E.MeV.u,
									DE.MeV.u           = DE.MeV.u,
									particle.no        = particle.no,
									fluence.cm2        = fluence.cm2,
									n.bins.read        = n.bins.read,
									PACKAGE            = "libamtrack")
		
		df   <- data.frame(  depth.step             = res$depth.step,
							 depth.g.cm2            = res$depth.g.cm2,
							 E.MeV.u                = res$E.MeV.u,
							 DE.MeV.u               = res$DE.MeV.u,
							 particle.no            = res$particle.no,
							 fluence.cm2            = res$fluence.cm2)

		if (compress == TRUE){
			df  <- df[df$fluence.cm2 != 0,]
		}	
	}
	return(df)
}

