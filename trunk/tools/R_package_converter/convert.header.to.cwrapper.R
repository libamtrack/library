################################################################################################
# C wrapper generator
################################################################################################
# Copyright 2006, 2010 The libamtrack team
# 
# This file is part of the AmTrack program (libamtrack.sourceforge.net).
#
#    Created on: 13.09.2010
#    Creator: kleinf
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

rm(list = ls())

source("type.conversion.R")

load("functions.ssd")

write("// Automatically created header file\n\n#include \"AmTrack.h\"\n#include <stdlib.h>\n#include <stdbool.h>\n", file = "./package/src/Rwrapper.h")

write("// Automatically created header and body file\n\n#include \"Rwrapper.h\"\n", file = "./package/src/Rwrapper.c")

# replacement for "grepl" function to ensure compatibilty with R <= 2.9.0
grep.bool	<-	function(pattern, x, ...){
	results	<-	grep(pattern, x, ...)
	indices	<-	1:length(x)
	bool.res	<-	is.element(indices, results)
}

for(i in 1:length(functions)){
	# i <- 7
	tmp <- functions[[i]]
	
	##############################
	# create header
	##############################
	header <- character(nrow(tmp$parameter))
	
#SG:	header[1] <- paste(tmp$type, " ", tmp$name, "_R( ", 
	header[1] <- paste("void ", tmp$name, "_R( ", 
					convert.c(tmp$parameter$type[1]), "\t",
					tmp$parameter$name[1], ",", sep = "")

	if(length(header) > 2){
		for(j in 2:(length(header)-1)){
			header[j] <- paste( "\t", convert.c(tmp$parameter$type[j]), " ",
							tmp$parameter$name[j], ",", sep = "")
		}
	}

	if(length(header) > 1){
		header[length(header)] <- paste( "\t", convert.c(tmp$parameter$type[length(header)]), " ",
							tmp$parameter$name[length(header)], ");", sep = "")
	}
	write(c(header, "\n"), file = "./package/src/Rwrapper.h", append = T)

	###########################
	# create function body
	###########################

	body <- gsub("\n", "", gsub(";", "", header), fixed = T)
	
	# open body
	body <- c(body, "{")

	# conversion and allocation of input parameters
	para <- functions[[i]]$parameter

	input <- grep.bool(pattern = "in", x = para$in.out)
	vector <- (para$length != 1)

	# create count variable i, if sum(vector) > 0
	if(sum(vector) > 0){
		body <- c(body, "  long i;")
	}

	# add the input parameters	
	ii <- which(input & !vector)
	if(length(ii) > 0){
		for(j in ii){
			body <- c(body, paste("  ", gsub("*", "", para$type[j], , fixed = T), " ", para$name[j], get.extension(para$type[j]), " = (",
					gsub("*", "", gsub("const ", "", para$type[j]), fixed = T), ")(*",  para$name[j], 
					");", sep = ""))
		}
	}
	body <- c(body, "")


	# get the length of the arrays and check for numbers (this will produce a warning!
	para.length <- para$length
	numbers <- !is.na(as.numeric(para.length)) 

	# add "_long" to all variables that are no numbers
	para.length[!numbers] <- paste(para.length[!numbers], "_long", sep = "")

	# write changed para.length back into data.frame
	para$length <- para.length

	para.length <- unique(para$length)[unique(para$length) != 1]

	for(l in para.length){
		
		jj <- which(input & vector & para$length == l)
		if(length(jj) != 0){
			body <- c(body, "\n//Allocate space for the input parameter.")
			for(k in jj){
				body <- c(body, paste("  ", gsub("const ", "", para$type[k]), "* ",
						para$name[k],  get.extension(para$type[k]),
						" = (", gsub("const ", "", para$type[k]), 
						"*)calloc(", l,  ",sizeof(", gsub("const ", "", para$type[k]),
						"));", sep = ""))
			}
			body <- c(body, "")
			# fill the data into the allocated space
			body <- c(body, "\n//Fill in the input parameter.")
			body <- c(body, paste("  for(i = 0 ; i < ", l, "; i++){", sep = ""))
			for(k in jj){
				body <- c(body, paste("\t", para$name[k], get.extension(para$type[k]),
						"[i] = (", gsub("const ", "", para$type[k]), ")", para$name[k],
						"[i];", sep = ""))
			}
			body <- c(body, "  }")
		}
		kk <-  which(!input & vector & para$length == l)
		if(length(kk) > 0){
			body <- c(body, "\n//Allocate space for the results.")
			for(j in kk){
				body <- c(body, paste("  ", gsub("const ", "", para$type[j]), "* ",
						para$name[j], get.extension(para$type[j]),
						" = (", gsub("const ", "", para$type[j]), 
						"*)calloc(", l, ",sizeof(", gsub("const ", "", para$type[j]),
						"));", sep = ""))
			}
		}

	}
	
	# function call
	if(tmp$type != "void"){
		return.var.txt	<-	paste(tmp$type, "returnValue_internal = \t")
	}else{
		return.var.txt	<-	""
	}

	body <- c(body, paste("\n  ", return.var.txt, tmp$name, "( ", para$name[1], 
				 get.extension(para$type[1]), ",", sep = ""))
	if(tmp$type != "void"){
		para.max	<-	nrow(para) - 1
	}else{
		para.max	<-	nrow(para)
	}

	for(j in 2:(para.max - 1)){
		if(length(grep("*", para$type[j], fixed = T)) == 0)
			body <- c(body, paste("\t", para$name[j], get.extension(para$type[j]), ",", sep = ""))
		else
			body <- c(body, paste("\t&", para$name[j], get.extension(para$type[j]), ",", sep = ""))
	} 			
	body <- c(body, paste("\t", para$name[para.max], get.extension(para$type[para.max]), ");", sep = ""))


	output <- grep.bool(pattern = "out", x = para$in.out) 

	body <- c(body, paste("\n//Results:"))
	if(tmp$type != "void"){
		body <- c(body, paste("\n\t*returnValue = (", gsub("*", "", gsub("\t", "", convert.c(para$type[nrow(para)]), fixed = T), fixed = T), ")returnValue_internal;"), sep = "")
	}
	
	kk <-  which(output & !vector)
	if(length(kk) > 0){
		for(j in kk){
			body <- c(body, paste("  *", para$name[j],	" = (", 
					gsub("\t", "", gsub("*", "", gsub("const ", "", convert.c(para$type[j])), fixed = T), fixed = T),			
					 ")", para$name[j], get.extension(para$type[j]), 
					";", sep = ""))
			body <- c(body, "")
		}
	}

	ll <-  which(output & vector)
	if(length(ll) > 0){
		for(l in para.length){
			kk <-  which(!input & vector & para$length == l)
			if(length(kk) > 0){
				# fill the data into the allocated space
				body <- c(body, paste("  for(i = 0 ; i < ", l, "; i++){", sep = ""))
				for(j in kk){
					body <- c(body, paste("\t", para$name[j],	"[i] = (", 
							gsub("\t", "", gsub("*", "", gsub("const ", "", convert.c(para$type[j])), fixed = T), fixed = T),			
							 ")", para$name[j], get.extension(para$type[j]),
							"[i];", sep = ""))
				}
				body <- c(body, "  }")
			}
		}
	}
	

	body <- c(body, "\n//Free allocated space")
	kk <- which(vector)	
	if(length(kk) > 0){
		for(j in kk){
			body <- c(body, paste("  free(", para$name[j], get.extension(para$type[j]), ");", sep = ""))			
		}
	}

	# close body
	body <- c(body, "}\n\n")
	write(body, file = "./package/src/Rwrapper.c", append = T)
}