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

source("../../../tools/automatic_wrapper_generator/R.type.conversion.R")

load("functions.sdd")

write("#ifndef AT_R_WRAPPER_H_\n#define AT_R_WRAPPER_H_\n// Automatically created header file\n\n#include <stdlib.h>\n#include <stdbool.h>\n", file = "AT_R_Wrapper.h")
# Add include files for functions used
used.header.files <- unique(sapply(functions, function(x){x$header.file.name}))
for (header.file in used.header.files){
   write(paste("#include\"",
               header.file,
               "\"", 
               sep = ""),
         file = "AT_R_Wrapper.h", append = TRUE)  
}
write("\n", file = "AT_R_Wrapper.h", append = TRUE)  

write("// Automatically created header and body file\n\n#include \"AT_R_Wrapper.h\"\n", file = "AT_R_Wrapper.c")

# replacement for "grepl" function to ensure compatibilty with R <= 2.9.0
grep.bool	<-	function(pattern, x, ...){
	results	<-	grep(pattern, x, ...)
	indices	<-	1:length(x)
	bool.res	<-	is.element(indices, results)
}

for(i in 1:length(functions)){
	# i <- 51
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
	write(c(header, "\n"), file = "AT_R_Wrapper.h", append = T)

	###########################
	# create function body
	###########################

	# get parameter information for current functions
	para         <- functions[[i]]$parameter

	# function declaration
	body         <- gsub("\n", "", gsub(";", "", header), fixed = T)
	body         <- c(body, "{")

	# select input parameters and arrays (vectors)
      input        <- grep.bool(pattern = "in", x = para$in.out)
	vector       <- (para$length != 1)

	# create count variable i, if sum(vector) > 0
	if(sum(vector) > 0){
		body        <- c(body, "  long i;")
	}

	# add the input parameters	
	ii           <- which(input & !vector)
	if(length(ii) > 0){
		for(j in ii){
			body <- c(body, paste("  ", gsub("*", "", para$type[j], , fixed = T), " ", para$name[j], get.extension(para$type[j]), " = (",
					gsub("*", "", gsub("const ", "", para$type[j]), fixed = T), ")(*",  para$name[j], 
					");", sep = ""))
		}
	}
	body         <- c(body, "")

	# check for those array lengths that are given by a variable 
	numbers               <- para$length%in%para$name
  para.length           <- para$length
  
	# add "_long" to all variables that are no numbers
	para.length[numbers]  <- paste(para.length[numbers], "_long", sep = "")

	# write changed para.length back into data.frame
	para$length           <- para.length
	para.length           <- unique(para$length)[unique(para$length) != 1]

	for(l in para.length){
		
		jj <- which(input & vector & para$length == l)
		if(length(jj) != 0){
			for(k in jj){
##
	            if(grepl("char", para$type[k]) != TRUE){
##
					body <- c(body, "\n//Allocate space for the input parameter.")
					body <- c(body, paste("  ", gsub("const ", "", para$type[k]), "* ",
							para$name[k],  get.extension(para$type[k]),
							" = (", gsub("const ", "", para$type[k]), 
							"*)calloc(", l,  ",sizeof(", gsub("const ", "", para$type[k]),
							"));", sep = ""))
##
				}
##
			}
			body <- c(body, "")
			# fill the data into the allocated space
			if(sum(grepl("char", para$type[jj]))==0){
        body <- c(body, "\n//Fill in the input parameter.")
			  body <- c(body, paste("  for(i = 0 ; i < ", l, "; i++){", sep = ""))
			  for(k in jj){
				  body <- c(body, paste("\t", para$name[k], get.extension(para$type[k]),
						  "[i] = (", gsub("const ", "", para$type[k]), ")", para$name[k],
						  "[i];", sep = ""))
			  }
			  body <- c(body, "  }")
		  }else{
#        body <- c(body, "\n// Copy strings\n")
#  		  for(k in jj){
#				  body <- c(body, paste("\tstrcpy(", 
#                                para$name[k], get.extension(para$type[k]),
#                                ",(*",
#                                para$name[k], 
#                                "));",
#                                sep = ""))
#			  }
		  }
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
      
      # Add definition of non-array pointers that are used as output variables
      idx.output.pointers     <- grep.bool(pattern = "out", x = para$in.out) & para$length == 1 & !grep.bool(pattern = "in.out", x = para$in.out)
      if(sum(idx.output.pointers) > 0){
      output.pointers         <- para[idx.output.pointers, ]
      body                    <- c(body, "\n//Define type-casted output variables")
      for (i in 1:nrow(output.pointers)){
	     # DEBUG: i <- 1
           output.pointers$type[i] <- gsub("*", "", output.pointers$type[i], fixed = T)
           body                    <- c( body, 
                                         paste( "\t", 
                                                output.pointers$type[i], 
                                                " ", 
                                                output.pointers$name[i], 
                                                get.extension(output.pointers$type[i]),
						            " = 0;", 
                                                sep = ""))
      }
	}

	########################
      # Generate function call
	########################
	if(tmp$type != "void"){
		return.var.txt	<-	paste(tmp$type, "returnValue_internal = \t")
	}else{
		return.var.txt	<-	""
	}

   
	if(grepl("char", para$type[1], fixed = TRUE) == TRUE & grepl("char*", para$type[1], fixed = TRUE) == FALSE){
		body <- c(body, paste("\n  ", return.var.txt, tmp$name, "( ", para$name[1], 
					 get.extension(para$type[1]), "[0],", sep = ""))
	}else{
		body <- c(body, paste("\n  ", return.var.txt, tmp$name, "( ", para$name[1], 
					 get.extension(para$type[1]), ",", sep = ""))
	}
	if(tmp$type != "void"){
		para.max	<-	nrow(para) - 1
	}else{
		para.max	<-	nrow(para)
	}

	for(j in 2:(para.max - 1)){
		if((length(grep("*", para$type[j], fixed = TRUE)) == 0) | (grepl("char*", para$type[j], fixed = TRUE) == TRUE)){
				if(grepl("char", para$type[j], fixed = TRUE) == TRUE & grepl("char*", para$type[j], fixed = TRUE) == FALSE){
				    body <- c(body, paste("\t", para$name[j], get.extension(para$type[j]), "[0],", sep = ""))
				}else{
				    body <- c(body, paste("\t", para$name[j], get.extension(para$type[j]), ",", sep = ""))
				}
			}
		else
			body <- c(body, paste("\t&", para$name[j], get.extension(para$type[j]), ",", sep = ""))
	} 			
	
	if((length(grep("*", para$type[para.max], fixed = TRUE)) == 0) | (grepl("char*", para$type[para.max], fixed = TRUE) == TRUE)){
	  	if(grepl("char", para$type[para.max], fixed = TRUE) == TRUE & grepl("char*", para$type[para.max], fixed = TRUE) == FALSE){
	    	body <- c(body, paste("\t", para$name[para.max], get.extension(para$type[para.max]), "[0]);", sep = ""))
	    }else{
	    	body <- c(body, paste("\t", para$name[para.max], get.extension(para$type[para.max]), ");", sep = ""))
	    }
      }else{
	      body <- c(body, paste("\t&", para$name[para.max], get.extension(para$type[para.max]), ");", sep = ""))
      }

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
#			kk <-  which(!input & vector & para$length == l)
#			kk <-  which(vector & para$length == l)
		kk <-  which(output & vector & para$length == l)
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
	write(body, file = "AT_R_Wrapper.c", append = T)
}

write("#endif\n", file = "AT_R_Wrapper.h", append = T)