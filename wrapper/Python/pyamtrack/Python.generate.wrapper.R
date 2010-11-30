################################################################################################
# Python wrapper generator
################################################################################################
# Copyright 2006, 2010 The libamtrack team
# 
# This file is part of the AmTrack program (libamtrack.sourceforge.net).
#
#    Created on: 29.11.2010
#    Creator: Steffen Greilich
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

write("# Automatically created wrapper file\n", file = "./libamtrack.R")

# replacement for "grepl" function to ensure compatibilty with R <= 2.9.0
grep.bool	<-	function(pattern, x, ...){
	results	<-	grep(pattern, x, ...)
	indices	<-	1:length(x)
	bool.res	<-	is.element(indices, results)
}

for(i in 1:length(functions)){
	#i<-3
	tmp <- functions[[i]]

	##############################
	# create header
	##############################
	para          <- tmp$parameter

	# Create header with the para[in] parameter, skip array sizes that can be derived from
      # given arrays

	const                    <- grep.bool(pattern = "in", x = tmp$parameter.comment$type) 
	pos.in                   <- which(const) 
      pos.in                   <- pos.in[!(pos.in %in% para$array.size.defined.by.variable.no)]

	# get the position of the parameters from the position in the comments
	pos.in        <- match(tmp$parameter.comment$name[pos.in], para$name)
	# replace "_" with "."
	para$name     <- gsub("_", ".", para$name, fixed = T)

	header <- character(length(pos.in))
	if(length(header) == 1){
           header[1] <- paste(gsub("_", ".", tmp$name), " <- function( ", 
                                   para$name[pos.in[1]], "){\n", sep = "")
	}else{
           header[1] <- paste(gsub("_", ".", tmp$name), " <- function( ", 
                                   para$name[pos.in[1]], ",", sep = "")
      }

      if(length(header) > 2){
		for(j in 2:(length(header)-1)){
			header[j] <- paste( "\t\t\t", para$name[pos.in[j]], ",", sep = "")
		}
	}
	if(length(header) > 1){
		header[length(header)] <- paste( "\t\t\t", para$name[pos.in[length(pos.in)]], 
								"){\n", sep = "")
	}

	###########################
	# create wrapper body
	###########################
      
  	# Derive array sizes
      if(sum(para$derivable.array.size.variable) > 0){
           for(j in 1:length(para$derivable.array.size.variable)){
                # DEBUG: j <- 1
                if(para$derivable.array.size.variable[j]){
                     idx     <- grep(j, para$array.size.defined.by.variable.no)
                     header  <- c( header, 
                                   paste( "\t",
                                          para$name[j], 
                                          "\t<- length(", 
                                          para$name[idx[1]], 
						                  ")", 
                                          sep = ""))
                     if(length(idx) > 1){
                          for(k in 2:length(idx)){
                               # DEBUG: k <- 2
                               header  <- c( header, 
                                             paste( "\tif(",
                                                    para$name[j], 
                                                    " != length(", 
                                                    para$name[idx[k]], 
		          		                      ")){cat(\"",
                                                    paste(  "Array size mismatch for \'", 
                                                            para$name[j],
                                                            "\'!\\n",
                                                            sep = ""),
                                                    "\")\n\t\treturn}\n",
                                                    sep = ""))
                          }
                     }
                }
           }
	}

      # Add missing zeros to arrays that should have a defined number of elements
      idx.fixed.arrays    <- which(!para$derivable.array.size.variable & para$fixed.size.array)
      if(length(idx.fixed.arrays) > 0){
           # DEBUG: j <- 1
           for(j in 1:length(idx.fixed.arrays)){
                header  <- c( header, 
                              paste( "\t",
                                     para$name[idx.fixed.arrays[j]],
                                     " <- c(",
                                     para$name[idx.fixed.arrays[j]],
                                     ", rep(0, ",
                                     para$length[idx.fixed.arrays[j]],
                                     " - length(",
                                     para$name[idx.fixed.arrays[j]],
                                     ")))\n",
                                     sep = ""))                                     
           }
      }

	# assign the results (from para["out"])
	pos.out <- grep("out", tmp$parameter.comment$type)
	# get the position of the parameters from the position in the comments
	pos.out <- match(tmp$parameter.comment$name[pos.out], gsub(".", "_", para$name, fixed = T))

	# check if output parameter is also input
	check <- match(pos.out, pos.in)	
	pos.only.out <- pos.out[is.na(check)]
	
	if(length(pos.only.out) > 0){
		# j <- 6
		for(j in 1:(length(pos.only.out))){
			header <- c(header, paste( "\t", para$name[pos.only.out[j]], 
								" <- numeric(", gsub("_", ".", para$length[pos.only.out[j]], fixed = T),
								")", sep = ""))
		}
	}
	
	# add return variable if existing
	if(tmp$type != "void"){
		header <- c(header, paste( "\t", "returnValue = numeric(1)", sep = ""))
	}
	header <- c(header, "")

	# function call
	header <- c(header, paste("\tres <- .C(\"", tmp$name, "_R\", ", sep = ""))

	if(nrow(para) > 1)
	for(j in 1:(nrow(para) - 1)){
		header <- c(header, paste("\t\t\t", gsub("_", ".", para$name[j]), " = ", convert.R(para$type[j]), "(",  gsub("_", ".", para$name[j]), "),", sep = ""))
	}
	
	header <- c(header, paste("\t\t\t", gsub("_", ".", para$name[nrow(para)]), " = ", convert.R(para$type[nrow(para)]), "(",  gsub("_", ".", para$name[nrow(para)]), "),PACKAGE=\"libamtrack\")", sep = ""))

	header <- c(header, "")

	# SG: always return list
	if (tmp$type != "void"){
		n.return.elements	<- 1
	}else{
		n.return.elements	<- 0
	}
	n.return.elements	<-	 n.return.elements + length(pos.out)
	if(n.return.elements>= 1){
		header <- c(header, "\t return.list <- NULL")
		j <- 0
            names.in.return.list   <- "c("
            if(length(pos.out) > 0){
			for(j in 1:length(pos.out)){
			# j <- 1
                        header                <- c(header, paste("\t return.list[[", j, "]] <- res$", para$name[pos.out[j]], "", sep = ""))
                        names.in.return.list  <- paste(names.in.return.list, "\"", para$name[pos.out[j]], "\"", sep = "")
                        if (j < length(pos.out)){
                             names.in.return.list  <- paste(names.in.return.list, ",", sep = "")
                        }
			}
		}
		if(tmp$type != "void"){
			header                <- c(header, paste("\t return.list[[", j+1, "]] <- res$returnValue", sep = ""))
                  if(n.return.elements == 1){ # skip comma if returnValue is the only parameter in list
                       names.in.return.list  <- paste(names.in.return.list, "\"returnValue\"", sep = "")
                  }else{
                       names.in.return.list  <- paste(names.in.return.list, ", \"returnValue\"", sep = "")
                  }
		}

            names.in.return.list  <- paste(names.in.return.list, ")", sep = "")
		
            # Add names to return list            
            header   <- c(header, paste("\t names(return.list) <- ", names.in.return.list, sep = ""))
            # Add "return" statement
            header <- c(header, "\t return(return.list)")
		header <- c(header, "}")
	}


	write(c(header, "\n"), file = "./libamtrack.R", append = T)

}
#tmp
#header
