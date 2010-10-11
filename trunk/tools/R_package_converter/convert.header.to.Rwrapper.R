################################################################################################
# R wrapper generator
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

write("// Automatically created wrapper file\n", file = "./package/R/wrapper.R")


for(i in 1:length(functions)){
	tmp <- functions[[i]]

	##############################
	# create header
	##############################
	para <- tmp$parameter

	# create header out of the para[in] parameter without length n (works every where???)
	const <- grepl("in", tmp$parameter.comment$type) 
	pos.in <- which(const & tmp$parameter.comment$name != "n")
	# get the position of the parameters from the position in the comments
	pos.in <- match(tmp$parameter.comment$name[pos.in], para$name)
	# replace "_" with "."
	para$name <- gsub("_", ".", para$name, fixed = T)

	header <- character(length(pos.in))
	header[1] <- paste(gsub("_", ".", tmp$name), " <- function( ", 
					para$name[pos.in[1]], ",", sep = "")
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

	# assign length n of vectors
	if(sum(para$name == "n") > 0){
		ii <- which(para$length == "n")
		header <- c(header, paste("\tn\t<- length(", para$name[ii[1]], 
						")", sep = ""))
	}


	# assign the results (from para["out"])
	pos.out <- grep("out", tmp$parameter.comment$type)
	# get the position of the parameters from the position in the comments
	pos.out <- match(tmp$parameter.comment$name[pos.out], gsub(".", "_", para$name, fixed = T))

	# check if output parameter is also input
	check <- match(pos.out, pos.in)	
	pos.only.out <- pos.out[is.na(check)]
	
	if(length(pos.only.out) > 0){
		for(j in 1:(length(pos.only.out))){
			header <- c(header, paste( "\t", para$name[pos.only.out[j]], 
								" = numeric(", para$length[pos.only.out[j]],
								")", sep = ""))
		}
	}
	
	header <- c(header, "")

	# function call
	header <- c(header, paste("\tres <- .c(\"", tmp$name, "_R\", n = as.integer(n),", sep = ""))

	if(nrow(para) > 2)
	for(j in 2:(nrow(para) - 1)){
		header <- c(header, paste("\t\t\t", gsub("_", ".", para$name[j]), " = ", convert.R(para$type[j]), "(",  gsub("_", ".", para$name[j]), "),", sep = ""))
	}
	
	header <- c(header, paste("\t\t\t", gsub("_", ".", para$name[nrow(para)]), " = ", convert.R(para$type[j]), "(",  gsub("_", ".", para$name[nrow(para)]), "))", sep = ""))

	header <- c(header, "")

	# create list to return, if more than one parameter is returned
	if(length(pos.out) == 1){
		header <- c(header, paste("\t return(res$", para$name[pos.out[length(pos.out)]], ")", sep = ""))
		header <- c(header, "}")
	} else if (length(pos.out) > 1){
		header <- c(header, "\t return.list <- NULL")
		for(j in 1:length(pos.out)){
			header <- c(header, paste("\t return.list[[", j, "]] <- res$", para$name[pos.out[j]], "", sep = ""))
		}
		header <- c(header, "\t return(return.list)")
		header <- c(header, "}")
	}


	write(c(header, "\n"), file = "./package/R/wrapper.R", append = T)

}
#tmp
#header
