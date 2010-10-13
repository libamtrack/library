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

for(i in 1:length(functions)){
	# i <- 8
	tmp <- functions[[i]]

	# replace "_" by "."
	tmp$name <- gsub("_", ".", tmp$name)
	tmp$parameter.comment$name <- gsub("_", ".", tmp$parameter.comment$name)

	# remove "[]"
	tmp$parameter.comment$name <- gsub("[]", "", tmp$parameter.comment$name, fixed = T)

	# paste all input parameters together
	ii <- which(grepl("in", tmp$parameter.comment$type)) # SG commented out: & tmp$parameter.comment$name != "n")
	parameter.list <- character(0)
	if(length(ii) > 0){ 
		if(length(ii) > 1){ 
			for(i in ii[1:(length(ii) - 1)]){
				parameter.list <- paste(parameter.list, tmp$parameter.comment$name[i], ", ", sep = "") 
				if(which(i == ii)%%5 == 0) parameter.list <- paste(parameter.list, "\n\t", sep = "") 

			}
		}
		parameter.list <- paste(parameter.list, tmp$parameter.comment$name[ii[length(ii)]], sep = "") 
	}

	write("% Automatically created Rd file\n", file = paste("./package/man/", tmp$name, ".Rd", sep = ""))

	##############################
	# create header
	##############################
	header <- character(0)

	header <- paste("% TODO File path/", tmp$name, ".Rd", sep = "")

	header <- c(header, paste("\\name{", tmp$name, "}", sep = ""))

	header <- c(header, paste("\\alias{", tmp$name, "}", sep = ""))

	header <- c(header, paste("\\title{", tmp$name, "}", sep = ""))

	header <- c(header, paste("\\description{", tmp$description, "}", sep = ""))

	header <- c(header, paste("\\usage{", tmp$name, "(", parameter.list, ")\n}", sep = ""))

	# start arguments
	header <- c(header, "\\arguments{")

	para <- tmp$parameter.comment
	para.in <- which(grepl("in", para$type)) # SG commented out: & para$name != "n")
	if(length(para.in) > 0){
		for(i in para.in){
			header <- c(header, paste("  \\item{", para$name[i], "}{", para$comment[i], "}", sep = ""))
		}
	}
	# end arguments
	header <- c(header, "}")

	# start value
	header <- c(header, "\\value{\n% TODO proper return definition of lists!!!)")

	para.out <- grep("out", para$type)
	if(length(para.out) > 0){
		for(i in para.out){
			header <- c(header, paste("  \\item{", para$name[i], "}{", para$comment[i], "}", sep = ""))
		}
	}

	if(tmp$type != "void"){
		para.return <- grep("return", para$type)
		if(length(para.return) > 0){
			for(i in para.return){
				header <- c(header, paste("  \\item{", para$name[i], "}{", para$name[i], "}", sep = ""))
			}
		}
	}

	# end value
	header <- c(header, "}")

	write(c(header, "\n"), file = paste("./package/man/", tmp$name, ".Rd", sep = ""), append = T)
}

#tmp
#header



