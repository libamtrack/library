try(setwd("C:/Data/__Studium/01-R_test/wrapper"), silent = T)

rm(list = ls())

source("type.conversion.R")

load("functions.ssd")

for(i in 1:length(functions)){
	tmp <- functions[[i]]$comment

	# replace "_" by "."
	tmp$name <- gsub("_", ".", tmp$name)
	tmp$parameter$name <- gsub("_", ".", tmp$parameter$name)

	# remove "[]"
	tmp$parameter$name <- gsub("[]", "", tmp$parameter$name, fixed = T)

	# paste all input parameters together
	ii <- which(grepl("in", tmp$parameter$type) & tmp$parameter$name != "n")
	parameter.list <- character(0)
	if(length(ii) > 0){ 
		if(length(ii) > 1){ 
			for(i in ii[1:(length(ii) - 1)]){
				parameter.list <- paste(parameter.list, tmp$parameter$name[i], ", ", sep = "") 
				if(which(i == ii)%%5 == 0) parameter.list <- paste(parameter.list, "\n\t", sep = "") 

			}
		}
		parameter.list <- paste(parameter.list, tmp$parameter$name[ii[length(ii)]], sep = "") 
	}

	write("// Automatically created Rd file\n", file = paste(tmp$name, ".Rd", sep = ""))

	##############################
	# create header
	##############################
	header <- character(0)

	header <- paste("% TODO File path/", tmp$name, ".Rd", sep = "")

	header <- c(header, paste("\\name{LibAmTrack}", sep = ""))

	header <- c(header, paste("\\alias{", tmp$name, "}", sep = ""))

	header <- c(header, paste("\\title{", tmp$name, "}", sep = ""))

	header <- c(header, paste("\\description{", tmp$description, "}", sep = ""))

	header <- c(header, paste("\\usage{", tmp$name, "(", parameter.list, ")\n}", sep = ""))

	# start arguments
	header <- c(header, "\\arguments{")

	para <- tmp$parameter
	para.in <- which(grepl("in", para$type) & para$name != "n")
	if(length(para.in) > 0){
		for(i in para.in){
			header <- c(header, paste("  \\item{", para$name[i], "}{", para$comment[i], "}", sep = ""))
		}
	}
	# end arguments
	header <- c(header, "}")

	# start value
	header <- c(header, "% TODO proper return definition of lists!!!)")
	header <- c(header, "\\value{")

	para.out <- grep("out", para$type)
	if(length(para.out) > 0){
		for(i in para.out){
			header <- c(header, paste("  \\item{", para$name[i], "}{", para$comment[i], "}", sep = ""))
		}
	}
	# end value
	header <- c(header, "}")

	write(c(header, "\n"), file = paste(tmp$name, ".Rd", sep = ""), append = T)
}
tmp
header



