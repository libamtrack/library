try(setwd("C:/Data/__Studium/01-R_test/wrapper"), silent = T)

rm(list = ls())

source("type.conversion.R")

load("functions.ssd")

write("// Automatically created header file\n", file = "header.h")

write("// Automatically created header and body file\n", file = "header.c")

for(i in 1:length(functions)){
	tmp <- functions[[i]]
	
	cat(tmp[[1]],"\n")

	##############################
	# create header
	##############################
	header <- character(nrow(tmp$parameter))
	
	header[1] <- paste(tmp$type, " ", tmp$name, "_R( ", 
					convert.c(tmp$parameter$type[1]), "\t",
					tmp$parameter$name[1], ",", sep = "")

	if( tmp[[1]] == "AT_SuccessiveConvolutions" )
		cat( "\t",convert.c(tmp$parameter$type[1])," (","\t",tmp$parameter$type[1],")","\t",tmp$parameter$name[1], "\n" )

	if(length(header) > 2){
		for(j in 2:(length(header)-1)){
			header[j] <- paste( "\t", convert.c(tmp$parameter$type[j]), " ",
							tmp$parameter$name[j], ",", sep = "")
							
			if( tmp[[1]] == "AT_SuccessiveConvolutions" )
				cat( "\t",convert.c(tmp$parameter$type[j])," (","\t",tmp$parameter$type[j],")","\t",tmp$parameter$name[j], "\n" )
		}
	}

	if(length(header) > 1){
		header[length(header)] <- paste( "\t", convert.c(tmp$parameter$type[length(header)]), " ",
							tmp$parameter$name[length(header)], ");", sep = "")

		if( tmp[[1]] == "AT_SuccessiveConvolutions" )
			cat( "\t",convert.c(tmp$parameter$type[length(header)])," (","\t",tmp$parameter$type[length(header)],")","\t",tmp$parameter$name[length(header)], "\n" )

	}
	write(c(header, "\n"), file = "header.h", append = T)

	###########################
	# create function body
	###########################

	body <- gsub("\n", "", gsub(";", "", header), fixed = T)
	
	# open body
	body <- c(body, "{")

	# conversion and allocation of input parameters
	para <- functions[[i]]$parameter

	input <- grepl("const", para$type, fixed = T)
	vector <- (para$length != 1)

	# create count variable i, if sum(vector) > 0
	if(sum(vector) > 0){
		body <- c(body, "  long i;")
	}

	# add the input parameters
	if( tmp[[1]] == "AT_SuccessiveConvolutions" ){
		ii <- which(input)
		if(length(ii) > 0){
			for(j in ii){
					cat( "\tpara.name=_", para$name[j], "_\ttype=_", para$type[j] ,"_\n", sep="")
			}
		}
	}
	
	
	ii <- which(input & !vector)
	if(length(ii) > 0){
		for(j in ii){
			body <- c(body, paste("  ", para$type[j], " ", para$name[j], "_", 
					gsub("const ", "", para$type[j]), " = (",
					gsub("const ", "", para$type[j]), ")(*",  para$name[j], 
					");", sep = ""))
			
#			if( tmp[[1]] == "AT_SuccessiveConvolutions" )
#				cat( "\tpara.name=_", para$name[j], "_\ttype=_", para$type[j] ,"_\n", sep="")
		}
	}
	body <- c(body, "")

	para.length <- unique(para$length)[unique(para$length) != 1]

	for(l in para.length){
		
		jj <- which(input & vector & para$length == l)
		if(length(jj) != 0){
			body <- c(body, "\n//Allocate space for the input parameter.")
			for(k in jj){
				body <- c(body, paste("  ", gsub("const ", "", para$type[k]), "* ",
						para$name[k], "_", gsub("const ", "", para$type[k]),
						" = (", gsub("const ", "", para$type[k]), 
						"*)calloc(", l, ",sizeof(", gsub("const ", "", para$type[k]),
						"));", sep = ""))
			}
			body <- c(body, "")
			# fill the data into the allocated space
			body <- c(body, "\n//Fill in the input parameter.")
			body <- c(body, paste("  for(i = 0 ; i < ", l, "; i++){", sep = ""))
			for(k in jj){
				body <- c(body, paste("\t", para$name[k], "_", gsub("const ", "", para$type[k]),
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
						para$name[j], "_", gsub("const ", "", para$type[j]),
						" = (", gsub("const ", "", para$type[j]), 
						"*)calloc(", l, ",sizeof(", gsub("const ", "", para$type[j]),
						"));", sep = ""))
			}
		}

	}
	
	# function call
	body <- c(body, paste("\n  ", tmp$name, "( ", para$name[1], 
				"_", gsub("const ", "", para$type[1]), ",", sep = ""))
	for(j in 2:(nrow(para) - 1)){
		body <- c(body, paste("\t", para$name[j], "_", 
				gsub("const ", "", para$type[j]), ",", sep = ""))
	} 			
	body <- c(body, paste("\t", para$name[nrow(para)], "_",
			gsub("const ", "", para$type[nrow(para)]), ");", sep = ""))


	body <- c(body, paste("\n//Results:"))
	kk <-  which(!input & !vector)
	if(length(kk) > 0){
		for(j in kk){
			body <- c(body, paste("  *", para$name[j],	" = (", 
					gsub("*\t\t", "", convert.c(para$type[j]), fixed = T),			
					 ")", para$name[j], "_", 
					gsub("const ", "", para$type[j]), 
					";", sep = ""))
			body <- c(body, "")
		}
	}

	ll <-  which(!input & vector)
	if(length(ll) > 0){
		for(l in para.length){
			kk <-  which(!input & vector & para$length == l)
			if(length(kk) > 0){
				# fill the data into the allocated space
				body <- c(body, paste("  for(i = 0 ; i < ", l, "; i++){", sep = ""))
				for(j in kk){
					body <- c(body, paste("\t", para$name[j],	"[i] = (", 
							gsub("*\t\t", "", convert.c(para$type[j]), fixed = T),			
							 ")", para$name[j], "_", 
							gsub("const ", "", para$type[j]), 
							"[i];", sep = ""))
				}
			}
		}
	}
	body <- c(body, "  }")

	body <- c(body, "\n//Free allocated space")
	kk <- which(vector)	
	if(length(kk) > 0){
		for(j in kk){
			body <- c(body, paste("  free(", para$name[j],  "_", 
					gsub("const ", "", para$type[j]), ");", sep = ""))			
		}
	}

	# close body
	body <- c(body, "}\n\n")
	write(body, file = "header.c", append = T)
}

