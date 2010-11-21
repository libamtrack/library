################################################################################################
# Rd wrapper generator
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


# Clean workspace
rm(list = ls())

# Load the function information extracted by 'read.cur.description.R'
load("functions.ssd")

# Replacement for "grepl" function to ensure compatibility with R <= 2.9.0
grep.bool	<-	function(pattern, x, ...){
	results	<-	grep(pattern, x, ...)
	indices	<-	1:length(x)
	bool.res	<-	is.element(indices, results)
      return(bool.res)
}

# List of hard-coded Rd files for variable descriptions
hardcoded.variable.descriptions    <- c("E.MeV.u", "particle.no", "material.no", "er.model", "rdd.model", "fluence.cm2.or.dose.Gy", "number.of.field.components")

# Read in hard-coded examples
hardcoded.examples                 <- scan( "./package/man/examples.txt", 
                                            what  = "character",
                                            sep   = "\n") 

# Find first line of examples and read the names
pos.start.examples                 <- grep("::", hardcoded.examples)
names.examples                     <- substring( hardcoded.examples[pos.start.examples],
                                                 3,
                                                 nchar(hardcoded.examples[pos.start.examples]) - 2)

# Read in links to sources
links.to.sources                   <- scan( "./package/man/links.to.sources.txt", 
                                            what  = "character",
                                            sep   = "\n") 

# Find first line of links and read the names
idx.links.to.sources               <- grep("::", links.to.sources) + 1
names.links                        <- substring( links.to.sources[idx.links.to.sources - 1],
                                                 3,
                                                 nchar(links.to.sources[idx.links.to.sources - 1]) - 2)

# Read in type conversion table for C to R
source("type.conversion.R")


for(i in 1:length(functions)){
	# i <- 4
	cur.function <- functions[[i]]

	# replace "_" by "."
	cur.function$name <- gsub("_", ".", cur.function$name)
	cur.function$parameter.comment$name <- gsub("_", ".", cur.function$parameter.comment$name)

	# remove "[]"
	cur.function$parameter.comment$name <- gsub("[]", "", cur.function$parameter.comment$name, fixed = T)

      # Remove derivable array.size arguments from in-list
	ii <- which(grep.bool(pattern = "in", x = cur.function$parameter.comment$type)) 
      if(length(ii) > 0){
           derivable.array.size.names <- gsub("_", ".", cur.function$parameter$name[which(cur.function$parameter$derivable.array.size.variable)])
           for (j in 1:length(ii)){
                #DEBUG: j <- 1
                if(cur.function$parameter.comment$name[j] %in% derivable.array.size.names){
                     ii[j]  <- 0
                }
           }
           ii    <- ii[ii!=0]
      }

	# paste all input parameters together
	parameter.list <- character(0)
	if(length(ii) > 0){ 
		if(length(ii) > 1){ 
			for(i in ii[1:(length(ii) - 1)]){
				parameter.list <- paste(parameter.list, cur.function$parameter.comment$name[i], ", ", sep = "") 
				if(which(i == ii)%%5 == 0) parameter.list <- paste(parameter.list, "\n\t", sep = "") 

			}
		}
		parameter.list <- paste(parameter.list, cur.function$parameter.comment$name[ii[length(ii)]], sep = "") 
	}

	write("% Automatically created Rd file\n", file = paste("./package/man/", cur.function$name, ".Rd", sep = ""))

	# Create Rd description
      cur.description   <- character(0)
	cur.description   <- paste("% TODO File path/", cur.function$name, ".Rd", sep = "")
	cur.description   <- c(cur.description, paste("\\name{", cur.function$name, "}", sep = ""))
	cur.description   <- c(cur.description, paste("\\alias{", cur.function$name, "}", sep = ""))
	cur.description   <- c(cur.description, paste("\\title{", cur.function$name, "}", sep = ""))
	cur.description   <- c(cur.description, paste("\\description{", cur.function$description, "}", sep = ""))
	cur.description   <- c(cur.description, paste("\\usage{", cur.function$name, "(", parameter.list, ")\n}", sep = ""))

	# Add arguments
	cur.description   <- c(cur.description, "\\arguments{")
	para              <- cur.function$parameter.comment
	para.in           <- which(grep.bool(pattern = "in", x = para$type)) 
      # TODO: produce common code for skipping derivable array sizes here and above
	if(length(para.in) > 0){
           for (j in 1:length(para.in)){
                #DEBUG: j <- 1
                if(cur.function$parameter.comment$name[j] %in% derivable.array.size.names){
                     para.in[j]  <- 0
                }
           }
           para.in <- para.in[para.in!=0]
      }

	if(length(para.in) > 0){
            for(j in para.in){
		      # DEBUG: j <- para.in[1]
                  line.to.add    <- NULL
                  # Check if domuented arguments are found in comment. If so, cross-reference to them
                  para.comment   <- gsub("_", ".", para$comment[j], fixed = T)
                  for (k in 1:length(para$name)){
                       # DEBUG: k <- 1
                       name.match      <- regexpr(para$name[k], para.comment)
                       if(name.match > 0 & (para$name[k] %in% hardcoded.variable.descriptions)){          # if hardcoded description exists for an argument mentioned in another argument's description
                            para.comment       <- paste( substring(para.comment, 1, name.match - 1),      # paste a link to it into the description
                                                         paste("\\code{\\link{", para$name[k], "}}", sep = ""),
                                                         substring(para.comment, name.match + attr(name.match, "match.length"), nchar(para.comment)),
                                                         sep = "")
                       }
                  }
                  if (para$name[j] %in% hardcoded.variable.descriptions){
                       cur.description         <- c(cur.description, paste("  \\item{", para$name[j], "}{", para.comment, " (see also \\code{\\link{", para$name[j], "}}).}", sep = ""))
                  }else{
                       cur.description         <- c(cur.description, paste("  \\item{", para$name[j], "}{", para.comment, ".}", sep = ""))
                  }
		}
	}
	cur.description <- c(cur.description, "}")

	# Values
	cur.description <- c(cur.description, "\\value{\n% TODO proper return definition of lists!!! ADD NUMBER_OF_FIELD_COMPONENT_DESCRIBTION AGAIN!!!)")

	para.out <- grep("out", para$type)
	if(length(para.out) > 0){
		for(i in para.out){
			cur.description <- c(cur.description, paste("  \\item{", para$name[i], "}{", gsub("_", ".", para$comment[i]), "}", sep = ""))
		}
	}

	if(cur.function$type != "void"){
		para.return <- grep("return", para$type)
		if(length(para.return) > 0){
			for(i in para.return){
				cur.description <- c(cur.description, paste("  \\item{", para$name[i], "}{", para$name[i], "}", sep = ""))
			}
		}
	}
	cur.description <- c(cur.description, "}")

    # Add link to source if exists
      if(cur.function$name %in% names.links){
           idx.link              <- idx.links.to.sources[grep(cur.function$name, names.links)]
           cur.description       <- c(cur.description, "\\seealso{\nView the C source code here:\n")
           cur.description       <- c(cur.description, paste("\\url{", links.to.sources[idx.link]), "}", sep = "")
           cur.description       <- c(cur.description, "}")
           
      }

	  # Add example(s) if exists
      if(cur.function$name %in% names.examples){
           idx.start.example         <- pos.start.examples[grep(cur.function$name, names.examples)] + 1
           # if last example do not try to find end index but set it to last index
           if(grep(cur.function$name, names.examples) == length(pos.start.examples)){
               idx.end.example            <- length(hardcoded.examples)
           }else{
               idx.end.example            <- pos.start.examples[grep(cur.function$name, names.examples) + 1] - 1
           }

           # extract example text and paste into Rd description
           cur.description       <- c(cur.description, "\\examples{")
           cur.description       <- c(cur.description, paste(hardcoded.examples[idx.start.example:idx.end.example]))
           cur.description       <- c(cur.description, "}")
      }
     
    # Write current description to file
	write(c(cur.description, "\n"), file = paste("./package/man/", cur.function$name, ".Rd", sep = ""), append = T)
}




