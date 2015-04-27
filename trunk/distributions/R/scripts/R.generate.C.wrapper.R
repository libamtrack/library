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

# Clean workspace
rm(list = ls())

#Arguments passed by shell (arg1: path and arg2: library name)
args <- commandArgs(TRUE)

# Read type conversion functions (R <-> C types)
source(paste(args[1], "/R.type.conversion.R", sep = ""))

# Read function information from parsed doxygen comments
load("functions.sdd")

# Write first part of header wrapper file, macros to avoid multi include,
# include standard headers
write("#ifndef AT_R_WRAPPER_H_\n#define AT_R_WRAPPER_H_\n// Automatically created header file\n\n#include <stdlib.h>\n#include <stdbool.h>\n", 
      file = "AT_R_Wrapper.h")

# Add include statements for the functions used
used.header.files <- unique( sapply( functions, 
                                     function(x){x$header.file.name}))

for (header.file in used.header.files){
   write(paste("#include \"",
               header.file,
               "\"", 
               sep = ""),
         file = "AT_R_Wrapper.h", 
         append = TRUE)  
}

write("\n", file = "AT_R_Wrapper.h", append = TRUE)  


# Create source wrapper file
write("// Automatically created header and body file\n\n#include \"AT_R_Wrapper.h\"\n", 
      file = "AT_R_Wrapper.c")

##############################################
# Loop over all functions from doxygen parsing
for(i in 1:length(functions)){
  # i <- 1
	
  # Extract current function, parameters and number of parameters
  curFun    <- functions[[i]]
  para      <- functions[[i]]$parameter
  n.para    <- nrow(para)          

  cat( "Processing function no. ",
       i,
       " (",
       curFun$name,
       ")\n",
       sep = "")
  
  ################################
  # A. create function declaration
  ################################
  header    <- character(n.para + 1)
  # function name
  header[1] <- paste( "void ",         
                      curFun$name, 
                      "_R( ", 
	              sep = "")
  
  # function parameters
  if(n.para>0){            
    for(j in 1:n.para){
    # j <- 1
      
      # first line
      start.line <- "\t\t"
      if(j == 1){
        start.line <- ""
      }  

      # last line
      end.line    <- ""
      if(j != n.para){
        end.line    <- ","
      }

      # add j-th argument
      header[j] <- paste(   header[j],
                            start.line,
                            convert.c( curFun$parameter$type[j] ), 
                            "\t",
  			    curFun$parameter$name[j], 
                            end.line, 
                            sep = "") 
    }
  }
  
  # end of function declaration
  header[length(header)] <- ");"      
    
  if(curFun$dependency != ""){
    header <- c(paste("#ifdef ", curFun$dependency, sep = ""),
                header, 
                "#endif")  
  }
  
  write( c(header, "\n"), 
         file   = "AT_R_Wrapper.h", 
         append = TRUE)

    
  ###########################
  # B. create function body
  ###########################

  ##############################################
  # copy function declaration, replace semicolon
  body               <- header
  if(curFun$dependency != ""){                                      # if dependency is present, header is 
    body <- body[-length(body)]                                     # one line longer than expected
  }
  body[length(body)] <- "){"

  # get nature of parameters (in / out, array / no array, char)
  input.para.idx         <- grepl(pattern = "in",                   # parameter is input argument
                                  x       = para$in.out)
  output.para.idx        <- grepl(pattern = "out",                  # parameter is output argument
                                  x       = para$in.out)
  inoutput.para.idx      <- grepl(pattern = "in.out",               # parameter is both
                                  x       = para$in.out)

  array.idx              <- para$length != 1                        # parameter is array (i.e. length > 1)
  pointer.idx            <- grepl("*", para$type, fixed = TRUE)     # parameter is pointer
  char.idx               <- grepl("char", para$type)                # parameter is char
  length.by.variable.idx <- para$length%in%para$name                # array length for parameter is given by other parameter (e.g. "n")

  
  # add a count variable "i", if arrays are present
  if(sum(array.idx) > 0){
    body  <- c( body, 
                "  long i;")
  }

  ###########################################################
  # add type convert statements for non-array input variables
  # except chars
  idx <- which(input.para.idx & !array.idx & !char.idx)
  if(length(idx) > 0){
    for(j in idx){
    # j <- idx[1]
      	
      body <- c( body, 
                 paste( "  ", 
                        type.no.pointer( para$type[j] ), 
                        " ", 
                        para$name[j], 
                        get.extension( para$type[j]), 
                        " = (",
                        type.no.pointer.no.const ( para$type[j] ),
                        ")(*",  
                        para$name[j], 
			");", 
                        sep = ""))
    }
  }

  ###########################################################
  # add type convert statements for char variables
  #
  # chars can never be arrays in the C-wrapper
  # this has to be intercepted by loop on the R side
  idx <- which(input.para.idx & char.idx)
  if(length(idx) > 0){
    for(j in idx){
    # j <- idx[1]
      	
      body <- c( body, 
                 paste( "  ", 
                        type.no.pointer( para$type[j] ), 
                        "* ", 
                        para$name[j], 
                        get.extension( para$type[j]), 
                        " = (",
                        type.no.pointer.no.const ( para$type[j] ),
                        "*)(*",  
                        para$name[j], 
			");", 
                        sep = ""))
    }
  }
	
  ###########################################################
  # add conversion statements for arrays
  
  # also input parameters that are used to
  # give the size of input array, thus
  # change the parameter names for those
  # ! assumed that they are all int/long
  para$length[length.by.variable.idx]  <- paste( para$length[length.by.variable.idx], 
                                                 "_long", 
                                                 sep = "")

  # loop through all array lengths
  # this way, we can group arrays of same length in
  # one type casting loop
  array.lengths <- unique(para$length)[unique(para$length) != 1]
  for(array.length in array.lengths){
  # array.length <- array.lengths[1]
		
    #################
    # i. INPUT ARRAYS

    # find those input arrays that match in length
    cur.array.idxs <- which( input.para.idx & 
                             array.idx &
                             para$length == array.length)
    if(length(cur.array.idxs) > 0){
  	  
      # Allocate array space
      body <- c(body, "\n//Allocate space for the input parameter.")
      for(k in cur.array.idxs){
      # k <- cur.array.idxs[1]
        body <- c(body,
                  paste("  ", 
                        type.no.const( para$type[k] ), 
                        "* ",
				  		          para$name[k],
                        get.extension( para$type[k] ),
							          " = (",
                        type.no.const( para$type[k] ), 
                        "*)calloc(", 
                        array.length,  
                        ",sizeof(", 
                        type.no.const( para$type[k] ), 
                        "));", 
                        sep = ""))
			}

      body <- c(body, "")
			
      # Copy  data into the allocated space
      body <- c(body, 
                "\n//Fill in the input parameter.")
			body <- c(body, 
                paste( "  for(i = 0 ; i < ", 
                       array.length, 
                       "; i++){", 
                       sep = ""))
			for(k in cur.array.idxs){
				  body <- c(body, 
                    paste( "\t", 
                           para$name[k], 
                           get.extension( para$type[k] ),
						               "[i] = (", 
                           type.no.const( para$type[k] ),
                           ")", 
                           para$name[k],
						               "[i];",
                           sep = ""))
			}
			body <- c(body, "  }")
		  
		}

    ###################
    # ii. OUTPUT ARRAYS
    cur.array.idxs <- which( !input.para.idx & 
                             array.idx &
                             para$length == array.length)
		if(length(cur.array.idxs) > 0){
			body <- c(body, "\n//Allocate space for the results.")
			for(k in cur.array.idxs){
				body <- c( body, 
                   paste( "  ", 
                          type.no.const( para$type[k] ),
                          "* ",
						              para$name[k], 
                          get.extension(para$type[k]),
						              " = (", 
                          type.no.const(  para$type[k] ), 
						              "*)calloc(", 
                          array.length, 
                          ",sizeof(", 
                          type.no.const( para$type[k] ),
						              "));", 
                          sep = ""))
			}
		}
	}
      
  ###########################################################
  # add conversion statements for non-array pointers 
  # that are used as output-only variables
  idx.output.pointers     <- output.para.idx &
                             pointer.idx &
                             !char.idx &
                             !inoutput.para.idx
  if(sum(idx.output.pointers) > 0){
    output.pointers         <- para[idx.output.pointers, ]
    body                    <- c(body, "\n//Define type-casted output variables")
    for (j in 1:nrow(output.pointers)){
     # DEBUG: j <- 1
       body                    <- c( body, 
                                     paste( "\t", 
                                     type.no.pointer( output.pointers$type[j] ), 
                                     " ", 
                                     output.pointers$name[j], 
                                     get.extension( type.no.pointer( output.pointers$type[j] )),
				                             " = 0;", 
                                     sep = ""))
    }
  }


  ###########################
  # C. Generate function call
  ###########################

  # In case function a return value
  if(curFun$type == "void"){
  	return.var.txt	<-	""
  	para.max	      <-	nrow(para)
	}else{
  	return.var.txt	<-	paste( curFun$type, 
                               "returnValue_internal = \t")
  	para.max	      <-	nrow(para) - 1
	}

  # Parameter characteristics (char will be treated differently)
  char.but.no.pointer    <- input.para.idx &
                            char.idx
  
  array.or.char.pointer  <- array.idx |
                            !pointer.idx |
                            char.idx
    
  for(j in 1:para.max){
  # j <- 1
    
    # Write first line
    start.line <- "\t"
    if(j == 1){
      start.line <- paste( "\n  ", 
                           return.var.txt, 
                           curFun$name, 
                           "( ",
                           sep = "")
    }  
    
    # Write last line
    end.line <- ","
    if(j == para.max){
      end.line <- ");"
    }
     
    # Write arguments
    if(array.or.char.pointer[j]){                             # array
	body <- c( body,             
        paste( start.line,
               para$name[j], 
               get.extension(para$type[j]), 
               end.line,
               sep = ""))
    }else{						       # no array
	body <- c( body,                    
        paste( start.line,
               "&", 
               para$name[j], 
               get.extension( para$type[j] ), 
               end.line,
               sep = ""))
    } 			
  }
  
  ############################
  # Cast and copy back results
  ############################

  body <- c( body, 
             paste("\n//Results:"))
	
  # Cast return value, if any
  if(curFun$type != "void"){
		body <- c( body, 
               paste( "\n\t*returnValue = (", 
                      gsub( "*", 
                            "", 
                            gsub( "\t", 
                                  "", 
                                  convert.c(para$type[nrow(para)]), 
                                  fixed = TRUE), 
                            fixed = TRUE), 
                      ")returnValue_internal;"), 
                      sep = "")
	}
	
	# Cast non-array variables
  kk <-  which(output.para.idx & !array.idx & !char.idx)
	if(length(kk) > 0){
		for(j in kk){
			# j <- kk[1]
      body <- c( body, 
                 paste( "  *", 
                        para$name[j],	
                        " = (", 
					              gsub( "\t", 
                              "", 
                              gsub( "*", 
                                    "", 
                                    gsub( "const ", 
                                          "", 
                                          convert.c(para$type[j])), 
                                    fixed = TRUE), 
                              fixed = TRUE),			
					              ")", 
                        para$name[j], 
                        get.extension(para$type[j]), 
					              ";", 
                        sep = ""))
			body <- c( body, 
                 "")
		}
	}

	# Cast and copy arrays
  ll <-  which( output.para.idx & array.idx & !char.idx)
	if(length(ll) > 0){
		for(l in unique(para$length)){
		  kk <-  which( output.para.idx & 
                    array.idx & 
                    para$length == l)
		  if(length(kk) > 0){
				# fill the data into the allocated space
				body <- c( body, 
                   paste( "  for(i = 0 ; i < ", 
                          l, 
                          "; i++){", 
                          sep = ""))
				for(j in kk){
					body <- c( body, 
                     paste( "\t", 
                            para$name[j],	
                            "[i] = (", 
							              gsub( "\t", 
                                  "", 
                                  gsub( "*", 
                                        "", 
                                        gsub( "const ", 
                                              "", 
                                              convert.c(para$type[j])), 
                                        fixed = TRUE), 
                                  fixed = TRUE),			
							              ")", 
                            para$name[j], 
                            get.extension(para$type[j]),
							              "[i];", 
                            sep = ""))
				}
				body <- c( body, 
                   "  }")
			}
		}
	}
	
  # Free space allocated for arrays
	body <- c(body, "\n//Free allocated space")
	kk <- which(array.idx & !char.idx)	
	if(length(kk) > 0){
		for(j in kk){
				body <- c( body, 
                   paste( "  free(", 
                          para$name[j], 
                          get.extension(para$type[j]), 
                          ");", 
                          sep = ""))
    }
  }

  # write body to source file
  body <- c( body, 
             "}\n")
  
  if(curFun$dependency != ""){
    body <- c(body, 
              "#endif")
  }
  
  write( c(body, "\n"), 
         file   = "AT_R_Wrapper.c", 
         append = TRUE)

}# loop: functions from doxygen parsing

# write header to file
write( "#endif\n", 
       file    = "AT_R_Wrapper.h", 
       append  = TRUE)

