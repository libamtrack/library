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


types <- c("bool", "long", "int", "double", "char")
conversion.c <- c("int*\t", "int*\t", "int*\t", "float*", "char**")
conversion.R <- c("as.integer", "as.integer", "as.integer", "as.single", "as.character")
extension <- c("_bool", "_long", "_long", "_double", "_char")

conversion.table <- data.frame(	type = c(types, paste("const", types),  paste(types, "*", sep = "")),
						conversion.c = c(paste(conversion.c, "\t", sep = ""),
									paste("const", conversion.c),
									paste(conversion.c, "\t", sep = "")),
						conversion.R = c(conversion.R, conversion.R, conversion.R),
						extension = c(extension, extension, extension),
						stringsAsFactors = F
					)

# add additional "\t\t" to "float*"
conversion.table$conversion.c[conversion.table$conversion.c == "float*\t"] <- "float*\t\t" 


convert.c <- function(type){
				ii <- type == conversion.table$type
				if( sum(ii) > 0 )
					return(conversion.table$conversion.c[ii])
				else{
					print(paste(type, "is an unknown variable type!"))
					return(type)
				}				
			}

get.extension <- function(type){
				ii <- type == conversion.table$type
#				cat("type=_",type,"_ ",conversion.table$type,ii,"\n",sep="")
				if( sum(ii) > 0 )
					return(conversion.table$extension[ii])
				else{
					print(paste(type, "is an unknown variable type!"))
					return(type)
				}				
			}

convert.R <- function(type){
				ii <- type == conversion.table$type
				if( sum(ii) > 0 )
					return(conversion.table$conversion.R[ii])
				else{
					print(paste(type, "is an unknown variable type!"))
					return(type)
				}
			}

# Function that return types but without some selected features
type.no.pointer          <- function(x){
  return( gsub( "*", 
                "", 
                x, 
                fixed = TRUE))
}

type.no.const            <- function(x){
  return( gsub( "const ", 
                "", 
                x, 
                fixed = TRUE))
}

type.no.pointer.no.const <- function(x){
  return( gsub( "*", 
                "", 
                type.no.const(x),
                fixed = TRUE))
}
