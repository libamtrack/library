################################################################################################
# C wrapper generator
################################################################################################
# Copyright 2006, 2009 Steffen Greilich / the libamtrack team
# 
# This file is part of the AmTrack program (libamtrack.sourceforge.net).
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


types <- c("bool", "long", "int", "long", "double")
conversion.c <- c("int*\t", "int*\t", "int*\t", "long*\t", "float*")
conversion.R <- c("as.integer", "as.integer", "as.integer", "as.integer", "as.single")
extension <- c("_long", "_long", "_long", "_long", "_double")

conversion.table <- data.frame(	type = c(types, paste("const", types), paste(types,"*",sep = "")),
						conversion.c = c(paste(conversion.c, "\t", sep = ""),
									 paste("const", conversion.c),
									 paste(conversion.c)),
						conversion.R = c(conversion.R, conversion.R, conversion.R),
						extension = c(extension, extension, extension),
						stringsAsFactors = F
					)

# add additional "\t\t" to "float*"
conversion.table$conversion.c[conversion.table$conversion.c == "float*\t"] <- "float*\t\t" 


convert.c <- function(type){
				ii <- type == conversion.table$type
#				cat("type=_",type,"_ ",conversion.table$type,ii,"\n",sep="")
				if( sum(ii) > 0 )
					return(conversion.table$conversion.c[ii])
				else
					return(type)				
			}

convert.R <- function(type){
				ii <- type == conversion.table$type
				return(conversion.table$conversion.R[ii])
			}