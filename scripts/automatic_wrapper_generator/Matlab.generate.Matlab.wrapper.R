################################################################################################
# Matlab wrapper generator
################################################################################################
# Copyright 2006, 2011 The libamtrack team
# 
# This file is part of the AmTrack program (libamtrack.sourceforge.net).
#
#    Created on: 21.03.2011
#    Creator: greilich
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

# replacement for "grepl" function to ensure compatibilty with R <= 2.9.0
grep.bool	<-	function(pattern, x, ...){
	results	<-	grep(pattern, x, ...)
	indices	<-	1:length(x)
	bool.res	<-	is.element(indices, results)
	return(bool.res)
}

for(i in 1:length(functions)){
	#i<-1
	tmp <- functions[[i]]
	file.name <- paste(tmp$name, "m", sep = ".")

	write("% Automatically created wrapper file", file = file.name)
	
	para   	<- tmp$parameter[grep.bool(pattern = "in", x = tmp$parameter$in.out)]
	n.para 	<- nrow(para)
	para.string <- ""
	help.string <- ""
	for(j in 1:n.para){
		if (j != 1){
			separator = ", "
		}else{
			separator = ""
		}
		para.string <- paste(para.string, para$name[j], sep = separator)
		help.string <- paste(help.string, "% ", para$name[j], ": ", tmp$parameter.comment$comment[j], "\n", sep = "")
	}
	
	write(paste("function ", tmp$name, "(", para.string, ")", sep = ""), file = file.name, append = TRUE)
	
	# Add help
	write(paste("% ", gsub("\n", "", tmp$description), sep = ""),					file = file.name, append = TRUE)
	write("%\n% List on input parameters\n% -------------------------",	file = file.name, append = TRUE)
	write(help.string,									file = file.name, append = TRUE)


	# Add loading of library if not yet done
	write("\nif not(libisloaded('libamtrack'))", 	file = file.name, append = TRUE)
	write("\thfile = \'libamtrack.h\';", 		file = file.name, append = TRUE)
	write("\tloadlibrary(\'libamtrack\', hfile);",	file = file.name, append = TRUE)
	write("end\n", 						file = file.name, append = TRUE)
  
	# Add function call
	write(paste("\ncalllib(\'libamtrack\', ", tmp$name, ", ", tmp$parameter$name[j], ", ", para.string, ")", sep = ""), file = file.name, append = TRUE)

	# End
	write("\nend", file = file.name, append = TRUE)
}
