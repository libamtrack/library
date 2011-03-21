################################################################################################
# Matlab header wrapper generator
################################################################################################
# Copyright 2006, 2011 The libamtrack team
# 
# This file is part of the AmTrack program (libamtrack.sourceforge.net).
#
#    Created on: 3.2.2011
#    Creator: sgreilich
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

write("#ifndef LIBAMTRACK_H_\n#define LIBAMTRACK_H_\n// Automatically created header file for use of libamtrack in Matlab\n\n", file = "libamtrack.h")

load("functions.sdd")

for(i in 1:length(functions)){
	# i <- 1
	
	#####################################
	# copy function comment & declaration
	#####################################

	comment          <-  functions[[i]]$raw.comment.text
      declaration      <-  functions[[i]]$raw.declaration.text

	write(c(comment, declaration, "\n\n"), file = "libamtrack.h", append = T)
}

write("#endif\n", file = "libamtrack.h", append = T)