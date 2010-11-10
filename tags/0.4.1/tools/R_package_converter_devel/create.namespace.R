################################################################################################
# Namespace generator for the libamtrack R package
################################################################################################
# This script will crawl all header-files in the /include director of the libamtrack trunk
# and create a file containing the names of all functions found. The user can then decide
# which of them to include in the R package by deleting / commenting out individual
# function names.
################################################################################################
# Copyright 2006, 2010 The libamtrack team
# 
# This file is part of the AmTrack program (libamtrack.sourceforge.net).
#
#    Created on: 18.10.2010
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

# Clear workspace
rm(list = ls())

# Save current working directory
cur.dir <- getwd()

# Navigate to include path within libamtrack trunk, stop if fails
if (try(setwd("../../include")) != TRUE){
    stop("Please start script in /tools/R_package_generator")
}

# Get the names of all header files, exclude "AT_Wrapper_R.h"
names.all.files <- list.files(".")
names.all.header.files <- names.all.files[grep(".h", names.all.files)]
names.all.header.files <- names.all.header.files[-grep("AT_Wrapper_R.h", names.all.header.files)]

# initialize function name vector
function.names <- NULL

# Loop through all header files
for(name.cur.header.file in names.all.header.files){
      # DEBUG: name.cur.header.file <- names.all.header.files[7]
	
      # Read in current header file
      cur.header.file <- scan(name.cur.header.file, what = "character", sep = "\n", strip.white = FALSE)
      print(paste("Read: ", name.cur.header.file))

      # Find end of doxygen comments AND end of functions
	# The end of a doxygen comment ("*/") is not unique and can be
      # used for regular comments too. Thus, we use only those end positions that
      # lie before the end of a function declaration (");")
	end.doxygen.comments <- grep("*/", cur.header.file, fixed = T)
	end.functions <- grep (");", cur.header.file, fixed = T)

      # If function are found, extract names and store them. Otherwise skip current header file
      if(length(end.functions) > 0){
            cur.function.names <- NULL
            for (i in 1:length(end.functions)){
                  # DEBUG: i <- 1
                  pos <- max(end.doxygen.comments[end.doxygen.comments < end.functions[i]])
	            cur.function.names <- c(cur.function.names, cur.header.file[pos+1])
            }

            # extract function names. This should be in between a white space (after return
            # variable format and a opening bracket for the parameters, e.g.
            # double AT_example_function( const long param1, ...
            
            # exclude static variable definitions
            indices.static        <- grep("static", cur.function.names)
            if (length(indices.static > 0)){
                 cur.function.names    <- cur.function.names[-grep("static", cur.function.names)]}

            # remove "inline" statements
            indices.inline        <- grep("inline", cur.function.names)
            cur.function.names[indices.inline]    <- substring(  text   = cur.function.names[indices.inline], 
                                                                 first  = regexpr(" ", cur.function.names[indices.inline]) + 1, 
                                                                 last   = nchar(cur.function.names[indices.inline]))
            # then extract name
            cur.function.names    <- substring(  text   = cur.function.names, 
                                                 first  = regexpr(" ", cur.function.names) + 1, 
                                                 last   = regexpr("\\(", cur.function.names) - 1)
            print("Extracted functions:")
            print(cur.function.names)
      
            # append function names for current header file to function name vextor
            function.names <- c(function.names, cur.function.names)
      }else{
            print("Header file contains no functions.")
      }
}	

# Sort function names in alphabetical order
function.names <- function.names[order(function.names)]

# Restore working directory
setwd(cur.dir)

write( function.names, 
             file = "NAMESPACE.all",
             sep = "\n")

# End of script

