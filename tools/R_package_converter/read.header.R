################################################################################################
# Script to extract all function declaration information for the libamtrack R package
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

# Save current working directory
cur.dir <- getwd()

# Read in namespace file
namespace <- scan(file = "NAMESPACE", what = "character")

# Navigate to include path within libamtrack trunk. Stop if this fails
if(try(setwd("../../include")) == FALSE){
     stop("Please start script from /tools/R_package_converter")
}

# Read in all header files
record <- list.files(".")
record <- record[grep(".h", record)]

# Remove the old style R-wrappers from the list of header files - they are obsolete here
to.remove <- c("AT_Wrapper_R.h")
pos.remove <- match(to.remove, record)
pos.remove <- c(pos.remove, grep("Data", record))
record <- record[-pos.remove]

# Initialize vectors to hold extracted information
current.function <- NULL
functions <- NULL
processed.functions <- NULL
function.no <- 1

# Replacement for "grepl" function to ensure compatibilty with R <= 2.9.0
grep.bool     <-     function(pattern, x, ...){
     results     <-     grep(pattern, x, ...)
     indices     <-     1:length(x)
     bool.res     <-     is.element(indices, results)
}

for(file in record){
     # DEBUG: file <- record[4]
     
     # read the complete *.h file
     include <- scan(file, what = "character", strip.white = T)
     print(paste("Read: ", file))

     # just keep the code of current.functions passed to R
     start.current.functions <- grep("/**", include, fixed = T)
     end.current.functions <- grep (");", include, fixed = T)

     # break and process next file in case no function are found
     if(length(start.current.functions) == 0 | length(end.current.functions) == 0){
          print("Skipping file - no functions found.")
          next
     }

     # remove the first entry of start.curent.functions from the doxygen file comment
     start.current.functions <- start.current.functions[2:length(start.current.functions)]
     
     # check if start.funtions has same length as end.current.functions
     keep <- NULL
     if(length(start.current.functions) == length(end.current.functions)){
          for(i in 1:length(start.current.functions)){
               keep <- c(keep, (start.current.functions[i]):(end.current.functions[i])) 
          }     
     } else if (length(start.current.functions) < length(end.current.functions)){
          # go through "/**" and just select the next ");" 
          j = 1
          for(i in 1:length(start.current.functions)){
               while(start.current.functions[i] > end.current.functions[j]) j = j+1
               keep <- c(keep, (start.current.functions[i]):(end.current.functions[j])) 
          }
     } else {
          # go through ");" and select the previous "/**" 
          j = length(start.current.functions)
          for(i in length(end.current.functions):1){
               while(end.current.functions[i] < start.current.functions[j]) j = j-1
               keep <- c(keep, (start.current.functions[j]):(end.current.functions[i])) 
          }
     }
     include <- include[keep]

     #############################################
     # separate the doxygen comment and the code
     #############################################

     # all vectors should be of same length
     comment.start <- grep("/**", include, fixed = T)
     comment.end <- grep("*/", include, fixed = T)
     current.function.end <- grep(";", include)

     for (i in 1:length(comment.start)){
          # DEBUG: i <- 1

          pos.comment <- comment.start[i]:comment.end[i] 
          pos.code <- (comment.end[i]+1):current.function.end[i]

          # extract the name, comment and the code
          name <- gsub("(", "", include[pos.code[2]], fixed = T)
          print(paste("Processing:", name))

          # check if name is in namespace before processing
          if(name %in% namespace){

               raw.comment <- include[pos.comment]
               raw.code <- include[pos.code]
          
               ########################
               # process the comment
               ########################
               current.function$name <- name     

               # get the "*" position
               breaks <- grep("*", raw.comment, fixed = T)

               # remove the first "/**" entry
               breaks <- breaks[-1] 
          
               # extract the description (text before the first parameter)
               first.parameter <- grep("@", raw.comment, fixed = T)

               if(length(first.parameter) > 0){
                    pos.description <- (3:(first.parameter[1]-2))
          
                    tmp <- raw.comment[pos.description[1]]
                    if(length(pos.description) >= 2){
                         for(k in 2:(length(pos.description))){
                              # DEBUG: k <- 2
                              if(raw.comment[pos.description[k]] != "*"){
                                   tmp <- paste(tmp, raw.comment[pos.description[k]])
                              }else{
                                   tmp <- paste(tmp, "\n", sep = "")}
                         }
                    }
 
                    current.function$description <- tmp          

                    # extract the parameters
                    breaks <- breaks[breaks > (first.parameter[1] - 2)]

                    # create data.frame to hold the parameter data
                    parameter <- data.frame( type = character(length(breaks) - 1),
                                             name = character(length(breaks) - 1),
                                             comment = character(length(breaks) - 1),
                                             array.size = 1,
                                             stringsAsFactors = F)          

                    for(j in 1:(length(breaks)-1)){
                         pos <- (breaks[j] + 1):(breaks[j+1] - 1)
                         parameter$type[j] <- raw.comment[pos[1]]
                         parameter$name[j] <- raw.comment[pos[2]]     

                         if(length(pos) > 2){
                              # paste the comment together
                              tmp <- NULL
                                   for(k in pos[3]:pos[length(pos)]){
                                   tmp <- paste(tmp, raw.comment[k])
                              }
                              parameter$comment[j] <- tmp     
                              }

                         if("array" %in% raw.comment[pos] | "(array" %in% raw.comment[pos]) {
                              array.size.pos1 <- grep("array", raw.comment[pos])
                              array.size.pos2 <- grep("size", raw.comment[pos])
                              if((array.size.pos1 +2) %in% array.size.pos2)      parameter$array.size[j] <- gsub(")", "", raw.comment[pos][array.size.pos2 + 1])
                         }
                    }
                    current.function$parameter.comment <- parameter
               } else current.function$parameter.comment <- "empty"

               ####################################
               # process the code
               ####################################

               current.function$type <- raw.code[1]
          
               # remove the first 2 entries
               raw.code <- raw.code[-c(1,2)]

               # get the line breaks
               breaks <- grep(",", raw.code, fixed = T)

               # remove additional comments in the code "//...."
               to.remove <- grep("//", raw.code, fixed = T)
               if(length(to.remove) > 0)
                    for(k in 1:length(to.remove)){
                         remove <- to.remove[k]:((breaks[breaks > to.remove[k]])[1] - 2)
                         if(raw.code[length(remove)] == "const") remove <- remove[1:(length(remove)-1)]
                         check <- raw.code[remove]
                         raw.code <- raw.code[-remove]

                    # adjust breaks by shifting of length(remove)
                    breaks[breaks > to.remove[k]] <- breaks[breaks > to.remove[k]] - length(remove)
               } # hier weiter

               # remove "," 
               raw.code <- gsub(",", "", raw.code)

               # create data.frame to hold the parameter data
               parameter <- data.frame(type = character(length(breaks) + 1),
                                   name = character(length(breaks) + 1),
                                   length = rep(1, length(breaks) + 1),
                                   stringsAsFactors = F )          


               # first position
               pos <- 1:(breaks[1])
               if(length(pos) == 2){
                    parameter$type[1] <- raw.code[pos[1]]
                    parameter$name[1] <- raw.code[pos[2]]
               } else if(length(pos) > 2) {
                    parameter$type[1] <- paste(raw.code[pos[1]], raw.code[pos[2]])
                    parameter$name[1] <- raw.code[pos[3]]
               }

               # second to second to last position
               if(length(breaks) >= 2) for(j in 2:(length(breaks))){
                    pos <- (breaks[j - 1] + 1):(breaks[j])
                    if(length(pos) == 2){
                         parameter$type[j] <- raw.code[pos[1]]
                         parameter$name[j] <- raw.code[pos[2]]
                    } else if(length(pos) > 2) {
                         parameter$type[j] <- paste(raw.code[pos[1]], raw.code[pos[2]])
                         parameter$name[j] <- raw.code[pos[3]]
                    }
               }

               # last position
               pos <- (breaks[length(breaks)] + 1):length(raw.code)
               if(length(pos) == 2){
                    parameter$type[length(breaks) + 1] <- raw.code[pos[1]]
                    parameter$name[length(breaks) + 1] <- gsub(");", "", raw.code[pos[2]])
               } else if(length(pos) == 3) {
                    parameter$type[length(breaks) + 1] <- paste(raw.code[pos[1]], raw.code[pos[2]])
                    parameter$name[length(breaks) + 1] <- gsub(");", "", raw.code[pos[3]])
               }

               # find the vectors and remove "[x]" from the name
               vectors <- grep("[", parameter$name, fixed = T)
               if(length(vectors) > 0){               
                    parameter$name[vectors] <- unlist(strsplit(parameter$name[vectors], "[", fixed = T))[(seq(2, length(vectors)*2, by = 2) - 1)]     
               
                    # get the array.size from the doxgen comment
                    pos.array.size <- match(parameter$name[vectors], current.function$parameter.comment$name)
                    parameter$length[vectors] <- current.function$parameter.comment$array.size[pos.array.size]
               }
               
               # remove possible NA from parameter$length by TODO
               parameter$length[is.na(parameter$length)] <- "TODO"

               # add parameters to list
               current.function$parameter <- parameter

               # add input output information from doxygen file 
               pos.in <- grep.bool(pattern = "in", x = current.function$parameter.comment$type) 
               pos.out <- grep.bool(pattern = "out", x = current.function$parameter.comment$type) 

               # get the position of the parameters from the position in the comments
               pos.in.para <- match(current.function$parameter.comment$name[pos.in & !pos.out], current.function$parameter$name)
               pos.out.para <- match(current.function$parameter.comment$name[pos.out & !pos.in], current.function$parameter$name)
               pos.in.out.para <- match(current.function$parameter.comment$name[pos.in & pos.out], current.function$parameter$name)

               current.function$parameter$in.out <- "TODO"
               current.function$parameter$in.out[pos.in.para] <- "in" 
               current.function$parameter$in.out[pos.out.para] <- "out" 
               current.function$parameter$in.out[pos.in.out.para] <- "in.out" 

               # SG: add return parameters if existing
               if(current.function$type != "void"){
                    if(length(grep("return", current.function$parameter$name))>0){
                         print("!ERROR: name conflict with return variable!")
                    }
                    current.function$parameter     <-     rbind.data.frame(     current.function$parameter,
                                                                 data.frame(     type     =  include[pos.code[1]],
                                                                           name     =  "returnValue",
                                                                           length = 1,
                                                                           in.out = "return"))
               }

               functions[[function.no]] <- current.function
               processed.functions <- c(processed.functions, name)
               function.no <- function.no + 1

          }else{ # name %in% namespace
               print(paste("Skipped as not in namespace.", name))
          }
     }
}     

# restore working directory
setwd(cur.dir)

save(functions, file = "functions.ssd")
