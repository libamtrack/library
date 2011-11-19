################################################################################################
# Script to extract all function declaration information for the libamtrack R package
################################################################################################
# Copyright 2006, 2010 The libamtrack team
# 
# This header.file.name is part of the AmTrack program (libamtrack.sourceforge.net).
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
# long with AmTrack (header.file.name: copying.txt).
# If not, see <http://www.gnu.org/licenses/>
################################################################################################

# Clean workspace
rm(list = ls())

# Save current working directory
cur.dir           <- getwd()

# Read in namespace header.file.name
namespace         <- scan(file = "NAMESPACE", what = "character", sep = "\n")

# find those function marked with "noR: " for which no R wrapper,
# only C wrapper will be created
noR               <- grepl("noR: ", namespace)
namespace         <- gsub("^noR: ", "", namespace)

# Navigate to include path within libamtrack trunk. Stop if this fails
if((try(setwd("../../../include")) == FALSE)&(try(setwd("../../include")) == FALSE)){
    stop("I do not see from current folder /include. Please move to other place !")
}

# Read in all header file names
header.file.names <- list.files(".")
header.file.names <- header.file.names[grep(".h", header.file.names)]

# Remove the old style R-wrappers from the list of header files - they are obsolete here. 
# Additionally header files for which doxygen comments are not ready yet can be excluded
to.remove         <- c("AT_NumericalRoutines.h")
pos.remove        <- match(to.remove, header.file.names)
#pos.remove        <- c(pos.remove, grep("Data", header.file.names))
header.file.names <- header.file.names[-pos.remove]

# Initialize vectors to hold extracted information
functions             <- NULL
processed.functions   <- NULL
function.no           <- 1

# Replacement for "grepl" function to ensure compatibilty with R <= 2.9.0
grep.bool     <-     function(pattern, x, ...){
     results     <-     grep(pattern, x, ...)
     indices     <-     1:length(x)
     bool.res    <-     is.element(indices, results)
}

for(header.file.name in header.file.names){
     # DEBUG: header.file.name <- header.file.names[21]
     
     #####################################################
     # A. Extract doxygen comment and function declaration
     #####################################################

     # Read the current header file (separated into 'words' (=chunks of strings that were separated by whitespaces, line
     # break etc. All whitespaces are removed during read in
     header.file.text                 <- scan(header.file.name, what = "character", strip.white = T)
     
     # Read the current header file (in lines)
     raw.header.file.text             <- scan(header.file.name, what = "character", sep = "\n")

     # Read the corresponding source file (in lines)
     source.file.name                 <- gsub(".h$", ".c", header.file.name)
     source.file.text                 <- scan( paste("../src/", source.file.name, sep = ""),
                                               what = "character", 
                                               sep = "\n", 
                                               blank.lines.skip = FALSE)

     print(paste("Read: ", header.file.name))

     # Find the start of doxygen comments ("/**") and the end of a declaration (");")
     # They should embrase a function declaration (but also others, like enumerators etc.)
     pos.start.doxygen.comment        <- grep("/**", header.file.text, fixed = T)
     pos.end.function.declaration     <- grep (");", header.file.text, fixed = T)
     raw.pos.start.doxygen.comment    <- grep("/**", raw.header.file.text, fixed = TRUE)
     raw.pos.end.function.declaration <- grep(");", raw.header.file.text, fixed = TRUE)

     # Break and process next header file in case no function are found
     if(length(pos.start.doxygen.comment) == 0 | length(pos.end.function.declaration) == 0){
          print("Skipping header file - no functions found.")
          next
     }

     # Remove the first entry of pos.start.doxygen.comment as this is (most likely!) the
     # brief description of the header file itself (@brief)
     # TODO: replace by more robust statement which looks for "@brief"
     pos.start.doxygen.comment     <- pos.start.doxygen.comment[-1]
     raw.pos.start.doxygen.comment <- raw.pos.start.doxygen.comment[-1]
     
     # Check whether both vectors have same length
     # If not, remove all those comment-starts that do not correlate with 
     # a function-end
     # first, initialize keep-variable which will contain all element indices of header.file.text to keep
     idx.keep     <- NULL
     raw.idx.keep <- NULL
     # If both vectors have the same length, keep everything but reformat
     if(length(pos.start.doxygen.comment) == length(pos.end.function.declaration)){
          for(i in 1:length(pos.start.doxygen.comment)){
               idx.keep     <- c(idx.keep,     (pos.start.doxygen.comment[i]):(pos.end.function.declaration[i])) 
               raw.idx.keep <- c(raw.idx.keep, (raw.pos.start.doxygen.comment[i]):(raw.pos.end.function.declaration[i])) 
          }     
     # if more function ends, just select those that correspond to a comment-start
     } else if (length(pos.start.doxygen.comment) < length(pos.end.function.declaration)){
          j <- 1
          for(i in 1:length(pos.start.doxygen.comment)){
               while(pos.start.doxygen.comment[i] > pos.end.function.declaration[j]){
                    j = j+1
               }
               idx.keep      <- c(idx.keep,     (pos.start.doxygen.comment[i]):(pos.end.function.declaration[j]))
               raw.idx.keep  <- c(raw.idx.keep, (raw.pos.start.doxygen.comment[i]):(raw.pos.end.function.declaration[j]))
          } 
     # if more comment starts, just select those that correspond to a function end
     } else {
          j <- length(pos.start.doxygen.comment)
          for(i in length(pos.end.function.declaration):1){
               while(pos.end.function.declaration[i] < pos.start.doxygen.comment[j]){
                    j = j-1
               }
               idx.keep     <- c(idx.keep,     (pos.start.doxygen.comment[j]):(pos.end.function.declaration[i])) 
               raw.idx.keep <- c(raw.idx.keep, (raw.pos.start.doxygen.comment[j]):(raw.pos.end.function.declaration[i])) 
          }
     }

     # Now, only the comment and declaration text of suitable functions are left
     reduced.header.file.text           <- header.file.text[idx.keep]
     raw.reduced.header.file.text       <- raw.header.file.text[raw.idx.keep]
     


     #####################################################
     # B. Process doxygen comment and function declaration
     #####################################################

     # Find start and end of doxygen comments in reduced header file
     # as well as end of declarations
     pos.start.doxygen.comment           <- grep("/**", reduced.header.file.text, fixed = T)
     pos.end.doxygen.comment             <- grep("*/", reduced.header.file.text, fixed = T)
     pos.end.function.declaration        <- grep(");", reduced.header.file.text, fixed = T)

     raw.pos.start.doxygen.comment       <- grep("/**", raw.reduced.header.file.text, fixed = TRUE)
     raw.pos.end.doxygen.comment         <- grep("*/", raw.reduced.header.file.text, fixed = TRUE)
     raw.pos.end.function.declaration    <- grep(");", raw.reduced.header.file.text, fixed = TRUE)

     # If length are different comments are corrupted, skip file and process next
     if((length(pos.start.doxygen.comment) != length(pos.end.doxygen.comment))|(length(pos.start.doxygen.comment) != length(pos.end.function.declaration))){
          print("Comments corrupted, skipping file.")
          next
     }

     # Loop through all functions found
     for (i in 1:length(pos.start.doxygen.comment)){
          # DEBUG: i <- 2
          current.function    <- NULL
          
          # Indices for reduced.header.file.text that hold comment and declaration
          idx.comment         <- pos.start.doxygen.comment[i]:pos.end.doxygen.comment[i] 
          idx.declaration     <- (pos.end.doxygen.comment[i]+1):pos.end.function.declaration[i]
          raw.idx.comment     <- raw.pos.start.doxygen.comment[i]:raw.pos.end.doxygen.comment[i] 
          raw.idx.declaration <- (raw.pos.end.doxygen.comment[i]+1):raw.pos.end.function.declaration[i]

          # Extract the function name (which is on second position of the declaration, opening bracket has to be removed
          name                <- gsub("(", "", reduced.header.file.text[idx.declaration[2]], fixed = T)
          print(paste("Processing:", name))

          # Check if name is in namespace otherwise skip
          if(name %in% namespace){

               # Extract comment and declaration text
               current.function$raw.comment.text        <- raw.reduced.header.file.text[raw.idx.comment]
               current.function$raw.declaration.text    <- raw.reduced.header.file.text[raw.idx.declaration]
               raw.comment.text                         <- reduced.header.file.text[idx.comment]
               raw.declaration.text                     <- reduced.header.file.text[idx.declaration]
          
               # Function name
               current.function$name                    <- name     

               # Store line in header file (for links to code in documentation)
               # TODO: To take the first appearance of the function name as the actual
               # source code and not a call to the function (grep will find both)
               # is a pretty sloppy approach and should be improved
               current.function$src.line.no             <- grep(name, source.file.text, fixed = TRUE)[1]
               current.function$src.file.name           <- source.file.name

               ########################################
               # B1. Process doxygen parameter comments
               ########################################

               # Get start positions of doyxgen comment lines in comment text
               # by grepping "*" 
               pos.comment.lines       <- grep("*", raw.comment.text, fixed = T)

               # Remove the first line ("/**")
               pos.comment.lines       <- pos.comment.lines[-1] 
          
               # Find start positions of doxygen entries ("@")
               pos.doxygen.entries     <- grep("@", raw.comment.text, fixed = T)

               # If no entries found, skip processing and store empty parameter comments
               # TODO: Check whether these are "@param" entries
               # TODO: The function comment should be processed even if there is no @ entry
               if(length(pos.doxygen.entries) > 0){

                    # The description of the function itself is located 
                    idx.doxygen.parameter.description    <- (3:(pos.doxygen.entries[1]-2))
                    # Walk through function description text and remove comment indicators ("*") but keep line breaks
                    function.description                 <- raw.comment.text[idx.doxygen.parameter.description[1]]
                    if(length(idx.doxygen.parameter.description) >= 2){
                         for(idx in 2:(length(idx.doxygen.parameter.description))){
                              if(raw.comment.text[idx.doxygen.parameter.description[idx]] != "*"){
                                   function.description <- paste(  function.description, 
                                                                   raw.comment.text[idx.doxygen.parameter.description[idx]])
                              }else{
                                   function.description <- paste(  function.description, 
                                                                   "\n", sep = "")}
                         }
                    }
                    # Store the description
                    current.function$description <- function.description          


                    # Extract the parameter comments
                    pos.parameter.comment.lines  <- pos.comment.lines[pos.comment.lines > (pos.doxygen.entries[1] - 2)]

                    # Create data.frame to hold the parameter information from doxygen comment
                    # As last entry is "@return", subtract one position
                    # TODO: Replace guess work on @doxygen entries by clear names (@param, @return, etc.)
                    parameter                    <- data.frame( type              = character(length(pos.parameter.comment.lines) - 1),
                                                                name              = character(length(pos.parameter.comment.lines) - 1),
                                                                comment           = character(length(pos.parameter.comment.lines) - 1),
                                                                array.size        = 1,
                                                                stringsAsFactors  = F)          

                    # Only process if number of parameter comments is >= 2 (i.e. one @param and one @return (always there))
                    # TODO: Awkward
                    if(length(pos.parameter.comment.lines) >= 2){
                         
                         # Loop through all @param comments, leave out @return (last one)
                         for(j in 1:(length(pos.parameter.comment.lines)-1)){
                              # DEBUG: j <- 1
                              
                              # Vector with all indices in raw.comment.text that belong to current @param comment
                              idx.param.comment <- (pos.parameter.comment.lines[j] + 1):(pos.parameter.comment.lines[j+1] - 1)
                              # Comment type (@param etc.) is first entry
                              parameter$type[j] <- raw.comment.text[idx.param.comment[1]]
                              # Parameter name is second entry
                              # TODO: MAYBE SECOND ENTRY...
                              parameter$name[j] <- raw.comment.text[idx.param.comment[2]]     

                              # Extract the actual comment (third to last position), paste together from text chunks
                              # Enter loop only if comment exists
                              if(length(idx.param.comment) > 2){
                                   parameter.comment <- NULL
                                   for(k in idx.param.comment[3]:idx.param.comment[length(idx.param.comment)]){
                                        parameter.comment <- paste(parameter.comment, raw.comment.text[k])
                                   }
                                   # Store comment
                                   parameter$comment[j] <- parameter.comment
                              }

                              # If key-word array appears in parameter comment
                              # try to find size information
                              if(("array" %in% raw.comment.text[idx.param.comment])|("(array" %in% raw.comment.text[idx.param.comment])){
                                   array.size.pos1       <- grep("array", raw.comment.text[idx.param.comment])
                                   array.size.pos2       <- grep("size", raw.comment.text[idx.param.comment])
                                   if((array.size.pos1 + 2) %in% array.size.pos2){
                                        parameter$array.size[j] <- gsub(")", "", raw.comment.text[idx.param.comment][array.size.pos2 + 1])
                                   }
                              }
                         } # @param loop
                         current.function$parameter.comment <- parameter
                    }else{ # length(pos.parameter.comment.lines) >= 2
                         current.function$parameter.comment <- "empty" # no parameter comments
                    }
               }else{ #length(pos.doxygen.entries) > 0
                    current.function$parameter.comment <- "empty"
               }

               ######################################
               # B2. Process the function declaration
               ######################################

               # Return parameter type, first position of declaration text
               current.function$type              <- raw.declaration.text[1]
          
               # remove the first 2 entries, i.e. return type and function name
               raw.declaration.text               <- raw.declaration.text[-c(1,2)]

               # Extract positions of declaration line breaks
               pos.declaration.line.breaks        <- grep(",", raw.declaration.text, fixed = T)

               # remove additional comments in the code "//...."
               if(length(pos.declaration.line.breaks) > 0){
                    to.remove                          <- grep("//", raw.declaration.text, fixed = T)
                    if(length(to.remove) > 0){
                         for(k in 1:length(to.remove)){
                              remove                        <- to.remove[k]:((pos.declaration.line.breaks[pos.declaration.line.breaks > to.remove[k]])[1] - 2)
                              if(raw.declaration.text[length(remove)] == "const"){
                                   remove <- remove[1:(length(remove)-1)]
                              }
                              check                         <- raw.declaration.text[remove]
                              raw.declaration.text          <- raw.declaration.text[-remove]
                         }
                         # adjust pos.declaration.line.breaks by shifting of length(remove)
                        pos.declaration.line.breaks[pos.declaration.line.breaks > to.remove[k]] <- pos.declaration.line.breaks[pos.declaration.line.breaks > to.remove[k]] - length(remove)
                    }
               }

               # remove line breaks (",")
               raw.declaration.text <- gsub(",", "", raw.declaration.text)

               # create data.frame to hold the parameter data
               parameter            <- data.frame(type             = character(length(pos.declaration.line.breaks) + 1),
                                                  name             = character(length(pos.declaration.line.breaks) + 1),
                                                  length           = rep(1, length(pos.declaration.line.breaks) + 1),
                                                  stringsAsFactors = F )          

               # a. First parameter
               if(length(pos.declaration.line.breaks) > 0){
                    pos <- 1:(pos.declaration.line.breaks[1])
                    if(length(pos) == 2){
                         parameter$type[1] <- raw.declaration.text[pos[1]]
                         parameter$name[1] <- raw.declaration.text[pos[2]]
                    } else if(length(pos) > 2) {
                         parameter$type[1] <- paste(raw.declaration.text[pos[1]], raw.declaration.text[pos[2]])
                         parameter$name[1] <- raw.declaration.text[pos[3]]
                    }

                    # b. Second to second to last parameter
                    if(length(pos.declaration.line.breaks) >= 2) for(j in 2:(length(pos.declaration.line.breaks))){
                         pos <- (pos.declaration.line.breaks[j - 1] + 1):(pos.declaration.line.breaks[j])
                         if(length(pos) == 2){
                              parameter$type[j] <- raw.declaration.text[pos[1]]
                              parameter$name[j] <- raw.declaration.text[pos[2]]
                         } else if(length(pos) > 2) {
                              parameter$type[j] <- paste(raw.declaration.text[pos[1]], raw.declaration.text[pos[2]])
                              parameter$name[j] <- raw.declaration.text[pos[3]]
                         }
                    }

                    # c. Last parameter
                    pos <- (pos.declaration.line.breaks[length(pos.declaration.line.breaks)] + 1):length(raw.declaration.text)
                    if(length(pos) == 2){
                         parameter$type[length(pos.declaration.line.breaks) + 1] <- raw.declaration.text[pos[1]]
                         parameter$name[length(pos.declaration.line.breaks) + 1] <- gsub(");", "", raw.declaration.text[pos[2]])
                    } else if(length(pos) == 3) {
                         parameter$type[length(pos.declaration.line.breaks) + 1] <- paste(raw.declaration.text[pos[1]], raw.declaration.text[pos[2]])
                         parameter$name[length(pos.declaration.line.breaks) + 1] <- gsub(");", "", raw.declaration.text[pos[3]])
                    }
               }else{ # if no line break, the second last entry *should* be the type, the last the name
                    parameter$name[1] <-raw.declaration.text[length(raw.declaration.text)-1]
                    parameter$type[1] <- gsub(");", "", raw.declaration.text[length(raw.declaration.text)])
               }

               # find the vectors and remove "[x]" from the name
               vectors        <- grep("[", parameter$name, fixed = T)
               if(length(vectors) > 0){               
                    parameter$name[vectors]   <- unlist(strsplit(parameter$name[vectors], "[", fixed = T))[(seq(2, length(vectors)*2, by = 2) - 1)]     
                    # get the array.size from the doxgen comment
                    pos.array.size            <- match(parameter$name[vectors], current.function$parameter.comment$name)
                    parameter$length[vectors] <- current.function$parameter.comment$array.size[pos.array.size]
               }
               
               # remove possible NA from parameter$length by TODO
               parameter$length[is.na(parameter$length)]          <- "TODO"

               # add parameters to list
               current.function$parameter  <- parameter

               # add input output information from doxygen header.file.name 
               pos.in                      <- grep.bool(pattern = "in", x = current.function$parameter.comment$type) 
               pos.out                     <- grep.bool(pattern = "out", x = current.function$parameter.comment$type) 

               # get the position of the parameters from the position in the comments
               pos.in.para                 <- match(current.function$parameter.comment$name[pos.in & !pos.out], current.function$parameter$name)
               pos.out.para                <- match(current.function$parameter.comment$name[pos.out & !pos.in], current.function$parameter$name)
               pos.in.out.para             <- match(current.function$parameter.comment$name[pos.in & pos.out], current.function$parameter$name)

               current.function$parameter$in.out                  <- "out"
               current.function$parameter$in.out[pos.in.para]     <- "in" 
               current.function$parameter$in.out[pos.out.para]    <- "out" 
               current.function$parameter$in.out[pos.in.out.para] <- "in.out" 

               # SG: add return parameters if existing
               if(current.function$type != "void"){
                    if(length(grep("return", current.function$parameter$name))>0){
                         print("!ERROR: name conflict with return variable!")
                    }
                    current.function$parameter     <-     rbind.data.frame(     current.function$parameter,
                                                                                data.frame(  type     =  reduced.header.file.text[idx.declaration[1]],
                                                                                             name     =  "returnValue",
                                                                                             length   = 1,
                                                                                             in.out   = "return"))
               }

               ##############################
               # Detect array size variables
               # for automatic array size
               # definition in R wrapper
               ##############################

               current.function$parameter$array.size.defined.by.variable.no     <- 0
               current.function$parameter$derivable.array.size.variable         <- FALSE

               for(j in 1:nrow(current.function$parameter)){
                  # j <- 12
                  if(length(regexpr("[[:alpha:]]", current.function$parameter$length[j]))>0){                   # only if non-digit character in length: numbers can never be variables
                     idx <- grep(paste( "^", current.function$parameter$length[j], "$", sep = ""), 
                                        current.function$parameter$name)                                        # look for exact match of length name in variable list
                     if(length(idx) == 1 & (current.function$parameter$in.out[j] == "in" | current.function$parameter$in.out[j] == "in.out")){
                          current.function$parameter$array.size.defined.by.variable.no[j] <- idx
                     }
                  }
               }
               for(j in 1:nrow(current.function$parameter)){
                  # j <- 1
                  idx    <- grep(paste("^", current.function$parameter$name[j], "$", sep = ""), current.function$parameter$length)  # use exact matching here to avoid ambiguities
                  if(length(idx) > 0){
					  # check if any of the dependent variables is input, otherwise skip automatic array size determination!
					  ii     <- (current.function$parameter$in.out[idx] == "in") | (current.function$parameter$in.out[idx] == "in.out")
					  if( sum(ii) > 0){
						  current.function$parameter$derivable.array.size.variable[j] <- TRUE
					  }
                  }
               }

               ################################
               # Detect input array that should 
               # have a fixed size > 1
               ################################

               current.function$parameter$fixed.size.array   <- FALSE

               for(j in 1:nrow(current.function$parameter)){
                  # j <- 1
                  if(regexpr("[[:alpha:]]", current.function$parameter$length[j]) == -1){                   # only length with digits
				                if(as.numeric(current.function$parameter$length[j])>1){
                             if(current.function$parameter$in.out[j] == "in"){
                                  current.function$parameter$fixed.size.array[j]     <- TRUE
                             }
                        }
                  }
               }
 
               # Add flag if both R and C wrappers or C wrappers only will be built
               current.function$wrapper.type <- "R and C"
               if(noR[match(name,namespace)]){
                  current.function$wrapper.type <- "C only"
               }
               
               # Add header file name for later include
               current.function$header.file.name <- header.file.name
               
               # Store processed function data
               functions[[function.no]]    <- current.function
               processed.functions         <- c(processed.functions, name)
               function.no                 <- function.no + 1
               print(paste("Successfully done.", name))
        }else{ # name %in% namespace
               print(paste("Skipped as not in namespace.", name))
        }
     } # end function loop --- for (i in 1:length(pos.start.doxygen.comment)
} # end header file loop --- for(header.file.name in header.file.names)   

# restore working directory
setwd(cur.dir)

save(functions, file = "functions.sdd")
