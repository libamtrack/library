# Script to create canonical namespace file for R package
# Mandatory from R 2.14.0
#
# Created: S. Greilich, 2011-11-16
##rev##

# Arguments
# 1: work dir
# 2: package name
# 3: package acronym
rm(list = ls())

args <- commandArgs(TRUE)

namespace.file          <- paste("# Namespace for package ", args[2], "\n#Created by R.create.package.namespace.R script\n\n", sep = "")

# Add imports
#namespace.file          <- paste( namespace.file,
#                                  "import(lattice)\n\n",
#                                  sep = "")

# Add loading hook for library
namespace.file          <- paste( namespace.file,
                                  "useDynLib(\"",args[2], "\", .registration = TRUE)\n\n",
                                  sep = "")

# Read functions from documentation directory. This is easiest as all functions have their corresponding
# Rd file (alternatively one could browse C code NAMESPACE and hardcoded wrappers).
# Remove Rd files for variables and other entities (which does not start with "AT.")
Code.namespace         <- list.files(path = paste(args[1], "/", args[2], "/man", sep = ""), 
                                     pattern = "*.Rd")
Code.namespace         <- gsub(".Rd", "", Code.namespace)
library.acronym        <- paste("^", args[3], ".", sep = "")
Code.namespace         <- Code.namespace[grepl(library.acronym, Code.namespace)]

namespace.file          <- paste( namespace.file,
                                  "export(\n",
                                  sep = "")

for (i in 1:length(Code.namespace)){
    # i <- 1
    delim <- ","
    if(i == length(Code.namespace)){
        delim <- ""
    }  
    namespace.file          <- paste( namespace.file,
                                      "\t\"",
                                      Code.namespace[i],
                                      "\"",
                                      delim,
                                      "\n",
                                      sep = "")
}

namespace.file          <- paste( namespace.file,
                                  ")\n\n",
                                  sep = "")

write(namespace.file, file = paste(args[1], "/", args[2], "/NAMESPACE", sep = ""))

