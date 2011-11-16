# Script to create canonical namespace file for R package
# Mandatory from R 2.14.0
#
# Created: S. Greilich, 2011-11-16
##rev##

rm(list = ls())

namespace.file          <- "# Namespace for package libamtrack\n#Created by R.create.package.namespace.R script\n\n"

# Add loading hook for library
namespace.file          <- paste( namespace.file,
                                  "useDynLib(\"libamtrack\", .registration = TRUE)\n\n",
                                  sep = "")

# Read functions from documentation directory. This is easiest as all functions have their corresponding
# Rd file (alternatively one could browse C code NAMESPACE and hardcoded wrappers).
# Remove Rd files for variables and other entities (which does not start with "AT.")
Code.namespace         <- list.files(path = "./libamtrack/man", 
                                     pattern = "*.Rd")
Code.namespace         <- gsub(".Rd", "", Code.namespace)
Code.namespace         <- Code.namespace[grepl("^AT.", Code.namespace)]

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

write(namespace.file, file = "./libamtrack/NAMESPACE")

