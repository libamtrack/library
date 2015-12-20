# Script to collect and add metainfo (version, date, svn revision, etc.)
# to R package during compilation
#
# Created: S. Greilich, 2011-09-04
##rev##

# Arg 1: Root source dir with version information in configure.ac / saved_svn_version.txt
# Arg 2: package name
# Arg 3: work dir
rm(list = ls())

#pass the name of the library from shell
args <- commandArgs(TRUE)

try(configure.ac.file <-scan(paste(args[1], "configure.ac", sep = "/"), what = character(), sep = "\n"))

date.long   <- Sys.time() # format(Sys.time(), "%a %b %d %H:%M:%S %Y")
date.short  <- Sys.Date()

get.argument <- function(x, arg.no = 2, filter = "[[:punct:]]|^[[:blank:]]"){
  return(gsub(filter, "", strsplit(x, ",")[[1]][arg.no]))
}

# Get version from AC_DEFINE
ii          <- grepl("AC_INIT", configure.ac.file)
code.version<- get.argument(configure.ac.file[ii], filter = "[[]|[]]|^[[:blank:]]")

# Get name, status from AC_DEFINE
ii          <- grepl("AC_DEFINE", configure.ac.file)
jj          <- grepl("CODE_NAME", configure.ac.file)
kk          <- grepl("CODE_STATUS", configure.ac.file)
code.name   <- get.argument(configure.ac.file[ii&jj]) 
code.status <- get.argument(configure.ac.file[ii&kk]) 

# Write to transient file to use in other scripts
meta.information <- list( code.version = code.version, 
            code.name    = code.name, 
            code.status  = code.status, 
            date.short   = date.short, 
            date.long    = date.long)
save( meta.information,
      file = "meta.information.sdd")

# Insert information
desc.file   <- scan(paste(args[3], args[2], "DESCRIPTION", sep = "/"), what = character(), sep = "\n")
desc.file   <- gsub("##VERSION##", code.version, desc.file)
desc.file   <- gsub("##DATE##", date.short, desc.file)
write(desc.file, file = paste(args[3], args[2], "DESCRIPTION", sep = "/"))

init.file   <- scan(paste(args[3], args[2], "R/initial.R", sep = "/"), what = character(), sep = "\n")
init.file   <- gsub("##VERSION##", code.version, init.file)
init.file   <- gsub("##NAME##", code.name, init.file)
init.file   <- gsub("##STATUS##, ##REVISION##, ##DATE##", date.short, init.file)
write(init.file, file = paste(args[3], args[2], "R/initial.R", sep = "/"))

docu.file   <- scan(paste(args[3], "/", args[2], "/man/", args[2], "-package.Rd", sep = ""), what = character(), sep = "\n")
docu.file   <- gsub("##VERSION##", code.version, docu.file)
docu.file   <- gsub("##DATE##", date.short, docu.file)
docu.file   <- gsub("##NAME##", code.name, docu.file)
write(docu.file, file = paste(args[3], "/", args[2], "/man/", args[2], "-package.Rd", sep = ""))
