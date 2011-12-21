# Script to collect and add metainfo (version, date, svn revision, etc.)
# to R package during compilation
#
# Created: S. Greilich, 2011-09-04
##rev##

rm(list = ls())

configure.ac.file <-scan("../../../configure.ac", what = character(), sep = "\n")

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

# Get svn revision from saved file (which works whether svnversion runs or not)
code.rev    <- scan("../../../saved_svn_version.txt", what = character())

# Write to transient file to use in other scripts
meta.information <- list( code.version = code.version, 
            code.name    = code.name, 
            code.status  = code.status, 
            date.short   = date.short, 
            date.long    = date.long)
save( meta.information,
      file = "meta.information.sdd")

# Insert information
desc.file   <- scan("./libamtrack/DESCRIPTION", what = character(), sep = "\n")
desc.file   <- gsub("##VERSION##", code.version, desc.file)
desc.file   <- gsub("##DATE##", date.short, desc.file)
write(desc.file, file = "./libamtrack/DESCRIPTION")

init.file   <- scan("./libamtrack/R/initial.R", what = character(), sep = "\n")
init.file   <- gsub("##VERSION##", code.version, init.file)
init.file   <- gsub("##NAME##", code.name, init.file)
if (code.status == "Release"){
   init.file   <- gsub("##STATUS##, ##REVISION##, ##DATE##", date.short, init.file)
}else{
   init.file   <- gsub("##STATUS##", code.status, init.file)
   init.file   <- gsub("##REVISION##", code.rev, init.file)
   init.file   <- gsub("##DATE##", date.long, init.file)
}
write(init.file, file = "./libamtrack/R/initial.R")

docu.file   <- scan("./libamtrack/man/libamtrack-package.Rd", what = character(), sep = "\n")
docu.file   <- gsub("##VERSION##", code.version, docu.file)
docu.file   <- gsub("##DATE##", date.short, docu.file)
docu.file   <- gsub("##NAME##", code.name, docu.file)
write(docu.file, file = "./libamtrack/man/libamtrack-package.Rd")
