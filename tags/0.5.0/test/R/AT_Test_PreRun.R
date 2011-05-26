#################################################
# Pre-run script for R test scripts in libamtrack
#
# This will trigger the build of the latest
# version of library/wrappers for direct R 
# access and load them if necessary
# 
# Should be called by "source" statement from
# corresponding test script
#
# If 'recompile' is set FALSE recompiling is
# skipped and the existing library is loaded.
#
# Started 2010-12-01, Steffen Greilich
#################################################

# Navigate to wrapper directory
start.dir           <- getwd()
setdir              <- try(setwd(file.path("../../wrapper/R/R_direct_access")))

if(setdir == FALSE){
   stop("Please start script from /test/R directory")    
}

# Trigger new build (default, even if variable is not set)
if(!exists("recompile")){
    recompile <- TRUE
}

if(recompile == TRUE){
    if(.Platform$OS.type == "windows"){
      system("create.direct.access.windows.bat")}
    if(.Platform$OS.type == "unix"){
      system("bash create.direct.access.linux.sh")}
}

# Load library and wrappers
if(.Platform$OS.type == "windows"){
     dyn.load("libamtrack.dll")}
if(.Platform$OS.type == "unix"){
     dyn.load("libamtrack.so")}

# Source wrappers
source("libamtrack.R")

# Return to original directory
setwd(start.dir)
