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
# Started 2010-12-01, Steffen Greilich
#################################################

# Navigate to wrapper directory
start.dir           <- getwd()

setdir_win <- try(setwd("..\\..\\wrapper\\R\\R_direct_access"))
setdir_lin <- try(setwd("../../wrapper/R/R_direct_access"))

if(setdir_win == FALSE && setdir_lin == FALSE){
   stop("Please start script from /test/R directory")    
}

# Trigger new build
try(system("create.direct.access.windows.bat"))
try(system("bash create.direct.access.linux.sh"))

# Load library and wrappers
try(dyn.load("libamtrack.dll"))
try(dyn.load("libamtrack.so"))
try(dyn.load("libamtrack.dylib"))

# Source wrappers
source("libamtrack.R")

# Return to original directory
setwd(start.dir)
