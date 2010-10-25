.First.lib <- function(lib, pkg){
 library.dynam("libamtrack", pkg, lib)
 cat("This is libamtrack 0.3 (2010-10-22). Type '?libamtrack' for help.\n")
}

.Last.lib <- function(libpath){
 try(library.dynam.unload("libamtrack", libpath))
}

