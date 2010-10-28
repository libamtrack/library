.First.lib <- function(lib, pkg){
 library.dynam("libamtrack", pkg, lib)
 cat("This is libamtrack 0.4.0 (2010-10-29). Type '?libamtrack' for help.\n")
}

.Last.lib <- function(libpath){
 try(library.dynam.unload("libamtrack", libpath))
}

