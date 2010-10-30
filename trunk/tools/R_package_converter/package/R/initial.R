.First.lib <- function(lib, pkg){
 library.dynam("libamtrack", pkg, lib)
 cat("This is libamtrack 0.4.1 'Blue Ocelot' (2010-10-31).\nType '?libamtrack' for help.\n")
}

.Last.lib <- function(libpath){
 try(library.dynam.unload("libamtrack", libpath))
}

