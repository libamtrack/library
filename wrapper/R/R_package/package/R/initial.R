.First.lib <- function(lib, pkg){
 library.dynam("libamtrack", pkg, lib)
 cat("This is libamtrack 0.5.1 'Blue Wombat' (devel, r997, 2011-06-15).\nType '?libamtrack' for help.\n")
}

.Last.lib <- function(libpath){
 try(library.dynam.unload("libamtrack", libpath))
}

