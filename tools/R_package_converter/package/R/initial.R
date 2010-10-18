.First.lib <- function(lib, pkg){
 library.dynam("libamtrack", pkg, lib)
 print("This is libamtrack 0.2 (2010-10-18).")
}

.Last.lib <- function(libpath){
 try(library.dynam.unload("libamtrack", libpath))
}

