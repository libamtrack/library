.First.lib <- function(lib, pkg){
 library.dynam("libamtrack", pkg, lib)
 print("This is libamtrack 0.2 (2010-10-18).")
}

.Last.lib <- function(lib, pkg){
 library.dynam.unload("libamtrack", pkg, lib)
 pos <- match("package:lattice", search())
 if(! is.na(pos)){
  detach(pos = pos)
 }
}

