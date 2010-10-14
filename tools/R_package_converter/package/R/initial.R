.First.lib <- function(lib, pkg){
 library.dynam("libamtrack", pkg, lib)
 print("This is libamtrack (2010-10-13).")
}

.Last.lib <- function(lib, pkg){
 library.dynam.unload("libamtrack", pkg, lib)
 pos <- match("package:lattice", search())
 if(! is.na(pos)){
  detach(pos = pos)
 }
}

