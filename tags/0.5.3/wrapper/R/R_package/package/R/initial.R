.onLoad <- function(lib, pkg){
 library.dynam("libamtrack", pkg, lib)
 packageStartupMessage("This is libamtrack ##VERSION## '##NAME##' (##STATUS##, ##REVISION##, ##DATE##).\nType '?libamtrack' for help.\n")
}

.Last.lib <- function(libpath){
 try(library.dynam.unload("libamtrack", libpath))
}

