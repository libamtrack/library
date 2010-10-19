Installation instructions for downloaded R packages

== Linux ==

* For local installation without root privileges (e.g. as user 'parryhotter' in a folder 'home/parryhotter/R'):
** In the console (not within R), run from the folder you downloaded the package into: 'R CMD INSTALL ./libamtrack.tar.gz -l home/parryhotter/R'
** To use libamtrack then in an R session type 'library(libamtrack, lib.loc="home/parryhotter/R")'
** To uninstall the package again run (in the console again): 'R CMD REMOVE libamtrack -l home/parryhotter/R'

* For installation using root privileges:
** From the folder you downloaded the package into, run (in the console) 'R CMD INSTALL ./libamtrack.tar.gz'
** This will install libamtrack into a system folder so any user can deploy it using 'library(libamtrack,)' during an R session.
** To uninstall the package again run (in the console) 'R CMD REMOVE libamtrack '


== Windows ==

* Start R (i.e. the R GUI). In the menu bar you will find “Packages” and “Install packages from local zip-file”. Navigate to the directory you downloaded the file to and select it.
* After successful installation, use 'library(libamtrack)' to start working. 

== General ==

* “AT.” and a double click on “tab” will give you a list of available functions
