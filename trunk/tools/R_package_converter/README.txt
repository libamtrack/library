In order to run package converter you have to:

1. Set the R working directory to this folder

2. Run "read.header.R" script

3. Run "convert.header.to.Rwrapper.R" script

4. Run "convert.header.to.cwrapper.R" script

5. Run "convert.doxygen.to.Rd.R" script

6. all required *.c and *.h files have to be copied in to ./package/src
	and optionally a Makevars ore Makefile to compiling and linking into 
	a shared object. (link to documentation: http://cran.r-project.org/doc/manuals/R-exts.pdf)

7. Test the makefile with "R CMD SHLIB" by creating a shared library and load it into R
	
8. Run the command "R CMD check package" to check for possible errors

9. If all errors are fixed, use "R CMD build package" to build package