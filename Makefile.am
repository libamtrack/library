SUBDIRS=src test example
DIST_SUBDIRS=src test example
ACLOCAL_AMFLAGS =-I m4

R-package:
	cd distributions/R
	./create.package.sh
	cp *.tar.gz ../..
	cd ../..

deb-package:
	echo $(VERSION)
	echo $(srcdir)

