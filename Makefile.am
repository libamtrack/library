SUBDIRS=thirdparty/cernlib src test example
DIST_SUBDIRS=thirdparty/cernlib src test example
ACLOCAL_AMFLAGS =-I m4

R-package:
	$(srcdir)/distributions/R/create.package.sh $(RTEST_FLAG)
	
deb-package:
	echo $(VERSION)
	echo $(srcdir)

