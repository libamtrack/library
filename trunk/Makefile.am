SUBDIRS=thirdparty/cernlib src test example
DIST_SUBDIRS=thirdparty/cernlib src test example
ACLOCAL_AMFLAGS =-I m4

R-package:
	$(srcdir)/distributions/R/create.package.sh

R-package-sandbox:
	$(srcdir)/distributions/R/create.package.sh --sandbox

R-forge:
	$(srcdir)/distributions/R/create.package.sh --sync /home/greilich/workspaces/libamtrack/R-Forge/pkg/

deb-package:
	echo $(VERSION)
	echo $(srcdir)

