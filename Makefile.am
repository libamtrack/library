SUBDIRS=src test example
DIST_SUBDIRS=src test example
ACLOCAL_AMFLAGS =-I m4

R-package:
	$(srcdir)/distributions/R/create.package.sh

R-forge:
	$(srcdir)/distributions/R/create.package.sh /home/greilich/workspaces/libamtrack/R-Forge/pkg/

deb-package:
	echo $(VERSION)
	echo $(srcdir)

