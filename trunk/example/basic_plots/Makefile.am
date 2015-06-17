bin_PROGRAMS = amtrack_plot
amtrack_plot_SOURCES = AT_plot.c
amtrack_plot_CFLAGS =-I../../include -I../../thirdparty/cernlib -Wall -std=gnu99 -fPIC $(GSL_CFLAGS)
amtrack_plot_LDADD = ../../src/libamtrack.la  $(GSL_LIBS)
