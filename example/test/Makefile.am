bin_PROGRAMS = amtrack_test
amtrack_test_SOURCES = AT_test.c
amtrack_test_CFLAGS =-I../../include -I../../thirdparty/cernlib -Wall -std=gnu99 -fPIC $(GSL_CFLAGS)
amtrack_test_LDADD = ../../src/libamtrack.la  $(GSL_LIBS)
