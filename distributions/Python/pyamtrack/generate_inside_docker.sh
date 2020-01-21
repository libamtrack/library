#!/bin/bash

# compile GSL sources
cd /io/gsl-2.6
./configure --prefix=/io/gsl
make -j2
make install

# install cBinder package from a local folder
/opt/python/cp36-cp36m/bin/python3 -m pip install /io/cBinder/

# generate extension and make wheel package
rm -rf /io/generated
CFLAGS='-std=c99' /opt/python/cp36-cp36m/bin/python3 -m cBinder pyamtrack \
-f /io/libamtrack/src/ \
-f /io/libamtrack/include/ \
-d /io/generated \
-es /io/symbols_to_export.txt  \
-mono libAT \
compile \
-i /io/gsl/include/ \
-i /io/libamtrack/include/ \
-b /io/gsl/lib/ \
-l gsl -l gslcblas -l m

# generate package once more with adjusted setup.py
cp /io/setup.py /io/generated
cp /io/version.py /io/generated/pyamtrack/
cd /io/generated
rm -rf build
rm -rf dist
/opt/python/cp36-cp36m/bin/python3 setup.py bdist_wheel

# add manylinux1 tag
cd /io/generated/dist
auditwheel repair *whl
