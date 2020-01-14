#!/bin/bash

# compile GSL sources
cd /io/gsl-2.6
./configure --prefix=/io/gsl
make -j2
make install

for PYVER in cp36-cp36m cp37-cp37m; do

  # install cBinder package from a local folder
  "/opt/python/${PYVER}/bin/python" -m pip install /io/cBinder/

  # generate extension and make wheel package
  rm -rf /io/generated
  CFLAGS='-std=c99' "/opt/python/${PYVER}/bin/python" -m cBinder pyamtrack \
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
  cd /io/generated
  rm -rf build
  rm -rf dist
  /opt/python/$PYVER/bin/python3 setup.py bdist_wheel

  # TODO check code below
  TAG=`echo ${PYVER} | cut -d'-' -f1`
  rename  's/-none-/${TAG}/' *whl

done

# add manylinux1 tag
#find /io/generated/dist -name "*.whl" -exec auditwheel repair {} \;
