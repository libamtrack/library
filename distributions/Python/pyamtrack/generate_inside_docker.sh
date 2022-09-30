#!/bin/bash

# compile GSL sources
cd /io/gsl-latest
./configure --prefix=/io/gsl
#  -j [N], --jobs[=N]          Allow N jobs at once; infinite jobs with no arg.
make -j
make install

# install cBinder package from a local folder
/opt/python/cp36-cp36m/bin/python3 -m pip install /io/cBinder/

# generate extension and make wheel package
rm -rf /io/generated
mkdir /io/generated
cd /io/generated
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

# fix git tags and create clean repo
cd /io/libamtrack
TAG=`git describe --tags`
rm -rf .git
git init
git add *
git config --global user.email "leszek.grzanka@ifj.edu.pl"
git config --global user.name "Leszek Grzanka"
git commit -m "Clean commit"
git tag -a $TAG -m "Tag"

# generate package once more with adjusted setup.py
# pure `cp` ask for confirmation to overwrite file even with `-f` option
# this may be due to the fact that `cp` was aliased to `cp -i` (interactive)
# therefore we avoid this problem by using pure `/bin/cp`
/bin/cp --force /io/setup.py /io/generated
cd /io/generated
rm -rf build
rm -rf dist
/opt/python/cp36-cp36m/bin/python3 setup.py bdist_wheel

# add manylinux1 tag
cd /io/generated/dist
auditwheel repair *whl
