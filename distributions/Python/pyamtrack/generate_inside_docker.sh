#!/bin/bash

# compile GSL sources
cd /io/gsl-latest || exit 1
./configure --prefix=/io/gsl || exit 1
#  -j [N], --jobs[=N]          Allow N jobs at once; infinite jobs with no arg.
make -j || exit 1
make install || exit 1

# install cBinder package from a local folder
/opt/python/cp36-cp36m/bin/python3 -m pip install /io/cBinder/ || exit 1

# generate extension and make wheel package
rm -rf /io/generated || exit 1
mkdir /io/generated || exit 1
cd /io/generated || exit 1
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
-l gsl -l gslcblas -l m || exit 1

# fix git tags and create clean repo
cd /io/libamtrack || exit 1
ls -alh . || exit 1
git describe --tags --abbrev=0 || exit 1
TAG=`git describe --tags --abbrev=0`
rm -rf .git || exit 1
git init || exit 1
git add * || exit 1
git config --global user.email "leszek.grzanka@ifj.edu.pl" || exit 1
git config --global user.name "Leszek Grzanka" || exit 1
git commit -m "Clean commit" || exit 1
git tag -a $TAG -m "Tag" || exit 1

# generate package once more with adjusted setup.py
# pure `cp` ask for confirmation to overwrite file even with `-f` option
# this may be due to the fact that `cp` was aliased to `cp -i` (interactive)
# therefore we avoid this problem by using pure `/bin/cp`
/bin/cp --force /io/setup.py /io/generated || exit 1
cd /io/generated || exit 1
rm -rf build || exit 1
rm -rf dist || exit 1
/opt/python/cp36-cp36m/bin/python3 setup.py bdist_wheel || exit 1

# add manylinux1 tag
cd /io/generated/dist || exit 1
auditwheel repair *whl || exit 1
