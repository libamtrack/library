#!/bin/bash

echo "Start"

cd ../
rm -rf test
mkdir -p test
cd test
mkdir -p libamtrack-0.5.3
cd libamtrack-0.5.3
echo `pwd`
mkdir -p include
cp -r ../../../../../include/*.h include
mkdir -p src
cp -r ../../../../../src/Makefile.am src 
cp -r ../../../../../src/*.c src
mkdir -p example/basic_plots
cp -r ../../../../../example/Makefile.am example 
cp -r ../../../../../example/basic_plots/Makefile.am example/basic_plots 
cp -r ../../../../../example/basic_plots/*.c example/basic_plots
mkdir -p example/demo
cp -r ../../../../../example/demo/Makefile.am example/demo 
cp -r ../../../../../example/demo/*.c example/demo
mkdir -p test/C
cp -r ../../../../../test/C/Makefile.am test/C 
cp -r ../../../../../test/Makefile.am test 
cp -r ../../../../../test/C/*.c test/C

cp -r ../../../../../configure.ac .
cp -r ../../../../../Makefile.am .
cp -r ../../../../../NEWS .
cp -r ../../../../../README .
cp -r ../../../../../AUTHORS .
cp -r ../../../../../ChangeLog .

mkdir m4
autoreconf --force --install

cd ../
tar -zcf libamtrack-0.5.3.tar.gz libamtrack-0.5.3

cd libamtrack-0.5.3

dh_make -s -c gpl2 -p libamtrack -e Leszek.Grzanka@ifj.edu.pl -f ../libamtrack-0.5.3.tar.gz
cp -f ../../debian/control debian/control

dpkg-buildpackage -rfakeroot

echo "End"