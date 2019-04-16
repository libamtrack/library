#!/bin/bash
cd $GSL_PATH/gsl-latest/gsl && emconfigure ./configure --prefix=$PWD/usr && emmake make && emmake make install && cd ../../

cp $GSL_ROOT_DIR/lib/libgsl.a .
cp $GSL_ROOT_DIR/lib/libgslcblas.a .
