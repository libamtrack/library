export GSL_PATH=$PWD
echo $GSL_PATH && cd $GSL_PATH && wget "http://ftpmirror.gnu.org/gnu/gsl/gsl-latest.tar.gz"
mkdir $GSL_PATH/gsl-latest  && tar -xzf gsl-latest.tar.gz -C $GSL_PATH/gsl-latest
mv $GSL_PATH/gsl-latest/** $GSL_PATH/gsl-latest/gsl
rm -r $GSL_PATH/gsl-latest.tar.gz

export GSL_ROOT_DIR=$GSL_PATH/gsl-latest/gsl/usr
export GSL_INCLUDE_DIR=$GSL_ROOT_DIR/include
export GSL_LIBRARY=$GSL_ROOT_DIR/lib/libgsl.a
export GSL_CBLAS_LIBRARY=$GSL_ROOT_DIR/lib/libgslcblas.a