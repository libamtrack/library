
echo $GSL_PATH && cd $GSL_PATH && wget "http://ftpmirror.gnu.org/gnu/gsl/gsl-latest.tar.gz"
mkdir $GSL_PATH/gsl-latest  && tar -xzf gsl-latest.tar.gz -C $GSL_PATH/gsl-latest
mv $GSL_PATH/gsl-latest/** $GSL_PATH/gsl-latest/gsl
rm -r $GSL_PATH/gsl-latest.tar.gz
