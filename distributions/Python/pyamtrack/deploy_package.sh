#!/bin/bash

set -x # Print command traces before executing command

set -e # Exit immediately if a simple command exits with a non-zero status.

set -o pipefail # Return value of a pipeline as the value of the last command to
                # exit with a non-zero status, or zero if all commands in the
                # pipeline exit successfully.


write_pypirc() {
PYPIRC=~/.pypirc

if [ -e "${PYPIRC}" ]; then
    rm ${PYPIRC}
fi

touch ${PYPIRC}
cat <<pypirc >${PYPIRC}
[distutils]
index-servers =
    pypi
[pypi]
username: ${PYPIUSER}
password: ${PYPIPASS}
pypirc

if [ ! -e "${PYPIRC}" ]; then
    echo "ERROR: Unable to write file ~/.pypirc"
    exit 1
fi
}

# write .pypirc file with pypi repository credentials
set +x
write_pypirc
set -x

echo "User" $PYPIUSER
pip3 install -U setuptools wheel --user

# TODO travis runs now on Ubuntu 16.04 with Python 3.5
# newer version of twine require at least Python 3.6, hence we install some older version here
pip3 install -U "twine<2" --user
twine --version

# upload only if tag present
if [[ $TRAVIS_TAG ]]; then
  python3 -m twine upload $TRAVIS_BUILD_DIR/distributions/Python/pyamtrack/generated/dist/wheelhouse/*.whl
fi