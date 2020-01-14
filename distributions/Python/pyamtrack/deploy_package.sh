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

# make a source package
pip3 install -U twine --user

# upload only if tag present
if [[ $TRAVIS_TAG != "" ]]; then
  twine upload $TRAVIS_BUILD_DIR/distributions/Python/pyamtrack/generated/dist/wheelhouse/*.whl
fi