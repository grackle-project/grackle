#!/bin/bash

set -eo pipefail

if [ "$#" -ne 1 ]; then
  echo "USAGE: $0 <PROJECT_PATH>"
fi
PROJECT_PATH="$1"

# check that all libraries included in the wheel have a license entry
# and confirm that we didn't include any unnecessary licenses
python ${PROJECT_PATH}/scripts/wheels/check_packaged_dependency_licenses.py

# launch the tests
cd ${PROJECT_PATH} && pytest --color=yes
