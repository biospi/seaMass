#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
unzip -n -qq bayestraq_R.zip
rm bayestraq_R.zip
export R_ROOT_DIR=$DIR/bayestraq_R
export RHOME=${R_ROOT_DIR}
bayestraq_R/bin/Rscript "$@"
