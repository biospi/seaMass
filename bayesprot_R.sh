#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
unzip -n -qq bayesprot_R.zip
rm bayesprot_R.zip
export R_ROOT_DIR=$DIR/bayesprot_R
export RHOME=${R_ROOT_DIR}
bayesprot_R/bin/Rscript "$@"
