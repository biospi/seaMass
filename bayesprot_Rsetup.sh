#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
pushd $DIR
export R_ROOT_DIR=$DIR/bayesprot_R
export RHOME=${R_ROOT_DIR}
mkdir bayesprot_R
wget http://mirrors.ebi.ac.uk/CRAN/src/base/R-3/R-3.2.5.tar.gz
tar xzf R-3.2.5.tar.gz
rm R-3.2.5.tar.gz
pushd R-*
./configure --with-x=no --with-readline=no --prefix=$DIR/bayesprot_R
make
make install
popd
pushd bayesprot_R
sed 's;'"${R_ROOT_DIR}"';\${R_ROOT_DIR};g' < bin/R > R.tmp
chmod 755 R.tmp
cp R.tmp bin/R
mv R.tmp lib*/R/bin/R
bin/Rscript $DIR/bayesprot_Rsetup.R
popd
zip -rq bayesprot_R.zip bayesprot_R
popd
