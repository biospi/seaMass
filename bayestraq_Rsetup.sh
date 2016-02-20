#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
pushd $DIR
export R_ROOT_DIR=$DIR/bayestraq_R
export RHOME=${R_ROOT_DIR}
mkdir bayestraq_R
wget http://mirrors.ebi.ac.uk/CRAN/src/base/R-latest.tar.gz
tar xzf R-latest.tar.gz
rm R-latest.tar.gz
pushd R-*
./configure --with-x=no --with-readline=no --prefix=$DIR/bayestraq_R
make
make install
popd
pushd bayestraq_R
sed 's;'"${R_ROOT_DIR}"';\${R_ROOT_DIR};g' < bin/R > R.tmp
chmod 755 R.tmp
cp R.tmp bin/R
mv R.tmp lib*/R/bin/R
bin/Rscript $DIR/bayestraq_Rsetup.R
popd
zip -rq bayestraq_R.zip bayestraq_R
popd
