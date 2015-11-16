#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
mkdir /tmp/bayestraq
pushd /tmp/bayestraq
wget http://mirrors.ebi.ac.uk/CRAN/src/base/R-latest.tar.gz	
tar xzf R-latest.tar.gz
pushd R-*
./configure --with-x=no --with-readline=no --prefix=/tmp/bayestraq/R
make
make install
popd
R/bin/Rscript $DIR/bayestraq_Rsetup.R
rm $DIR/bayestraq_R.zip
zip -r $DIR/bayestraq_R.zip R 
popd
