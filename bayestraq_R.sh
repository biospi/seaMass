#!/bin/bash
mkdir /tmp/bayestraq
unzip -n -qq bayestraq_R.zip -d /tmp/bayestraq
/tmp/bayestraq/R/bin/Rscript "$@"
rm -rf /tmp/bayestraq
