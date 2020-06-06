#!/bin/bash

echo "Downloading the data appendix from wiley.com..."

cd $(dirname $0)
wget https://media.wiley.com/product_ancillary/53/11195028/DOWNLOAD/Data%20Appendix-Wei.zip -O appendix.zip
unzip appendix.zip
mv "Data Appendix-Wei" "data"
