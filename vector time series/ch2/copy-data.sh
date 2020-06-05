#!/bin/bash

cd $(dirname $0)

cat ../../data/WW2a.csv | \
sed '1 s/.*/Period,1Automobile,2Building Materials,3Clothing,4Beer wine and liquor,5furnitures,6General Merchandise,7Grocery,8household Appliances/' > WW2a.csv
