#!/bin/bash

####Parameters
input=$1
output=$2

####Main
awk -F "\t" 'BEGIN {FS="\t";OFS="\t"}{if($4!="na"){print $4,$5,$6,$2}}' "$input" > "$output"