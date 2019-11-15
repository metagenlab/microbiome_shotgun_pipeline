input=$1
output=$2
awk 'BEGIN {FS="\t";OFS="\t"}{if($4!="na" && $5!="na"){ print $4,$5,$6,$2 }}' "$input" > "$output"