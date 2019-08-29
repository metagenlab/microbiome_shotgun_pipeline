#!/bin/bash
input_file=$1
echo "input file: $input_file"
covsum=$(awk '{s+=$3}END{print s}' "$input_file")
if [ "$covsum" == 0 ]
then
  echo "$covsum coverage found, removing file"
  rm "$input_file"
else
  echo "found $covsum bases coverage, cat into $2"
  cat "$input_file" >> "$2"
fi