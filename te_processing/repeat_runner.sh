#!/usr/bin/env bash

file_inputs=$1
file_pulled_from=$2

while IFS= read -r line
do
  IFS='_' read -r -a array <<< "$line"
  if [ ! -d "${array[0]}" ];
    then
    mkdir ${array[0]}
  fi
  if [ ! -d "${array[0]}/${array[1]}" ];
    then
    mkdir "${array[0]}/${array[1]}"
  fi
  mkdir "${array[0]}/${array[1]}/${array[2]}"
  gunzip "$file_pulled_from/$line" "${array[0]}/${array[1]}/${array[2]}"
  fasta=${line%".gz"}
  db=${line%".fa"}
  cd "${array[0]}/${array[1]}/${array[2]}" || exit
  bsub -o "db.out" -e "db.err" -M200 -R 'select[mem>200] rusage[mem=200]' -J make_$db  "BuildDatabase -name $db $fasta"
  bsub -q long -o "mod.out" -e "mod.err" -M10000 -R 'select[mem>10000] rusage[mem=10000] span[hosts=1]' -n 16 -J study_$db -w 'ended("make_$db")' "RepeatModeler -database $db -pa 4 -LTRStruct >& run.out"
  cd ../../..
done < "$file_inputs"
