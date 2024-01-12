#!/bin/bash
FILES=$(ls ~/p-sy58-0/remapping/meth_extract/human_XY/)
for F in $FILES; do
  ID=$(echo $F | sed 's/CpG_context_//')
  if [ "${#ID}" -lt "${#F}" ]; then
    ID=$(echo $ID | sed 's/_processed_bismark_bt2_pe.deduplicated.txt.gz//')
    echo "$ID"
    qsub -v IDARG=$ID 2_convert_bedGraph.pbs
    sleep 5s
  fi
done
