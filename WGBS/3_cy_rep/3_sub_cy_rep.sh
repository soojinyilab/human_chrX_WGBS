#!/bin/bash
FILES=$(ls ~/p-sy58-0/remapping/meth_extract/human_XY/)
for F in $FILES; do
  ID=$(echo $F | sed 's/_processed_bismark_bt2_pe.deduplicated.bedgraph.gz.bismark.cov.gz//')
  if [ "${#ID}" -lt "${#F}" ]; then
    echo "$ID"
    qsub -v IDARG=$ID 3_make_cy_rep.pbs
    sleep 5s
  fi
done
