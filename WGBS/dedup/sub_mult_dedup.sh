#!/bin/bash
FILES=$(ls ~/p-sy58-0/remapping/bs_mapped/human_XY/)
for F in $FILES; do
 ID=${F%%_processed_bismark_bt2_pe.bam*}
 if [ "${#ID}" -lt "${#F}" ]; then
   echo "$ID"
   qsub -v IDARG=$ID dedup.pbs
   sleep 5s
  fi
done
