#!/bin/bash
LETTERS=$(ls ~/p-sy58-0/remapping/dedup/human_XY/)
for L in $LETTERS; do
  ID=${L%%_processed_bismark_bt2_pe.deduplicated.bam*}
  if [ "${#ID}" -lt "${#L}" ]; then
    echo "$ID"
    qsub -v IDARG=$ID 1_meth_ex.pbs
    sleep 5s
  fi
done
