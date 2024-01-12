#!/bin/bash
FILES=$(ls ~/p-sy58-0/remapping/fasta_files/)
for F in $FILES; do
  ID=${F%%_R1_phiXremoved.fq.gz*}
  if [ "${#ID}" -lt "${#F}" ]; then
    echo $ID
    qsub -v IDARG=$ID mapping.pbs
    sleep 5s
  fi
done
