import numpy as np
import pandas as pd

# Created from VCF with GATK VariantsToTable
snp_file = '~/p-sy58-0/1000G_phase1.snps.high_confidence.hg38.tsv'
dbsnp = pd.read_csv(snp_file, sep='\t', dtype={'AF': np.float})
dbsnp.dropna(inplace=True)

af_cutoff = 0.01

dbsnp = dbsnp[(dbsnp['CHROM'] == 'chr5') & (dbsnp['AF'] >= af_cutoff)]
dbsnp = dbsnp[((dbsnp['REF'] == 'C') & (dbsnp['ALT'] == 'T')) | ((dbsnp['REF'] == 'G') & (dbsnp['ALT'] == 'A'))]
dbsnp.to_csv(f'~/p-sy58-0/1000G_phase1.chr5_minAF_{af_cutoff}_BS_snps.hg38.tsv', sep='\t', index=False)

