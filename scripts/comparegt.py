#!/usr/bin/env python

import pandas as pd
import numpy as np
import sys

chrom = sys.argv[3]
norm = int(sys.argv[4])

nanovcf = pd.read_csv(sys.argv[1], sep='\t', comment="#", 
                      names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']).query(
                          "REF.str.len() == 1 & ALT.str.len() == 1 & FILTER == 'PASS'")[['CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'SAMPLE']]
nanovcf['SAMPLE'] = nanovcf['SAMPLE'].apply(lambda x: list(set(x.split(":")[0].replace("/", "|").split("|"))))
nanovcf['GT'] = nanovcf.apply(lambda df: [df['REF'] if x == "0" else df['ALT'] for x in df['SAMPLE']], axis=1)


mpileup = pd.read_csv(sys.argv[2], sep='\t', comment="#", 
                      usecols=[0,1,4], names=['CHROM', 'POS', 'NUC']).query("CHROM == @chrom")
mpileup['NUC'] = mpileup['NUC'].apply(lambda x: [ s for s in x.split(",") if s != "<*>"])


merged = pd.merge(nanovcf, mpileup, on=['CHROM', 'POS'], how='inner')
merged['concor'] = merged.apply(lambda df: all(item in df['GT'] for item in df['NUC']), axis=1)

if ( norm > 0 ):
    merged = merged.sample(norm, random_state=330)

print(merged['concor'].shape[0], 
      merged['concor'].value_counts().get(False), 
      merged['concor'].value_counts(True).get(False,0).round(5))