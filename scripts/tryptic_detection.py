#!/usr/bin/env python

import numpy as np
import pandas as pd
import sys

df1 = pd.read_csv(sys.argv[1], sep='\t', engine='python')
bion_cols = ['Sequence']
[bion_cols.append("LFQ intensity {}".format(x)) for x in range(16, 20)]

new = []

# Split tryptic peptides file into data frame and only keep peptides > 6 amino acids

with open(sys.argv[2]) as f:
    for lines in f:
        lines = lines.strip()
        if lines[:1] == '>':
            name = lines[1:]
            name = name.split(' ')[0]

        elif len(lines) > 6:
            new.append([lines, name])

df2 = pd.DataFrame(new, columns=['peptide', 'prot_id'])

# Go through each protein of interest and determine number of peptides detected in at least 3 Bion-treated samples and match tryptic peptides with predicted peptides

print("prot_id\thyp_trypt_pept\tmatched_trypt_pept\ttotal_trypt_pept")

for i in np.unique(df2.prot_id):
    matches = []
    peptides_found = []
    keep_rows = []
    
    for row in df1.index:
        prots = df1.loc[row, 'Proteins'].split(';')
        if i in prots:
            keep_rows.append(row)

    tmp1 = df1[df1.index.isin(keep_rows)][bion_cols]
    tmp2 = df2[df2.prot_id == i]

    seq_trypt = tmp2.peptide.to_numpy()


    for row in tmp1.index:
        hit = 0

        for col in tmp1.columns:
            if col != 'Sequence':
                if tmp1.loc[row, col] != 0:
                    hit += 1


        if hit >= 3:
            seq_found = tmp1.loc[row, 'Sequence']

            for tryptic in seq_trypt:
                if seq_found in tryptic:
                    matches.append(tryptic)
                    peptides_found.append(seq_found)
                
                elif tryptic in seq_found:
                    matches.append(tryptic)
                    peptides_found.append(seq_found)



    print("{}:\t{}\t{}\t{}".format(i, len(seq_trypt), len(np.unique(matches)), len(np.unique(peptides_found))))

df2.to_csv("./{}".format(sys.argv[3]), index=False, sep='\t')
