#!/usr/bin/env python


# python script to extract genome region of interest with inserted SNPs or indels from a vcf file

# NB: the genome sequence must be a linear string

import sys
import os
import glob
import pandas as pd
import numpy as np
import random
import re
import time


# sys.argv[1]: genome of interest
# sys.argv[2]: vcf file with SNPs and InDels


finish_counter = 0

lines_of_comments = 0

with open(sys.argv[2]) as f:
    for lines in f:
        if lines[:1] == '#':
            lines_of_comments += 1

print(lines_of_comments)

gene_variants = list()
new = list()

# genome of interest is opened and region of genome with gene of interest extracted

with open(sys.argv[1]) as f:
    for line in f:
        line = line.strip()
        if line[:1] == ">":

            chrom = re.search('^>([A-Za-z0-9\.]+) S', line).group(1)
            print(chrom)
            gene_variants.append('>' + str(chrom) + ' Speruvianum')

            os.system('grep \'' + chrom + '\' ' + sys.argv[2] + ' > temp.tsv')
            snps = pd.read_csv('temp.tsv', sep='\t', engine='python', skiprows = 1, names=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QsdUAL', 'FILTER', 'INFO', 'FORMAT', 'la1954_aln_pe_sorted.bam'])
            print(len(snps))


        else:   
            nb_snps = 0
            nb_counting_snps = 0
            indel_correction = 0
            counter = 0


# adjusting the pointer in the vcf file to the right chromosome

            if len(snps) == 0:
                gene_variants = gene_variants[:-1]
                finish_counter += 1

            else:

                region = line
                print(len(region))
                
                for j in range(len(snps)):

                    counter += 1

                    if counter%100000 == 0:
                        print(counter*100/len(snps), '%\t', time.process_time()/60, 'min.')


                    mq = re.search('MQ=([0-9\.]*)', snps['INFO'][j])
                    if mq != None:
                        if int(mq.group(1)) >= 40:

                            ref = str(snps['REF'][j])
                            len_ref = len(str(ref))
                            before = snps['POS'][j] - 1 + indel_correction
                            after = before + len(str(ref))


                            alt = str(snps['ALT'][j]).split(",")

                            if len(alt) > 1:
                                alt = random.choice(alt)

                            else:
                                alt = alt[0]

                            len_alt = len(alt)
                            
                            if region[before:after].lower() == ref.lower():

                                region = region[:before] + alt + region[before+len_ref:]
                                indel_correction = indel_correction + len_alt - len_ref
                                nb_counting_snps += 1
                                nb_snps += 1

                            else:
                                print('ERROR')


                       
                    if j == len(snps) - 1:
                        gene_variants.append(region)
                        print(len(region))
                        print('time for chromosome: ' + str(time.process_time()))
                        finish_counter += 1


thefile = open("Speruvianum_reference_Sl30_NCBI.fasta", 'w')
for item in gene_variants:
    thefile.write("%s\n" % item)
thefile.close()
print(time.process_time())
