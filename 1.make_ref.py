# -*- coding: utf-8 -*-
"""
Created on Sat Jul 30 22:22:19 2022

@author: congt
"""

#GGGATATTTCTCG CAGATCAA GAAAAAAAGCTTGCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGCAAGCAAAAAAAC TTGATCTG TGAGAAATATTCTTA
ini_seq = 'GGGATATTTCTCGCAGATCAAGAAAAAAAGCTTGCNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGCAAGCAAAAAAACTTGATCTGTGAGAAATATTCTTA'
tri_nu_comb = []
for nu1 in ['A','T','G','C']:
    for nu2 in ['A','T','G','C']:
        for nu3 in ['A','T','G','C']:
            tri_nu_comb.append(nu1+nu2+nu3)
seq = ''
seq_list = []   
for i in range(1,7):
    for comb1 in tri_nu_comb:
        for comb2 in tri_nu_comb:
            seq = ini_seq[:i+12] + comb1 + ini_seq[i+15:21] + ini_seq[21:81] + ini_seq[81:87-i] + comb2 + ini_seq[90-i:]
            if seq not in seq_list:
                seq_list.append(seq)
                
var_list = []
shRNA_list = []

import pandas as pd
df = pd.DataFrame()
for i, seq in enumerate(seq_list):
    var_list.append('>Variant_'+str(i+1).zfill(5)) #zfill outputs 00001 instead of 1
    shRNA_list.append(seq[:35] + seq[67:])
df['Variant'] = var_list
df['shRNA_list'] = shRNA_list

f = open('TLR10-15-reference.fa','w+')
for i,var in enumerate(var_list):
  f.write(var+'\n')
  f.write(shRNA_list[i]+'\n')
f.close()
