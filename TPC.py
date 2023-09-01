# -*- coding: utf-8 -*-
"""
Created on Tue May 10 15:25:44 2022

@author: zjy52
"""

import sys
import pandas as pd
from Bio import SeqIO
import numpy as np
File=sys.argv[1]
Out=sys.argv[2]
aa=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
dp=[]
for i in aa:
    for j in aa:
        for l in aa:
            dp.append(i+j+l)
name=[]
se=[]
for seq in SeqIO.parse(File,'fasta'):
    name.append(seq.id)
    se.append(str(seq.seq))
a=np.linspace(0,len(name)-1,len(name))
co=['id']+dp
data=pd.DataFrame(index=a,columns=co)
data['id']=name

for i in range(len(data)):
    d=[]
    for j in range(len(se[i])-2):
        d.append(se[i][j]+se[i][j+1]+se[i][j+2])
    for m in dp:
        n=0
        for o in range(len(d)):
            if d[o]==m:
                n=n+1
        p=(n/(len(se[i])-2))*100
        data.iloc[i,:][m]=p
data=data.set_index('id')
data.to_csv(Out+'/'+File.split('.')[0]+'_TPC.csv')
