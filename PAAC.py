# -*- coding: utf-8 -*-
"""
Created on Mon May  9 13:37:01 2022

@author: zjy52
"""
import sys
import pandas as pd
from Bio import SeqIO
import numpy as np
from collections import Counter
File=sys.argv[1]
File2=sys.argv[2]
Out=sys.argv[3]
aa=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
pa=[]
ps=[]
for i in aa:
    pa.append('PAAC_'+i)
for i in range(1,6):
    ps.append('PAAC_factor_'+str(i))
df=pd.read_csv(File2)
for i in range(len(df)):
    for j in range(1,9):
        df.iloc[i,j]=(df.iloc[i,j]-np.mean(df.iloc[:,j]))/np.std(df.iloc[:,j])
da=df.set_index('Amino acid')
name=[]
se=[]
for seq in SeqIO.parse(File,'fasta'):
    name.append(seq.id)
    se.append(str(seq.seq))
for i in range(len(se)):
    se[i]=se[i].upper()
a=np.linspace(0,len(name)-1,len(name))
to=['sequence order-correlated factors 1',
    'sequence order-correlated factors 2',
    'sequence order-correlated factors 3',
    'sequence order-correlated factors 4',
    'sequence order-correlated factors 5']
co=['id']+to+pa+ps
data=pd.DataFrame(index=a,columns=co)
data['id']=name
for i in range(len(data)):
    for l in range(1,6):
        for j in range(len(se[i])-l-1):
            a1=se[i][j]
            a2=se[i][j+l+1]
            for k in range(0,8):
                p=np.sum((da.loc[a1,:][k]-da.loc[a2,:][k])**2/8)/(len(se[i])-l-1)
                data.iloc[i,l]=p
         
for i in range(len(data)):
    for j in range(len(pa)):
        data[pa[j]][i]=Counter(se[i])[aa[j]]/(1+0.1*np.sum(data.iloc[i,1:6]))
    for j in range(len(ps)):
        data[ps[j]][i]=data.loc[:,to[j]][i]/(1+0.1*np.sum(data.iloc[i,1:6]))
data=data.set_index('id')
data.to_csv(Out+'/'+File.split('.')[0]+'_PAAC.csv')