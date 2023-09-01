# -*- coding: utf-8 -*-
"""
Created on Mon May  9 11:48:12 2022

@author: zjy52
"""
import sys
import pandas as pd
from Bio import SeqIO
import numpy as np
File=sys.argv[1]
File2=sys.argv[2]
Out=sys.argv[3]
aa=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
df=pd.read_csv(File2)
da=df
for i in range(len(df)):
    for j in range(1,9):
        da.iloc[i,j]=(df.iloc[i,j]-np.mean(df.iloc[0:19,j]))/np.std(df.iloc[0:19,j])
name=[]
se=[]
for seq in SeqIO.parse(File,'fasta'):
    name.append(seq.id)
    se.append(str(seq.seq))
for i in se:
    i=i.upper()
a=np.linspace(0,len(name)-1,len(name))
co=['Hydrophobicity','Hydrophilicity','Side Mass','pK1(a-CO2H)','pK2(NH3)','pI','Polarity','Volume']
co=['id']+co
data=pd.DataFrame(index=a,columns=co)
data['id']=name
for i in range(len(data)):
    n=0
    l=0
    k=0
    m=0
    p=0
    u=0
    s=0
    v=0
    for d in range(len(da)):
        for j in range(len(se[i])):
            if se[i][j]==df['Amino acid'][d]:
                n=n+da['Hydrophobicity'][d]
                l=l+da['Hydrophilicity'][d]
                k=k+da['Side Mass'][d]
                m=m+da['pK1(a-CO2H)'][d]
                p=p+da['pK2(NH3)'][d]
                u=u+da['pI'][d]
                s=s+da['Polarity'][d]
                v=v+da['Volume'][d]
    data.iloc[i,1]=n
    data.iloc[i,2]=l
    data.iloc[i,3]=k
    data.iloc[i,4]=m
    data.iloc[i,5]=p
    data.iloc[i,6]=u
    data.iloc[i,7]=s
    data.iloc[i,8]=v
    
data=data.set_index('id')
data.to_csv(Out+'/'+File.split('.')[0]+'_total_properties.csv')   