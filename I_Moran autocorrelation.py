# -*- coding: utf-8 -*-
"""
Created on Tue May 10 09:24:07 2022

@author: zjy52
"""


import pandas as pd
from Bio import SeqIO
import numpy as np
import sys
File=sys.argv[1]
File2=sys.argv[2]
Out=sys.argv[3]
aa=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
df=pd.read_csv(File2)
da=df.set_index('Amino acid')
for i in range(len(df)):
    for j in range(0,8):
        da.iloc[i,j]=(df.iloc[i,j+1]-np.mean(df.iloc[0:19,j+1]))/np.std(df.iloc[0:19,j+1])
name=[]
se=[]
for seq in SeqIO.parse(File,'fasta'):
    name.append(seq.id)
    se.append(str(seq.seq))
for i in range(len(se)):
    se[i]=se[i].upper()
a=np.linspace(0,len(name)-1,len(name))
co=['Hydrophobicity','Hydrophilicity','Side Mass','pK1(a-CO2H)','pK2(NH3)','pI','Polarity','Volume']
c=['id']
for j in range(1,10):
    for i in co:
        c.append('I_'+i+'_'+str(j))
data=pd.DataFrame(index=a,columns=c)
data['id']=name
##d=10
for j in c[1:len(c)]:
    d=int(j.split('_')[2])
    for i in range(len(data)):
        n=0
        t=0
        for m in range(len(se[i])-d):
           a1=da.loc[se[i][m],j.split('_')[1]]
           a2=da.loc[se[i][m+d],j.split('_')[1]]
           n=n+(a1-np.mean(da.loc['A':'Y',j.split('_')[1]]))*(a2-np.mean(da.loc['A':'Y',j.split('_')[1]]))
        for l in range(len(se[i])):
            a1=da.loc[se[i][m],j.split('_')[1]]
            t=t+(a1-np.mean(da.loc['A':'Y',j.split('_')[1]]))**2
        data.loc[:,j][i]=(n/(len(se[i])-d))/(t/len(se[i]))
data=data.set_index('id')        
data.to_csv(Out+'/'+File.split('.')[0]+'_I_Moran autocorrelation.csv')