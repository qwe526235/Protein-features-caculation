# -*- coding: utf-8 -*-
"""
Created on Tue May 10 11:42:14 2022

@author: zjy52
"""


import pandas as pd
from Bio import SeqIO
import numpy as np
from collections import Counter
import sys
File=sys.argv[1]
File2=sys.argv[2]
Out=sys.argv[3]
aa=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
ex=['X','B','Z','J','O','U']
ax=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','X']
co=['id']
for i in ax:
    co.append('QSO_'+i)
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
l=list(da.columns.values)
st=pd.DataFrame(index=['std'],columns=l)
for j in l:
    n=[]
    for i in range(len(da.loc['A':'Y',:])):
       n.append(abs(da.loc[:,j][i]-da.loc[:,j][i+1]))
    st.loc['std',j]=np.std(n)
        
##distance #d=10 w=0.1
dt=pd.DataFrame(index=ax,columns=ax)
for i in ax:
    for j in ax:
        dt.loc[i,j]=((abs(da.loc[i,'Polarity']-da.loc[j,'Polarity'])/st.loc['std','Polarity'])*2+(abs(da.loc[i,'Volume']-da.loc[j,'Volume'])/st.loc['std','Volume'])*2)*0.05     
data=pd.DataFrame(index=a,columns=co)
data['id']=name
for i in range(len(data)):
    s=list(se[i])
    for j in range(len(s)):
        for m in ex:
            if s[j]==m:
                s[j]='X'
    s=''.join(s)
    n=0
    for k in range(11):
        for l in range(len(s)-k):
            n=n+(dt.loc[s[l],s[l+k]])**2
    for o in ax:
        data.iloc[i,:]['QSO_'+o]=Counter(s)[o]/(1+0.1*n)
        
data=data.set_index('id')
data.to_csv(Out+'/'+File.split('.')[0]+'_Quasi-sequence order.csv')