# -*- coding: utf-8 -*-
"""
Created on Mon May  9 08:08:19 2022

@author: zjy52
"""
import os
import pandas as pd
from Bio import SeqIO
import numpy as np
import sys
File=sys.argv[1]
Out=sys.argv[2]
aa=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
bb=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
##composition_class
name=[]
se=[]
for seq in SeqIO.parse(File,'fasta'):
    name.append(seq.id)
    se.append(str(seq.seq))
a=np.linspace(0,len(name)-1,len(name))
h_class=['T_Hydrophobicity_Polar','T_Hydrophobicity_Neutral','T_Hydrophobicity_Hydrophobicity']
NWV_class=['T_NWV_Class1','T_NWV_Class2','T_NWV_Class3']
Polarity_class=['T_Polarity_Class1','T_Polarity_Class2','T_Polarity_Class3']
Polarizability_class=['T_Polarizability_Class1','T_Polarizability_Class2','T_Polarizability_Class3']
Charge_class=['T_Positive','T_Neutral','T_Negative']
Secondary_structure=['T_Helix','T_Strand','T_Coil']
Solvent_accessibility=['T_Buried','T_Exposed','T_Intermediate']

#dataframe
co=['id']+h_class+NWV_class+Polarity_class+Polarizability_class+Charge_class+Secondary_structure+Solvent_accessibility
data=pd.DataFrame(index=a,columns=co)
data['id']=name
##h_class
c0=['R','K','E','D','Q','N']
c1=['G','A','S','T','P','H','Y']
c2=['C','L','V','I','M','F','W']
dp1=[]##12 21
dp2=[]##13 31
dp3=[]##23 32
for i in c0:
    for j in c1:
        dp1.append(i+j)
for i in c0:
    for j in c2:
        dp2.append(i+j)
for i in c1:
    for j in c2:
        dp3.append(i+j)
for i in c1:
    for j in c0:
        dp1.append(i+j)
for i in c2:
    for j in c0:
        dp2.append(i+j)
for i in c2:
    for j in c1:
        dp3.append(i+j)
for i in range(len(data)):
    n=0
    d=[]
    for j in range(len(se[i])-1):
        d.append(se[i][j]+se[i][j+1])
    for m in dp1:
        for o in range(len(d)):
            if d[o]==m:
                n=n+1
    p=(n/(len(se[i])-1))*100
    data.iloc[i,:][h_class[0]]=p
for i in range(len(data)):
    n=0
    d=[]
    for j in range(len(se[i])-1):
        d.append(se[i][j]+se[i][j+1])
    for m in dp2:
        for o in range(len(d)):
            if d[o]==m:
                n=n+1
    p=(n/(len(se[i])-1))*100
    data.iloc[i,:][h_class[1]]=p
for i in range(len(data)):
    n=0
    d=[]
    for j in range(len(se[i])-1):
        d.append(se[i][j]+se[i][j+1])
    for m in dp3:
        for o in range(len(d)):
            if d[o]==m:
                n=n+1
    p=(n/(len(se[i])-1))*100
    data.iloc[i,:][h_class[2]]=p
    
##NWV_class
c0=['G','A','S','T','P','D']
c1=['N','V','E','Q','I','L']
c2=['M','H','K','F','R','Y','W']
dp1=[]##12 21
dp2=[]##13 31
dp3=[]##23 32
for i in c0:
    for j in c1:
        dp1.append(i+j)
for i in c0:
    for j in c2:
        dp2.append(i+j)
for i in c1:
    for j in c2:
        dp3.append(i+j)
for i in c1:
    for j in c0:
        dp1.append(i+j)
for i in c2:
    for j in c0:
        dp2.append(i+j)
for i in c2:
    for j in c1:
        dp3.append(i+j)
for i in range(len(data)):
    n=0
    d=[]
    for j in range(len(se[i])-1):
        d.append(se[i][j]+se[i][j+1])
    for m in dp1:
        for o in range(len(d)):
            if d[o]==m:
                n=n+1
    p=(n/(len(se[i])-1))*100
    data.iloc[i,:][NWV_class[0]]=p
for i in range(len(data)):
    n=0
    d=[]
    for j in range(len(se[i])-1):
        d.append(se[i][j]+se[i][j+1])
    for m in dp2:
        for o in range(len(d)):
            if d[o]==m:
                   n=n+1
    p=(n/(len(se[i])-1))*100
    data.iloc[i,:][NWV_class[1]]=p
for i in range(len(data)):
    n=0
    d=[]
    for j in range(len(se[i])-1):
        d.append(se[i][j]+se[i][j+1])
    for m in dp3:
        for o in range(len(d)):
            if d[o]==m:
                n=n+1
    p=(n/(len(se[i])-1))*100
    data.iloc[i,:][NWV_class[2]]=p

##Polarity_class
c0=['L','I','F','W','C','M','V','Y']
c1=['P','A','T','G','S']
c2=['H','Q','R','K','N','E','D']
dp1=[]##12 21
dp2=[]##13 31
dp3=[]##23 32
for i in c0:
    for j in c1:
        dp1.append(i+j)
for i in c0:
    for j in c2:
        dp2.append(i+j)
for i in c1:
    for j in c2:
        dp3.append(i+j)
for i in c1:
    for j in c0:
        dp1.append(i+j)
for i in c2:
    for j in c0:
        dp2.append(i+j)
for i in c2:
    for j in c1:
        dp3.append(i+j)
for i in range(len(data)):
    n=0
    d=[]
    for j in range(len(se[i])-1):
        d.append(se[i][j]+se[i][j+1])
    for m in dp1:
        for o in range(len(d)):
            if d[o]==m:
                n=n+1
    p=(n/(len(se[i])-1))*100
    data.iloc[i,:][Polarity_class[0]]=p
for i in range(len(data)):
    n=0
    d=[]
    for j in range(len(se[i])-1):
        d.append(se[i][j]+se[i][j+1])
    for m in dp2:
        for o in range(len(d)):
            if d[o]==m:
                n=n+1
    p=(n/(len(se[i])-1))*100
    data.iloc[i,:][Polarity_class[1]]=p
for i in range(len(data)):
    n=0
    d=[]
    for j in range(len(se[i])-1):
        d.append(se[i][j]+se[i][j+1])
    for m in dp3:
        for o in range(len(d)):
            if d[o]==m:
                n=n+1
    p=(n/(len(se[i])-1))*100
    data.iloc[i,:][Polarity_class[2]]=p

##Polarizability_class
c0=['G','A','S','D','T']
c1=['C','P','N','V','E','Q','I','L']
c2=['K','M','H','F','R','Y','W']
dp1=[]##12 21
dp2=[]##13 31
dp3=[]##23 32
for i in c0:
    for j in c1:
        dp1.append(i+j)
for i in c0:
    for j in c2:
        dp2.append(i+j)
for i in c1:
    for j in c2:
        dp3.append(i+j)
for i in c1:
    for j in c0:
        dp1.append(i+j)
for i in c2:
    for j in c0:
        dp2.append(i+j)
for i in c2:
    for j in c1:
        dp3.append(i+j)
for i in range(len(data)):
    n=0
    d=[]
    for j in range(len(se[i])-1):
        d.append(se[i][j]+se[i][j+1])
    for m in dp1:
        for o in range(len(d)):
            if d[o]==m:
                n=n+1
    p=(n/(len(se[i])-1))*100
    data.iloc[i,:][Polarizability_class[0]]=p
for i in range(len(data)):
    n=0
    d=[]
    for j in range(len(se[i])-1):
        d.append(se[i][j]+se[i][j+1])
    for m in dp2:
        for o in range(len(d)):
            if d[o]==m:
                n=n+1
    p=(n/(len(se[i])-1))*100
    data.iloc[i,:][Polarizability_class[1]]=p
for i in range(len(data)):
    n=0
    d=[]
    for j in range(len(se[i])-1):
        d.append(se[i][j]+se[i][j+1])
    for m in dp3:
        for o in range(len(d)):
            if d[o]==m:
                n=n+1
    p=(n/(len(se[i])-1))*100
    data.iloc[i,:][Polarizability_class[2]]=p

##Charge_class
c0=['K','R']
c1=['A','N','C','Q','G','H','I','L','M','F','P','S','T','W','Y','V']
c2=['D','E']
dp1=[]##12 21
dp2=[]##13 31
dp3=[]##23 32
for i in c0:
    for j in c1:
        dp1.append(i+j)
for i in c0:
    for j in c2:
        dp2.append(i+j)
for i in c1:
    for j in c2:
        dp3.append(i+j)
for i in c1:
    for j in c0:
        dp1.append(i+j)
for i in c2:
    for j in c0:
        dp2.append(i+j)
for i in c2:
    for j in c1:
        dp3.append(i+j)
for i in range(len(data)):
    n=0
    d=[]
    for j in range(len(se[i])-1):
        d.append(se[i][j]+se[i][j+1])
    for m in dp1:
        for o in range(len(d)):
            if d[o]==m:
                n=n+1
    p=(n/(len(se[i])-1))*100
    data.iloc[i,:][Charge_class[0]]=p
for i in range(len(data)):
    n=0
    d=[]
    for j in range(len(se[i])-1):
        d.append(se[i][j]+se[i][j+1])
    for m in dp2:
        for o in range(len(d)):
            if d[o]==m:
                n=n+1
    p=(n/(len(se[i])-1))*100
    data.iloc[i,:][Charge_class[1]]=p
for i in range(len(data)):
    n=0
    d=[]
    for j in range(len(se[i])-1):
        d.append(se[i][j]+se[i][j+1])
    for m in dp3:
        for o in range(len(d)):
            if d[o]==m:
                n=n+1
    p=(n/(len(se[i])-1))*100
    data.iloc[i,:][Charge_class[2]]=p
    
##Secondary_structure
c0=['E','A','L','M','Q','K','R','H']
c1=['V','I','Y','C','W','F','T']
c2=['G','N','P','S','D']
dp1=[]##12 21
dp2=[]##13 31
dp3=[]##23 32
for i in c0:
    for j in c1:
        dp1.append(i+j)
for i in c0:
    for j in c2:
        dp2.append(i+j)
for i in c1:
    for j in c2:
        dp3.append(i+j)
for i in c1:
    for j in c0:
        dp1.append(i+j)
for i in c2:
    for j in c0:
        dp2.append(i+j)
for i in c2:
    for j in c1:
        dp3.append(i+j)
for i in range(len(data)):
    n=0
    d=[]
    for j in range(len(se[i])-1):
        d.append(se[i][j]+se[i][j+1])
    for m in dp1:
        for o in range(len(d)):
            if d[o]==m:
                n=n+1
    p=(n/(len(se[i])-1))*100
    data.iloc[i,:][Secondary_structure[0]]=p
for i in range(len(data)):
    n=0
    d=[]
    for j in range(len(se[i])-1):
        d.append(se[i][j]+se[i][j+1])
    for m in dp2:
        for o in range(len(d)):
            if d[o]==m:
                n=n+1
    p=(n/(len(se[i])-1))*100
    data.iloc[i,:][Secondary_structure[1]]=p
for i in range(len(data)):
    n=0
    d=[]
    for j in range(len(se[i])-1):
        d.append(se[i][j]+se[i][j+1])
    for m in dp3:
        for o in range(len(d)):
            if d[o]==m:
                n=n+1
    p=(n/(len(se[i])-1))*100
    data.iloc[i,:][Secondary_structure[2]]=p
    
#Solvent_accessibility
c0=['A','L','F','C','G','I','V','W']
c1=['P','K','Q','E','N','D']
c2=['M','P','S','T','H','Y']
dp1=[]##12 21
dp2=[]##13 31
dp3=[]##23 32
for i in c0:
    for j in c1:
        dp1.append(i+j)
for i in c0:
    for j in c2:
        dp2.append(i+j)
for i in c1:
    for j in c2:
        dp3.append(i+j)
for i in c1:
    for j in c0:
        dp1.append(i+j)
for i in c2:
    for j in c0:
        dp2.append(i+j)
for i in c2:
    for j in c1:
        dp3.append(i+j)
for i in range(len(data)):
    n=0
    d=[]
    for j in range(len(se[i])-1):
        d.append(se[i][j]+se[i][j+1])
    for m in dp1:
        for o in range(len(d)):
            if d[o]==m:
                n=n+1
    p=(n/(len(se[i])-1))*100
    data.iloc[i,:][Solvent_accessibility[0]]=p
for i in range(len(data)):
    n=0
    d=[]
    for j in range(len(se[i])-1):
        d.append(se[i][j]+se[i][j+1])
    for m in dp2:
        for o in range(len(d)):
            if d[o]==m:
                n=n+1
    p=(n/(len(se[i])-1))*100
    data.iloc[i,:][Solvent_accessibility[1]]=p
for i in range(len(data)):
    n=0
    d=[]
    for j in range(len(se[i])-1):
        d.append(se[i][j]+se[i][j+1])
    for m in dp3:
        for o in range(len(d)):
            if d[o]==m:
                n=n+1
    p=(n/(len(se[i])-1))*100
    data.iloc[i,:][Solvent_accessibility[2]]=p
data=data.set_index('id')   
data.to_csv(Out+'/'+File.split('.')[0]+'_transition.csv')