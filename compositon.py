# -*- coding: utf-8 -*-
"""
Created on Mon May  9 07:25:45 2022

@author: zjy52
"""

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
h_class=['C_Hydrophobicity_Polar','C_Hydrophobicity_Neutral','C_Hydrophobicity_Hydrophobicity']
NWV_class=['C_NWV_Class1','C_NWV_Class2','C_NWV_Class3']
Polarity_class=['C_Polarity_Class1','C_Polarity_Class2','C_Polarity_Class3']
Polarizability_class=['C_Polarizability_Class1','C_Polarizability_Class2','C_Polarizability_Class3']
Charge_class=['C_Positive','C_Neutral','C_Negative']
Secondary_structure=['C_Helix','C_Strand','C_Coil']
Solvent_accessibility=['C_Buried','C_Exposed','C_Intermediate']

#dataframe
co=['id']+h_class+NWV_class+Polarity_class+Polarizability_class+Charge_class+Secondary_structure+Solvent_accessibility
data=pd.DataFrame(index=a,columns=co)
data['id']=name
##h_class
c0=['R','K','E','D','Q','N']
c1=['G','A','S','T','P','H','Y']
c2=['C','L','V','I','M','F','W']
for i in range(len(data)):
    n=0
    for j in c0:
        for m in range(len(se[i])):
            if se[i][m]==j:
                n=n+1
    s=n/len(se[i])*100
    data.iloc[i,:][h_class[0]]=s
for i in range(len(data)):
    n=0
    for j in c1:
        for m in range(len(se[i])):
            if se[i][m]==j:
                n=n+1
    s=n/len(se[i])*100
    data.iloc[i,:][h_class[1]]=s
for i in range(len(data)):
    n=0
    for j in c2:
        for m in range(len(se[i])):
            if se[i][m]==j:
                n=n+1
    s=n/len(se[i])*100
    data.iloc[i,:][h_class[2]]=s
##NWV_class
c0=['G','A','S','T','P','D']
c1=['N','V','E','Q','I','L']
c2=['M','H','K','F','R','Y','W']
for i in range(len(data)):
    n=0
    for j in c0:
        for m in range(len(se[i])):
            if se[i][m]==j:
                n=n+1
    s=n/len(se[i])*100
    data.iloc[i,:][NWV_class[0]]=s
for i in range(len(data)):
    n=0
    for j in c1:
        for m in range(len(se[i])):
            if se[i][m]==j:
                n=n+1
    s=n/len(se[i])*100
    data.iloc[i,:][NWV_class[1]]=s
for i in range(len(data)):
    n=0
    for j in c2:
        for m in range(len(se[i])):
            if se[i][m]==j:
                n=n+1
    s=n/len(se[i])*100
    data.iloc[i,:][NWV_class[2]]=s

##Polarity_class
c0=['L','I','F','W','C','M','V','Y']
c1=['P','A','T','G','S']
c2=['H','Q','R','K','N','E','D']
for i in range(len(data)):
    n=0
    for j in c0:
        for m in range(len(se[i])):
            if se[i][m]==j:
                n=n+1
    s=n/len(se[i])*100
    data.iloc[i,:][Polarity_class[0]]=s
for i in range(len(data)):
    n=0
    for j in c1:
        for m in range(len(se[i])):
            if se[i][m]==j:
                n=n+1
    s=n/len(se[i])*100
    data.iloc[i,:][Polarity_class[1]]=s
for i in range(len(data)):
    n=0
    for j in c2:
        for m in range(len(se[i])):
            if se[i][m]==j:
                n=n+1
    s=n/len(se[i])*100
    data.iloc[i,:][Polarity_class[2]]=s
    
##Polarizability_class
c0=['G','A','S','D','T']
c1=['C','P','N','V','E','Q','I','L']
c2=['K','M','H','F','R','Y','W']
for i in range(len(data)):
    n=0
    for j in c0:
        for m in range(len(se[i])):
            if se[i][m]==j:
                n=n+1
    s=n/len(se[i])*100
    data.iloc[i,:][Polarizability_class[0]]=s
for i in range(len(data)):
    n=0
    for j in c1:
        for m in range(len(se[i])):
            if se[i][m]==j:
                n=n+1
    s=n/len(se[i])*100
    data.iloc[i,:][Polarizability_class[1]]=s
for i in range(len(data)):
    n=0
    for j in c2:
        for m in range(len(se[i])):
            if se[i][m]==j:
                n=n+1
    s=n/len(se[i])*100
    data.iloc[i,:][Polarizability_class[2]]=s

##Charge_class
c0=['K','R']
c1=['A','N','C','Q','G','H','I','L','M','F','P','S','T','W','Y','V']
c2=['D','E']
for i in range(len(data)):
    n=0
    for j in c0:
        for m in range(len(se[i])):
            if se[i][m]==j:
                n=n+1
    s=n/len(se[i])*100
    data.iloc[i,:][Charge_class[0]]=s
for i in range(len(data)):
    n=0
    for j in c1:
        for m in range(len(se[i])):
            if se[i][m]==j:
                n=n+1
    s=n/len(se[i])*100
    data.iloc[i,:][Charge_class[1]]=s
for i in range(len(data)):
    n=0
    for j in c2:
        for m in range(len(se[i])):
            if se[i][m]==j:
                n=n+1
    s=n/len(se[i])*100
    data.iloc[i,:][Charge_class[2]]=s

##Secondary_structure
c0=['E','A','L','M','Q','K','R','H']
c1=['V','I','Y','C','W','F','T']
c2=['G','N','P','S','D']
for i in range(len(data)):
    n=0
    for j in c0:
        for m in range(len(se[i])):
            if se[i][m]==j:
                n=n+1
    s=n/len(se[i])*100
    data.iloc[i,:][Secondary_structure[0]]=s
for i in range(len(data)):
    n=0
    for j in c1:
        for m in range(len(se[i])):
            if se[i][m]==j:
                n=n+1
    s=n/len(se[i])*100
    data.iloc[i,:][Secondary_structure[1]]=s
for i in range(len(data)):
    n=0
    for j in c2:
        for m in range(len(se[i])):
            if se[i][m]==j:
                n=n+1
    s=n/len(se[i])*100
    data.iloc[i,:][Secondary_structure[2]]=s

#Solvent_accessibility
c0=['A','L','F','C','G','I','V','W']
c1=['P','K','Q','E','N','D']
c2=['M','P','S','T','H','Y']
for i in range(len(data)):
    n=0
    for j in c0:
        for m in range(len(se[i])):
            if se[i][m]==j:
                n=n+1
    s=n/len(se[i])*100
    data.iloc[i,:][Solvent_accessibility[0]]=s
for i in range(len(data)):
    n=0
    for j in c1:
        for m in range(len(se[i])):
            if se[i][m]==j:
                n=n+1
    s=n/len(se[i])*100
    data.iloc[i,:][Solvent_accessibility[1]]=s
for i in range(len(data)):
    n=0
    for j in c2:
        for m in range(len(se[i])):
            if se[i][m]==j:
                n=n+1
    s=n/len(se[i])*100
    data.iloc[i,:][Solvent_accessibility[2]]=s
data=data.set_index('id')
data.to_csv(Out+'/'+File.split('.')[0]+'_composition.csv')