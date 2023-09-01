# -*- coding: utf-8 -*-
"""
Created on Mon May  9 08:28:45 2022

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
##distibution_class
name=[]
se=[]
for seq in SeqIO.parse(File,'fasta'):
    name.append(seq.id)
    se.append(str(seq.seq))
a=np.linspace(0,len(name)-1,len(name))
h_class=['D_Hydrophobicity_Polar1','D_Hydrophobicity_Polar2', 'D_Hydrophobicity_Polar3','D_Hydrophobicity_Polar4',
         'D_Hydrophobicity_Polar5',
         'D_Hydrophobicity_Neutral1','D_Hydrophobicity_Neutral2','D_Hydrophobicity_Neutral3', 'D_Hydrophobicity_Neutral4',
         'D_Hydrophobicity_Neutral5',
         'D_Hydrophobicity_Hydrophobicity1',
         'D_Hydrophobicity_Hydrophobicity2',
         'D_Hydrophobicity_Hydrophobicity3',     
         'D_Hydrophobicity_Hydrophobicity4',
         'D_Hydrophobicity_Hydrophobicity5']
NWV_class=['D_NWV_Class1_1','D_NWV_Class1_2','D_NWV_Class1_3','D_NWV_Class1_4','D_NWV_Class1_5',
           'D_NWV_Class2_1','D_NWV_Class2_2','D_NWV_Class2_3','D_NWV_Class2_4','D_NWV_Class2_5',
           'D_NWV_Class3_1',
           'D_NWV_Class3_2',
           'D_NWV_Class3_3',
           'D_NWV_Class3_4',
           'D_NWV_Class3_5']
Polarity_class=['D_Polarity_Class1_1','D_Polarity_Class1_2','D_Polarity_Class1_3','D_Polarity_Class1_4','D_Polarity_Class1_5',
                'D_Polarity_Class2_1','D_Polarity_Class2_2','D_Polarity_Class2_3','D_Polarity_Class2_4','D_Polarity_Class2_5',
                'D_Polarity_Class3_1',
                'D_Polarity_Class3_2',
                'D_Polarity_Class3_3',
                'D_Polarity_Class3_4',
                'D_Polarity_Class3_5']
Polarizability_class=['D_Polarizability_Class1_1','D_Polarizability_Class1_2','D_Polarizability_Class1_3','D_Polarizability_Class1_4',
                     'D_Polarizability_Class1_5',
                     'D_Polarizability_Class2_1','D_Polarizability_Class2_2','D_Polarizability_Class2_3','D_Polarizability_Class2_4',
                   'D_Polarizability_Class2_5',   
                   'D_Polarizability_Class3_1',
                   'D_Polarizability_Class3_2',
                   'D_Polarizability_Class3_3',
                    'D_Polarizability_Class3_4',
                   'D_Polarizability_Class3_5']
Charge_class=['D_Positive1',
              'D_Positive2',
              'D_Positive3',
              'D_Positive4',
              'D_Positive5',
              'D_Neutral1',
              'D_Neutral2',
              'D_Neutral3',
              'D_Neutral4',
              'D_Neutral5',
              'D_Negative1',
              'D_Negative2',
              'D_Negative3',
              'D_Negative4',
              'D_Negative5']
Secondary_structure=['D_Helix1',
                     'D_Helix2',
                     'D_Helix3',
                     'D_Helix4',
                     'D_Helix5',
                     'D_Strand1',
                     'D_Strand2',
                     'D_Strand3',
                     'D_Strand4',
                     'D_Strand5',
                     'D_Coil1',
                     'D_Coil2',
                     'D_Coil3',
                     'D_Coil4',
                     'D_Coil5']
Solvent_accessibility=['D_Buried1',
                       'D_Buried2',
                       'D_Buried3',
                       'D_Buried4',
                       'D_Buried5',
                       'D_Exposed1',
                       'D_Exposed2',
                       'D_Exposed3',
                       'D_Exposed4',
                       'D_Exposed5',
                       'D_Intermediate1',
                       'D_Intermediate2',
                       'D_Intermediate3',
                       'D_Intermediate4',
                       'D_Intermediate5']

co=['id']+h_class+NWV_class+Polarity_class+Polarizability_class+Charge_class+Secondary_structure+Solvent_accessibility
data=pd.DataFrame(index=a,columns=co)
data['id']=name

#h_class
c0=['R','K','E','D','Q','N']
c1=['G','A','S','T','P','H','Y']
c2=['C','L','V','I','M','F','W']
for i in range(len(data)):
    s=[0,0.25,0.5,0.75,1]
    N1=[]
    N2=[]
    N3=[]
    for j in range(len(se[i])):
        k=se[i][j]
        for m in c0:
            if m==k:
                N1.append(j)
    for o in range(len(se[i])):
         k=se[i][o]
         for m in c1:
             if m==k:
                 N2.append(o)
    for t in range(len(se[i])):
         k=se[i][t]
         for m in c2:
             if m==k:
                 N3.append(t)
    if len(N1)>0:          
        for l in range(5):
                    data.iloc[i,:][h_class[l]]=int((N1[int(s[l]*(len(N1)-1))]+1)/len(N1)*100)
    else:
        for l in range(5):
                    data.iloc[i,:][h_class[l]]=0
    if len(N2)>0:
        for y in range(5,10):
                    data.iloc[i,:][h_class[y]]=int((N2[int(s[y-5]*(len(N2)-1))]+1)/len(N2)*100)
    else:
        data.iloc[i,:][h_class[y]]=0
    if len(N3)>0:    
        for q in range(10,15):
                    data.iloc[i,:][h_class[q]]=int((N3[int(s[q-10]*(len(N3)-1))]+1)/len(N3)*100)
    else:
        data.iloc[i,:][h_class[q]]=0

##NWV_class
c0=['G','A','S','T','P','D']
c1=['N','V','E','Q','I','L']
c2=['M','H','K','F','R','Y','W']
for i in range(len(data)):
    s=[0,0.25,0.5,0.75,1]
    N1=[]
    N2=[]
    N3=[]
    for j in range(len(se[i])):
        k=se[i][j]
        for m in c0:
            if m==k:
                N1.append(j)
    for j in range(len(se[i])):
         k=se[i][j]
         for m in c1:
             if m==k:
                 N2.append(j)
    for j in range(len(se[i])):
         k=se[i][j]
         for m in c2:
             if m==k:
                 N3.append(j)
    if len(N1)>0:          
        for l in range(5):
                    data.iloc[i,:][NWV_class[l]]=int((N1[int(s[l]*(len(N1)-1))]+1)/len(N1)*100)
    else:
        for l in range(5):
                    data.iloc[i,:][NWV_class[l]]=0
    if len(N2)>0:
        for y in range(5,10):
                    data.iloc[i,:][NWV_class[y]]=int((N2[int(s[y-5]*(len(N2)-1))]+1)/len(N2)*100)
    else:
        data.iloc[i,:][NWV_class[y]]=0
    if len(N3)>0:    
        for q in range(10,15):
                    data.iloc[i,:][NWV_class[q]]=int((N3[int(s[q-10]*(len(N3)-1))]+1)/len(N3)*100)
    else:
        data.iloc[i,:][NWV_class[q]]=0
        
##Polarity_class
c0=['L','I','F','W','C','M','V','Y']
c1=['P','A','T','G','S']
c2=['H','Q','R','K','N','E','D']
for i in range(len(data)):
    s=[0,0.25,0.5,0.75,1]
    N1=[]
    N2=[]
    N3=[]
    for j in range(len(se[i])):
        k=se[i][j]
        for m in c0:
            if m==k:
                N1.append(j)
    for j in range(len(se[i])):
         k=se[i][j]
         for m in c1:
             if m==k:
                 N2.append(j)
    for j in range(len(se[i])):
         k=se[i][j]
         for m in c2:
             if m==k:
                 N3.append(j)
    if len(N1)>0:          
        for l in range(5):
                    data.iloc[i,:][Polarity_class[l]]=int((N1[int(s[l]*(len(N1)-1))]+1)/len(N1)*100)
    else:
        for l in range(5):
                    data.iloc[i,:][Polarity_class[l]]=0
    if len(N2)>0:
        for y in range(5,10):
                    data.iloc[i,:][Polarity_class[y]]=int((N2[int(s[y-5]*(len(N2)-1))]+1)/len(N2)*100)
    else:
        data.iloc[i,:][Polarity_class[y]]=0
    if len(N3)>0:    
        for q in range(10,15):
                    data.iloc[i,:][Polarity_class[q]]=int((N3[int(s[q-10]*(len(N3)-1))]+1)/len(N3)*100)
    else:
        data.iloc[i,:][Polarity_class[q]]=0
        
##Polarizability_class
c0=['G','A','S','D','T']
c1=['C','P','N','V','E','Q','I','L']
c2=['K','M','H','F','R','Y','W']
for i in range(len(data)):
    s=[0,0.25,0.5,0.75,1]
    N1=[]
    N2=[]
    N3=[]
    for j in range(len(se[i])):
        k=se[i][j]
        for m in c0:
            if m==k:
                N1.append(j)
    for j in range(len(se[i])):
         k=se[i][j]
         for m in c1:
             if m==k:
                 N2.append(j)
    for j in range(len(se[i])):
         k=se[i][j]
         for m in c2:
             if m==k:
                 N3.append(j)
    if len(N1)>0:          
        for l in range(5):
                    data.iloc[i,:][Polarizability_class[l]]=int((N1[int(s[l]*(len(N1)-1))]+1)/len(N1)*100)
    else:
        for l in range(5):
                    data.iloc[i,:][Polarizability_class[l]]=0
    if len(N2)>0:
        for y in range(5,10):
                    data.iloc[i,:][Polarizability_class[y]]=int((N2[int(s[y-5]*(len(N2)-1))]+1)/len(N2)*100)
    else:
        data.iloc[i,:][Polarizability_class[y]]=0
    if len(N3)>0:    
        for q in range(10,15):
                    data.iloc[i,:][Polarizability_class[q]]=int((N3[int(s[q-10]*(len(N3)-1))]+1)/len(N3)*100)
    else:
        data.iloc[i,:][Polarizability_class[q]]=0
        
##Charge_class
c0=['K','R']
c1=['A','N','C','Q','G','H','I','L','M','F','P','S','T','W','Y','V']
c2=['D','E']
for i in range(len(data)):
    s=[0,0.25,0.5,0.75,1]
    N1=[]
    N2=[]
    N3=[]
    for j in range(len(se[i])):
        k=se[i][j]
        for m in c0:
            if m==k:
                N1.append(j)
    for j in range(len(se[i])):
         k=se[i][j]
         for m in c1:
             if m==k:
                 N2.append(j)
    for j in range(len(se[i])):
         k=se[i][j]
         for m in c2:
             if m==k:
                 N3.append(j)
    if len(N1)>0:          
        for l in range(5):
                    data.iloc[i,:][Charge_class[l]]=int((N1[int(s[l]*(len(N1)-1))]+1)/len(N1)*100)
    else:
        for l in range(5):
                    data.iloc[i,:][Charge_class[l]]=0
    if len(N2)>0:
        for y in range(5,10):
                    data.iloc[i,:][Charge_class[y]]=int((N2[int(s[y-5]*(len(N2)-1))]+1)/len(N2)*100)
    else:
        data.iloc[i,:][Charge_class[y]]=0
    if len(N3)>0:    
        for q in range(10,15):
                    data.iloc[i,:][Charge_class[q]]=int((N3[int(s[q-10]*(len(N3)-1))]+1)/len(N3)*100)
    else:
            data.iloc[i,:][Charge_class[q]]=0

##Secondary_structure
c0=['E','A','L','M','Q','K','R','H']
c1=['V','I','Y','C','W','F','T']
c2=['G','N','P','S','D']
for i in range(len(data)):
    s=[0,0.25,0.5,0.75,1]
    N1=[]
    N2=[]
    N3=[]
    for j in range(len(se[i])):
        k=se[i][j]
        for m in c0:
            if m==k:
                N1.append(j)
    for j in range(len(se[i])):
         k=se[i][j]
         for m in c1:
             if m==k:
                 N2.append(j)
    for j in range(len(se[i])):
         k=se[i][j]
         for m in c2:
             if m==k:
                 N3.append(j)
    if len(N1)>0:          
        for l in range(5):
                    data.iloc[i,:][Secondary_structure[l]]=int((N1[int(s[l]*(len(N1)-1))]+1)/len(N1)*100)
    else:
        for l in range(5):
                    data.iloc[i,:][Secondary_structure[l]]=0
    if len(N2)>0:
        for y in range(5,10):
                    data.iloc[i,:][Secondary_structure[y]]=int((N2[int(s[y-5]*(len(N2)-1))]+1)/len(N2)*100)
    else:
        data.iloc[i,:][Secondary_structure[y]]=0
    if len(N3)>0:    
        for q in range(10,15):
                    data.iloc[i,:][Secondary_structure[q]]=int((N3[int(s[q-10]*(len(N3)-1))]+1)/len(N3)*100)
    else:
        data.iloc[i,:][Secondary_structure[q]]=0

#Solvent_accessibility
c0=['A','L','F','C','G','I','V','W']
c1=['P','K','Q','E','N','D']
c2=['M','P','S','T','H','Y']
for i in range(len(data)):
    s=[0,0.25,0.5,0.75,1]
    N1=[]
    N2=[]
    N3=[]
    for j in range(len(se[i])):
        k=se[i][j]
        for m in c0:
            if m==k:
                N1.append(j)
    for j in range(len(se[i])):
         k=se[i][j]
         for m in c1:
             if m==k:
                 N2.append(j)
    for j in range(len(se[i])):
         k=se[i][j]
         for m in c2:
             if m==k:
                 N3.append(j)
    if len(N1)>0:          
        for l in range(5):
                    data.iloc[i,:][Solvent_accessibility[l]]=int((N1[int(s[l]*(len(N1)-1))]+1)/len(N1)*100)
    else:
        for l in range(5):
                    data.iloc[i,:][Solvent_accessibility[l]]=0
    if len(N2)>0:
        for y in range(5,10):
                    data.iloc[i,:][Solvent_accessibility[y]]=int((N2[int(s[y-5]*(len(N2)-1))]+1)/len(N2)*100)
    else:
        data.iloc[i,:][Solvent_accessibility[y]]=0
    if len(N3)>0:    
        for q in range(10,15):
                    data.iloc[i,:][Solvent_accessibility[q]]=int((N3[int(s[q-10]*(len(N3)-1))]+1)/len(N3)*100)
    else:
        data.iloc[i,:][Solvent_accessibility[q]]=0

data=data.fillna(0)
data=data.set_index('id')       
data.to_csv(Out+'/'+File.split('.')[0]+'_distribution.csv')