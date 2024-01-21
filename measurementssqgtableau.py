# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 10:56:02 2022

@author: rafai
"""

import numpy as np 
import matplotlib.pyplot as plt
import networkx as nx

def tabfromstring(N,stringtabl):                                     #Stabilizers to binary representation
    tab=np.zeros((N,2*N+1),dtype=int)
   # if np.shape(stringtabl)[0]!=np.shape(stringtabl)[1]:
   #     raise "Not correct number of stabilizer"
    for i in range(N):
        for j in range(1,N+1):
            if stringtabl[i][j]=='X':   #Write the string of stabilizers in tableau formalism
                tab[i][j-1]=1 
                tab[i][j-1+N]=0
            if stringtabl[i][j]=='Z':
                tab[i][j-1]=0 
                tab[i][j-1+N]=1  
            if stringtabl[i][j]=='Y':
                tab[i][j-1]=1
                tab[i][j-1+N]=1
            if stringtabl[i][j]=='I':
                tab[i][j-1]=0
                tab[i][j-1+N]=0 
    for j in range(N):
     if stringtabl[j][0]=='+':
                tab[j][-1]=0
     if stringtabl[j][0]=='-':
                tab[j][-1]=1    
                
    return tab


def had(N,qub,tab):                                      #Apply Clifford gates Hadamard, CNOT, S gate
    for i in range(N):
        if (tab[i][qub]==0 and tab[i][qub+N]==1):
          tab[i][qub]=1^tab[i][qub]
          tab[i][qub+N]=1^tab[i][qub+N]
        elif  (tab[i][qub]==1 and tab[i][qub+N]==0):
           tab[i][qub]=1^tab[i][qub]
           tab[i][qub+N]=1^tab[i][qub+N] 
        elif (tab[i][qub]==1 and tab[i][qub+N]==1):
          tab[i][-1]=1^tab[i][-1]
        elif (tab[i][qub]==0 and tab[i][qub+N]==0):
          pass
    return tab 
def ss(N,qub,tab):
    for i in range(N):
        if (tab[i][qub]==0 and tab[i][qub+N]==1):
            pass
        elif  (tab[i][qub]==1 and tab[i][qub+N]==0):
             #tab[i][qub]=1^tab[i][qub]
             tab[i][qub+N]=1^tab[i][qub+N]
        elif (tab[i][qub]==1 and tab[i][qub+N]==1):
          tab[i][-1]=1^tab[i][-1]
          tab[i][qub+N]=1^tab[i][qub+N]
    return tab 


def cnot(N,con,tar,tab):
    for i in range(N):
       if  (tab[i][con]==0 and tab[i][con+N]==0) and   (tab[i][tar]==0 and tab[i][tar+N]==0):
           pass
         #  print(221)

       elif  (tab[i][con]==0 and tab[i][con+N]==0) and   (tab[i][tar]==1 and tab[i][tar+N]==0):
           pass
         #  print(222)
       
       elif  (tab[i][con]==0 and tab[i][con+N]==0) and   (tab[i][tar]==1 and tab[i][tar+N]==1):
           tab[i][con+N]=1^tab[i][con+N]
          # print(223)
           
       elif  (tab[i][con]==0 and tab[i][con+N]==0) and   (tab[i][tar]==0 and tab[i][tar+N]==1):
           tab[i][con+N]=1^tab[i][con+N]
          # print(224)
           
       #x operator first
       elif (tab[i][con]==1 and tab[i][con+N]==0) and   (tab[i][tar]==0 and tab[i][tar+N]==0):
           tab[i][tar]=1^tab[i][tar]
           
       elif (tab[i][con]==1 and tab[i][con+N]==0) and   (tab[i][tar]==1 and tab[i][tar+N]==0):
           tab[i][tar]=1^tab[i][tar]

       elif (tab[i][con]==1 and tab[i][con+N]==0) and   (tab[i][tar]==1 and tab[i][tar+N]==1):
           tab[i][con+N]=1^tab[i][con+N]
           tab[i][tar]=1^tab[i][tar]

       elif (tab[i][con]==1 and tab[i][con+N]==0) and   (tab[i][tar]==0 and tab[i][tar+N]==1):
           tab[i][con+N]=1^tab[i][con+N]
           tab[i][tar]=1^tab[i][tar]  
           tab[i][-1]=1^tab[i][-1]

       #y operator first"    
       elif (tab[i][con]==1 and tab[i][con+N]==1) and   (tab[i][tar]==0 and tab[i][tar+N]==0):    
                     tab[i][tar]=1^tab[i][tar]

       elif (tab[i][con]==1 and tab[i][con+N]==1) and   (tab[i][tar]==1 and tab[i][tar+N]==0):
           tab[i][tar]=1^tab[i][tar]

       elif (tab[i][con]==1 and tab[i][con+N]==1) and   (tab[i][tar]==1 and tab[i][tar+N]==1):
           tab[i][con+N]=1^tab[i][con+N]
           tab[i][tar]=1^tab[i][tar]
           tab[i][-1]=1^tab[i][-1]

       elif (tab[i][con]==1 and tab[i][con+N]==1) and   (tab[i][tar]==0 and tab[i][tar+N]==1):
           tab[i][con+N]=1^tab[i][con+N]
           tab[i][tar]=1^tab[i][tar]

       #z operator first"
       elif  (tab[i][con]==0 and tab[i][con+N]==1) and   (tab[i][tar]==0 and tab[i][tar+N]==0):
           pass

       elif  (tab[i][con]==0 and tab[i][con+N]==1) and   (tab[i][tar]==1 and tab[i][tar+N]==0):
           pass
       elif  (tab[i][con]==0 and tab[i][con+N]==1) and   (tab[i][tar]==1 and tab[i][tar+N]==1):
           tab[i][con+N]=1^tab[i][con+N]
       elif  (tab[i][con]==0 and tab[i][con+N]==1) and   (tab[i][tar]==0 and tab[i][tar+N]==1):
           tab[i][con+N]=1^tab[i][con+N]
           
    return tab 

def measx(phase,N,qub,tab):                     #  Pauli mesaurements in X,Y,Z basis
    allnonz=[]
    allnony=[]
    allnonx=[]
    new=[]
    if phase=='-':
       for i in range(2*N+1):
         if i==qub : 
          new.append(1)
         elif  i==qub+N : 
          new.append(0)
         elif i==2*N: 
          new.append(1)
         else:
            new.append(0)

        #print(222)
    else:
        for i in range(2*N+1):
         if i==qub : 
          new.append(1)
         elif  i==qub+N : 
          new.append(0)
         elif i==2*N: 
          new.append(0)
         else:
            new.append(0)
    for i in range(N):
        if tab[i][qub]==1 and  tab[i][qub+N]==0:
           allnonx.append(i)
        elif tab[i][qub]==0 and  tab[i][qub+N]==1:
            allnonz.append(i)
        elif tab[i][qub]==1 and  tab[i][qub+N]==1:
            allnony.append(i)    
    if len(allnonz)>1:
     for i in range(len(allnonz)):
        for j in range(i+1,len(allnonz)):
            tab[allnonz[j]]=tab[allnonz[j]]^tab[allnonz[i]]
     
    if len(allnonz)>0 and len(allnonx)>0:
     tab[allnonz[0]]=tab[allnonx[0]]^new 
     tab[allnonx[0]]=new 
    for i in range(N):
        if tab[i][i]==0:
           tab=had(N,i,tab)
    
    return tab

def measy(phase,N,qub,tab):
    allnonz=[]
    allnony=[]
    allnonx=[]
    new=[]
    if phase=='-':
       for i in range(2*N+1):
         if i==qub : 
          new.append(1)
         elif  i==qub+N : 
          new.append(1)
         elif i==2*N: 
          new.append(1)
         else:
            new.append(0)

        #print(222)
    else:
        for i in range(2*N+1):
         if i==qub : 
          new.append(1)
         elif  i==qub+N : 
          new.append(1)
         elif i==2*N: 
          new.append(0)
         else:
            new.append(0)
    for i in range(N):
        if tab[i][qub]==1 and  tab[i][qub+N]==0:
           allnonx.append(i)
        elif tab[i][qub]==0 and  tab[i][qub+N]==1:
            allnonz.append(i)
        elif tab[i][qub]==1 and  tab[i][qub+N]==1:
            allnony.append(i)    
   # print(allnonx,'listnx')
   # print(allnonz,'listnz')
    if len(allnonz)>1:
     for i in range(len(allnonz)):
        for j in range(i+1,len(allnonz)):
            tab[allnonz[j]]=tab[allnonz[j]]^tab[allnonz[i]]   
    if len(allnonz)>0 and len(allnonx)>0:

     tab[allnonz[0]]=tab[allnonz[0]]^tab[allnonx[0]]
     tab[qub]=new              
    for i in range(N):
        if tab[i][i]==0 and tab[i][i+N]==1:
            tab=had(N,i,tab)

    return tab

def measz(phase,N,qub,tab):
    allnonz=[]
    allnony=[]
    allnonx=[]
    new=[]
    if phase=='-':
       for i in range(2*N+1):
         if i==qub : 
          new.append(0)
         elif  i==qub+N : 
          new.append(1)
         elif i==2*N: 
          new.append(1)
         else:
            new.append(0)

        #print(222)
    else:
        for i in range(2*N+1):
         if i==qub : 
          new.append(0)
         elif  i==qub+N : 
          new.append(1)
         elif i==2*N: 
          new.append(0)
         else:
            new.append(0)
    for i in range(N):
        if tab[i][qub]==1 and  tab[i][qub+N]==0:
           allnonx.append(i)
        elif tab[i][qub]==0 and  tab[i][qub+N]==1:
            allnonz.append(i)
        elif tab[i][qub]==1 and  tab[i][qub+N]==1:
            allnony.append(i)  
    tab[qub]=new              
    return tab,allnonx,allnony,allnonz           


        
        
     
