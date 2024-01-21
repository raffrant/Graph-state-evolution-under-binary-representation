import numpy as np
import networkx as nx
import itertools
import matplotlib.pyplot as plt
import sys
#import measurementssqgtableau
import time
import random
from scipy.sparse import csr_matrix

def stabtoadj(x):
  ss=measurementssqgtableau.tabfromstring(len(x),x)
  n=np.shape(x)[0]
  adj=np.zeros((n,n),dtype=int)   
  for i in range(n):
        adj[i]=ss[i,n:2*n]
  return adj   
def adjtostab(adj):
    n=np.shape(adj)[0]
    stab=[['+'] for i in range(n)]
    for i in range(n):
     for j in range(n):
        if i==j: 
         stab[i].append('X')
        else:
            if adj[i][j]==0:
              stab[i].append('I')
            elif   adj[i][j]==1:
              stab[i].append('Z')
            else:
                pass
    for i in range(len(stab)):
       stab[i]= (''.join(stab[i][j] for j in range(len(stab)+1)))
    return stab
def randomstab(n):
    stab=[['+'] for i in range(n)]
    vv=['X','Z','Y','I']
    for j in range(n):
     for i in range(n):
        stab[i].append(random.choice(vv))
    for i in range(len(stab)):
       stab[i]= (''.join(stab[i][j] for j in range(len(stab)+1)))
    
    return stab 

ra=(randomstab(10))

def perm(N,tel,arx):
    return (csr_matrix(([1]*N, (tel, arx))).toarray())

print(ra)


ss=measurementssqgtableau.tabfromstring(len(ra),ra)
n=len(ss)

def rankcheck(a,j):                               #The function of gaussian elimination of the state
             n=np.array(a).shape[0]
             kkk=0
             vv=[]
             z=np.zeros(n)
             for n in range(n+1):    
                 vv.append([[1 if i == j else 0 for i in range(n)] for j in z])
             aas=[]
             for i in range(1,n+1):
                aas.append(vv[i][0])
             aas=(list(reversed(aas)))
             while(np.any(np.not_equal(a[j:,j], aas[j])) and kkk<=n):
              if (a[j,j]==1):
                     for k in range(j+1,n):
                         if a[k,:][j]==1:
                             a[k,:]=a[k,:]^a[j,:]
              if (a[j,j]==0):
                   a=measurementssqgtableau.had(n,j,a)  
                   ka=j
                   while ((ka<n-1) and (a[j,j]==0)):
                     a[[ka+1,j]]=a[[j,ka+1]]
                     ka+=1
              kkk+=1
           
             return a
start=time.time()       
#bb=allrank[5]  
for j in range(n):
  ss=rankcheck(ss, j)
  print(j)
print(ss)
end=time.time()
#print(end-start)
plt.imshow(ss)


def tabluetostab(x):
    n,m=np.shape(x)[0],np.shape(x)[1]
    stab=[[] for i in range(n)]

    for i in range(n):
        if x[i][-1]==0:
            stab[i].append('+')
        else:
            stab[i].append('-')
    for i in range(n):
     for j in range(n):
       if x[i,j]==0 and x[i,j+n]==0:
           stab[i].append('I')
       if x[i,j]==1 and x[i,j+n]==0:
           stab[i].append('X')
       if x[i,j]==0 and x[i,j+n]==1:
           stab[i].append('Z')
       if x[i,j]==1 and x[i,j+n]==1:
           stab[i].append('Y')    
    for i in range(len(stab)):
       stab[i]= (''.join(stab[i][j] for j in range(len(stab)+1)))
    return stab 
cc=(tabluetostab(ss))   


def echelonform(x):                                     #Put stabilizers in specific group, echelon form
   n,m=np.shape(x)[0],np.shape(x)[1]#n=len(x)
   for kk in range(1): 
    wh=[]
    xx=[]
    yy=[]
    zz=[]
    if kk==0:
        metritis=0
    else:
        pass
    for i in range(metritis,n):
        if ((x[i,kk]==1 and x[i,n+kk]==0) or (x[i,kk]==1 and x[i,n+kk]==1) or (x[i,kk]==0 and x[i,n+kk]==1)):
            if (x[i,kk]==1 and x[i,n+kk]==0):
                xx.append(i)
                wh.append(i)
            if (x[i,kk]==1 and x[i,n+kk]==1):
                yy.append(i)
                wh.append(i)    
            if (x[i,kk]==0 and x[i,n+kk]==1):
                zz.append(i)
                wh.append(i)
    print(wh,xx,yy,zz)
    if len(xx)>1:
     for i in range(1,len(xx)):
        x[xx[i]]=x[xx[i]]^x[xx[0]]
    if len(yy)>1:
     for i in range(1,len(yy)):
        x[yy[i]]=x[yy[i]]^x[yy[0]]
    if len(zz)>1:
     for i in range(1,len(zz)):
        x[zz[i]]=x[zz[i]]^x[zz[0]]    
    #print(x) 
    if True:
     for i in range(n):
        if ((x[i,kk]==1 and x[i,n+kk]==0) or (x[i,kk]==1 and x[i,n+kk]==1) or (x[i,kk]==0 and x[i,n+kk]==1)):
            
            for j in range(kk,n):
                if (x[j,kk]==0 and x[j,n+kk]==0):
                    x[[j,i]]=x[[i,j]]

   return x,list(reversed(tabluetostab(x)))
rref01b,ectab=(echelonform(ss))
telstab2=[]
for i in range(n):                              #Erase the phase of the stabilizers to calculate the height function
   telstab2.append(ectab[i].replace("+" or "-", ""))
  
def heig(x):
    n=len(x)
    heig=[]
    for i in range(n):
        heig.append(n-(i+1))
    numb=np.zeros(n,dtype=int)    
    
    for i in range(n):
     kk=0
     for j in range(n):
         if kk<1:
          if x[i][j]!="I":
             print(x[i][j]) 
             numb[i]=j+1
             kk+=1
    numberall=np.zeros(n,dtype=int)
    for i in range(len(numb)):
      for j in range(len(numb)):
          if numb[j]>i+1:
            numberall[i]+=1
    #print(numberall)        
    return heig-numberall  


print(heig(telstab2))

