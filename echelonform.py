def echelonform(x):
   n,m=np.shape(x)[0],np.shape(x)[1]#n=len(x)
   for kk in range(n-1): 
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
    for i in range(metritis,n):
        if ((x[i,kk]==1 and x[i,n+kk]==0) or (x[i,kk]==1 and x[i,n+kk]==1) or (x[i,kk]==0 and x[i,n+kk]==1)):
            
            for j in range(metritis,n):
                if (x[j,kk]==0 and x[j,n+kk]==0):
                    x[[j,i]]=x[[i,j]]
   
    for i in range(metritis,n):
        if ((x[i,kk]==1 and x[i,n+kk]==0) or (x[i,kk]==1 and x[i,n+kk]==1) or (x[i,kk]==0 and x[i,n+kk]==1)):
            metritis+=1 
   
   return x,list(reversed(tabluetostab(x)))
