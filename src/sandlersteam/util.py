# Author: Cameron F. Abrams <cfa22@drexel.edu>
import io
import pandas as pd
import numpy as np

def add_headers(tblstr,hdllist,strs):
    tbllns=tblstr.split('\n')
    for i in range(len(tbllns)):
        if tbllns[i].startswith(r'\begin{tabular}'):
            break
    if i<len(tbllns):
        i+=1
        tbllns.insert(i,r'\toprule')
        i+=1
        for ln,st in zip(hdllist,strs):
            lstr=' & '.join(ln)
            tbllns.insert(i,lstr+r'\\')
            i+=1
            if len(st)>0:
                tbllns.insert(i,st)
                i+=1
        tblstr='\n'.join(tbllns)
    return tblstr

"""
                                                                                                    1         1
          1         2         3         4         5         6         7         8         9         0         1
012345678901234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123    
 260  0.0012749  1127.9  1134.3   2.8830  0.0012645  1121.1   1133.7   2.8699  0.0012550  1114.6   1133.4   2.8576
(0,4),(6,15),    (17,24),(25,32), (34,40),(42,51),   (53,60), (62,69), (71,77),(79,88),(90,97),(99,106),(108,114)
"""
def my_split(data,hder,P,Tsat,fixw=False):
    ndfs=[]
    with io.StringIO(data) as f:
        if fixw:
            df=pd.read_fwf(f,colspecs=((0,4),(6,15), (17,24),(25,32),(34,40),(42,51), (53,60),(62,69),(71,77),(79,88),(90,97),(99,106),(108,114)),header=None,index_col=None)
        else:
            df=pd.read_csv(f,sep=r'\s+',header=None,index_col=None)
        df.columns=hder
        i=1
        for p,ts in zip(P,Tsat):
            ndf=pd.DataFrame({'T':df['T'].copy(),'P':np.array([p for _ in range(df.shape[0])])})
            if ndf.iloc[0,0]=='Sat.':
                ndf.iloc[0,0]=ts
            ndf['T']=ndf['T'].astype(float)
            tdf=df.iloc[:,i:i+4].copy()
            i+=4
            ndf=pd.concat((ndf,tdf),axis=1)
            ndf.dropna(axis=0,inplace=True)
            ndf.sort_values(by='T',inplace=True)
            ndfs.append(ndf)
    mdf=pd.concat(ndfs,axis=0)
    return mdf

def pformatter(x,max_places=6):
    i=0
    while not x.is_integer() and i <= max_places:
        x*=10
        i+=1
    fmtstr=r'{:.'+str(i)+r'f}'
    return fmtstr.format
