# Author: Cameron F Abrams <cfa22@drexel.edu>
import importlib
import os
import numpy as np
import pandas as pd
from .util import pformatter, add_headers
from scipy.interpolate import interp1d

class SATD:
    def __init__(self):
        with importlib.resources.path('sandlersteam','__init__.py') as f:
            inst_root=os.path.split(os.path.abspath(f.parent.parent))[0]
        pfn=['SandlerSatdSteamTableP1.txt','SandlerSatdSteamTableP2.txt']
        tfn=['SandlerSatdSteamTableT1.txt','SandlerSatdSteamTableT2.txt']
        punits=['kPa','MPa']
        colorder=['T','P','VL','VV','UL','DU','UV','HL','DH','HV','SL','DS','SV']
        self.colFirstTwoTexLabels={}
        self.colFirstTwoTexLabels['T']=['$T$ (C)', '$P$ (MPa)']
        self.colFirstTwoTexLabels['P']=['$P$ (MPa)', '$T$ (C)']
        self.colRestTexLabels=[r'$\hat{V}^L$', r'$\hat{V}^V$', 
                        r'$\hat{U}^L$', r'$\Delta\hat{U}$', r'$\hat{U}^V$',
                        r'$\hat{H}^L$', r'$\Delta\hat{H}$', r'$\hat{H}^V$',
                        r'$\hat{S}^L$', r'$\Delta\hat{S}$', r'$\hat{S}^V$'
                      ]
        self.colFirstTwoTexFormatters={}
        self.colFirstTwoTexFormatters['T']=[pformatter,pformatter]
        self.colFirstTwoTexFormatters['P']=[pformatter,pformatter]
        self.DF={}
        self.lim={}
        D=[]
        for f,pu in zip(pfn,punits):
            data_abs_path=os.path.join(inst_root,'data',f)
            df=pd.read_csv(data_abs_path,sep=r'\s+',header=0,index_col=None)
            if pu=='kPa':
                df['P']=df['P']/1000.0
            ndf=df[colorder]
            # print(ndf.info())
            # print(ndf.head())
            D.append(ndf)
        self.DF['P']=pd.concat(D,axis=0)
        self.DF['P'].sort_values(by='P',inplace=True)
        # print(self.DF['P'].head())
        # print(self.DF['P'].tail())
        # print(self.DF['P'].info())
        self.lim['P']=[self.DF['P']['P'].min(),self.DF['P']['P'].max()]
        # print(self.lim)
        D=[]
        for f,pu in zip(tfn,punits):
            data_abs_path=os.path.join(inst_root,'data',f)
            df=pd.read_csv(data_abs_path,sep=r'\s+',header=0,index_col=None)
            if pu=='kPa':
                df['P']=df['P']/1000.0
            ndf=df[colorder]
            # print(ndf.head())
            D.append(ndf)
        self.DF['T']=pd.concat(D,axis=0)
        # print(self.DF['T'].describe())
        self.DF['T'].sort_values(by='T',inplace=True)
        self.lim['T']=[self.DF['T']['T'].min(),self.DF['T']['T'].max()]

        self.interpolators={}
        for bp,cp in zip(['P','T'],['T','P']):
            self.interpolators[bp]={}
            X=np.array(self.DF[bp][bp].to_list())
            for p in [cp,'VL','VV','UL','UV','HL','HV','SL','SV']:
                Y=np.array(self.DF[bp][p].to_list())
                self.interpolators[bp][p]=interp1d(X,Y,kind='linear')
                
    def texprint(self,**kwargs):
        by=kwargs.get('by','T')
        cp='T' if by=='P' else 'P'
        assert by in 'PT'
        block=self.DF[by]
        if not block.empty:
            splits=[block[block['T']<97.0],block[block['T']>97.0]]
            splits[0]['P']=splits[0]['P']*1000 # kPa from MPa
            strsplits=[]
            for bs,pu in zip(splits,['kPa','MPa']):
                block_floatsplit=pd.DataFrame()
                cols=[by,cp,'VL','VV','UL','DU','UV','HL','DH','HV','SL','DS','SV']
                # fmts=r'r@{}l'*len(cols)
                fmts =r'>{\raggedleft}p{4mm}@{}p{4mm}>{\raggedleft}p{4mm}@{}p{4mm}' # T/P, P/T
                fmts+=r'>{\raggedleft}p{2mm}@{}p{10mm}' # VL
                fmts+=r'>{\raggedleft}p{4mm}@{}p{10mm}'  # VV
                fmts+=r'>{\raggedleft}p{5mm}@{}p{3mm}'  # UL
                fmts+=r'>{\raggedleft}p{6mm}@{}p{2mm}'  # DU
                fmts+=r'>{\raggedleft}p{6mm}@{}p{2mm}'  # UV
                fmts+=r'>{\raggedleft}p{5mm}@{}p{3mm}'  # HL
                fmts+=r'>{\raggedleft}p{6mm}@{}p{2mm}'  # DH
                fmts+=r'>{\raggedleft}p{6mm}@{}p{2mm}'  # HV
                fmts+=r'>{\raggedleft}p{2mm}@{}p{6mm}'  # SL
                fmts+=r'>{\raggedleft}p{2mm}@{}p{6mm}'  # DS
                fmts+=r'>{\raggedleft\arraybackslash}p{2mm}@{}p{6mm}'  # SV

                hdgs=[]
                for c in cols:
                    hdgs.append(c)
                    hdgs.append('~')
                    W=np.floor(bs[c])
                    F=bs[c]-W
                    FS=[f'{x:.8f}'[1:] for x in F]
                    PFS=[]
                    block_floatsplit[c+'w']=W
                    for w,f,fs in zip(W,F,FS):
                        v=w+f
                        # ss is the explicit decimal part, need to choose digits
                        if c=='P' and by=='T':
                            # pressure digit rules when table is indexed by T:
                            # if v<0.2: min of 5 dp
                            if v<0.2:
                                while len(fs)>6 and fs[-1]=='0': fs=fs[:-1]
                            # elif v<2.0: min of 4 dp
                            elif v<2.0:
                                while len(fs)>5 and fs[-1]=='0': fs=fs[:-1]
                            # elif v<20: min of 3 dp
                            elif v<20.0:
                                while len(fs)>4 and fs[-1]=='0': fs=fs[:-1]
                            else:
                                while len(fs)>3 and fs[-1]=='0': fs=fs[:-1]
                        elif c=='P' and by=='P':
                            # pressure digit rules when table is indexed by P
                            if pu=='kPa':
                                if v<1.0:
                                    while len(fs)>5 and fs[-1]=='0': fs=fs[:-1]
                                elif v<10:
                                    while len(fs)>2 and fs[-1]=='0': fs=fs[:-1]
                                else:
                                    fs=''
                            else:
                                # if v<0.4, min of 3 dp
                                if v<0.4:
                                    while len(fs)>4 and fs[-1]=='0': fs=fs[:-1]
                                elif v<4.0:
                                    while len(fs)>3 and fs[-1]=='0': fs=fs[:-1]
                                else:
                                    if f==0.0: fs=''
                                    else:
                                        while len(fs)>3 and fs[-1]=='0': fs=fs[:-1]
                        elif c=='T' and by=='T':
                            if f==0.0: fs=''
                            else:
                                while len(fs)>3 and fs[-1]=='0': fs=fs[:-1]
                        elif c=='T' and by=='P':
                            while len(fs)>3 and fs[-1]=='0': fs=fs[:-1]
                        else:
                            if c=='VL':
                                while len(fs)>7 and fs[-1]=='0': fs=fs[:-1]
                                fs=fs[:4]+' '+fs[4:]
                            elif c=='VV':
                                if v>10:
                                    while len(fs)>3 and fs[-1]=='0': fs=fs[:-1]
                                elif v>2:
                                    while len(fs)>4 and fs[-1]=='0': fs=fs[:-1]
                                elif v>0.2:
                                    while len(fs)>5 and fs[-1]=='0': fs=fs[:-1]
                                elif v>0.02:
                                    while len(fs)>6 and fs[-1]=='0': fs=fs[:-1]
                                elif v>0.002:
                                    while len(fs)>7 and fs[-1]=='0': fs=fs[:-1]
                                if len(fs)==6:
                                    fs=fs[:3]+' '+fs[3:]
                                elif len(fs)==7:
                                    fs=fs[:4]+' '+fs[4:]

                            elif c=='UL' or c=='HL':
                                if v<1400:
                                    while len(fs)>3 and fs[-1]=='0': fs=fs[:-1]
                                else:
                                    while len(fs)>2 and fs[-1]=='0': fs=fs[:-1]
                            elif 'S' in c:
                                while len(fs)>5 and fs[-1]=='0': fs=fs[:-1]
                            else:
                                while len(fs)>2 and fs[-1]=='0': fs=fs[:-1]
                        PFS.append(fs)
                    block_floatsplit[c+'d']=PFS
                strsplits.append(block_floatsplit)
            title=r'\begin{minipage}{\textwidth}'+'\n'+r'\tiny'+'\n'+r'\begin{center}'+'\n'
            ht1=[r'\multicolumn{2}{c}{~}',r'\multicolumn{2}{c}{~}',r'\multicolumn{4}{c}{Specific Volume}',r'\multicolumn{6}{c}{Internal Energy}',r'\multicolumn{6}{c}{Enthalpy}',r'\multicolumn{6}{c}{Entropy}']
            htst11=r'\cmidrule(lr){5-8}\cmidrule(lr){9-14}\cmidrule(lr){15-20}\cmidrule(lr){21-26}'
            if by=='T':
                first2=[r'\multicolumn{2}{c}{Temp.}',r'\multicolumn{2}{c}{Press.}']
            else:
                first2=[r'\multicolumn{2}{c}{Press.}',r'\multicolumn{2}{c}{Temp.}']
            ht2=first2+[r'\multicolumn{2}{c}{Sat.}',r'\multicolumn{2}{c}{Sat.}',
            r'\multicolumn{2}{c}{Sat.}',r'\multicolumn{2}{c}{~}',r'\multicolumn{2}{c}{Sat.}',
            r'\multicolumn{2}{c}{Sat.}',r'\multicolumn{2}{c}{~}',r'\multicolumn{2}{c}{Sat.}',
            r'\multicolumn{2}{c}{Sat.}',r'\multicolumn{2}{c}{~}',r'\multicolumn{2}{c}{Sat.}']
            if by=='T':
                first2=[r'\multicolumn{2}{c}{($^\circ$C)}',r'\multicolumn{2}{c}{(kPa)}']
            else:
                first2=[r'\multicolumn{2}{c}{(kPa)}',r'\multicolumn{2}{c}{($^\circ$C)}']
            ht3=first2+[r'\multicolumn{2}{c}{Liquid}',r'\multicolumn{2}{c}{Vapor}',
            r'\multicolumn{2}{c}{Liquid}',r'\multicolumn{2}{c}{Evap.}',r'\multicolumn{2}{c}{Vapor}',
            r'\multicolumn{2}{c}{Liquid}',r'\multicolumn{2}{c}{Evap.}',r'\multicolumn{2}{c}{Vapor}',
            r'\multicolumn{2}{c}{Liquid}',r'\multicolumn{2}{c}{Evap.}',r'\multicolumn{2}{c}{Vapor}']
            if by=='T':
                first2=[r'\multicolumn{2}{c}{$T$}',r'\multicolumn{2}{c}{$P$}']
            else:
                first2=[r'\multicolumn{2}{c}{$P$}',r'\multicolumn{2}{c}{$T$}']
            ht4=first2+[r'\multicolumn{2}{c}{$\hat{V}^L$}',r'\multicolumn{2}{c}{$\hat{V}^V$}',
            r'\multicolumn{2}{c}{$\hat{U}^L$}',r'\multicolumn{2}{c}{$\Delta\hat{U}$}',r'\multicolumn{2}{c}{$\hat{U}^V$}',
            r'\multicolumn{2}{c}{$\hat{H}^L$}',r'\multicolumn{2}{c}{$\Delta\hat{H}$}',r'\multicolumn{2}{c}{$\hat{H}^V$}',
            r'\multicolumn{2}{c}{$\hat{S}^L$}',r'\multicolumn{2}{c}{$\Delta\hat{S}$}',r'\multicolumn{2}{c}{$\hat{S}^V$}']
            
            tbl1=strsplits[0].to_latex(escape=False,header=False,column_format=fmts,index=False,float_format='%g')
            tbl1=add_headers(tbl1,[ht1,ht2,ht3,ht4],[htst11,'','',''])
            # tbl1=set_width(tbl1)
            if by=='T':
                first2=[r'\multicolumn{2}{c}{~}',r'\multicolumn{2}{c}{MPa}']
            else:
                first2=[r'\multicolumn{2}{c}{MPa}',r'\multicolumn{2}{c}{~}']
            ht3=first2+[r'\multicolumn{22}{c}{~}']
            tbl2=strsplits[1].to_latex(escape=False,header=False,column_format=fmts,index=False,float_format='%g')
            tbl2=add_headers(tbl2,[ht3],[''])
            # tbl2=set_width(tbl2)
            return title+tbl1+r'\\'+'\n'+tbl2+r'\end{center}'+'\n'+r'\end{minipage}'+'\n'
        else:
            return None