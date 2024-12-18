"""
Read steam tables from Sandler
"""
import pandas as pd
import os
import io
import importlib.resources
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import LinearNDInterpolator, interp1d, interp2d

def pformatter(x,max_places=6):
    i=0
    while not x.is_integer() and i <= max_places:
        x*=10
        i+=1
    fmtstr=r'{:.'+str(i)+r'f}'
    return fmtstr.format

class SATD:
    def __init__(self):
        with importlib.resources.path('Sandlersteam','__init__.py') as f:
            inst_root=os.path.split(os.path.abspath(f))[0]
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
            # print(df.head())
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
            # print(ndf.to_string())
            ndfs.append(ndf)
    mdf=pd.concat(ndfs,axis=0)
    return mdf

class SUPH:
    def __init__(self,phase='V'):
        with importlib.resources.path('Sandlersteam','__init__.py') as f:
            inst_root=os.path.split(os.path.abspath(f))[0]
        if phase=='V':
            ff='SandlerSuphSteamTables.txt'
        elif phase=='L':
            ff='SandlerSubcSteamTables.txt'
        fn=os.path.join(inst_root,'data',ff)
        # read in entire text file
        with open(fn,'r') as f:
            lines=f.read().split('\n')
        # identify header
        hder=lines[0].split()
        # identify lines with pressures
        plines=[]
        # save list of pressure values
        Pval=[]
        for i,l in enumerate(lines):
            if 'P = ' in l:
                plines.append(i)
        # extract and format all blocks as concurrent dataframes
        DFS=[]
        for i,(l,r) in enumerate(zip(plines[:-1],plines[1:])):
            # get pressures, and saturation temperatures, if present
            tokens=lines[l].split()
            kills=['P','=','MPa']
            for k in kills:
                while k in tokens:
                    tokens.remove(k)
            P=[]
            Tsat=[]
            for x in tokens:
                if x[0]=='(':
                    y=x[1:-1]
                    Tsat.append(float(y))
                else:
                    P.append(float(x))
            while len(Tsat)<len(P):
                Tsat.append(None)
            Pval.extend(P)
            data='\n'.join(lines[l+1:r])
            ndf=my_split(data,hder,P,Tsat,fixw=(phase=='L'))
            # ndf=ndf.sort_values(by='T')
            # print(ndf.to_string())
            DFS.append(ndf)
        tokens=lines[plines[-1]].split()
        kills=['P','=','MPa']
        for k in kills:
            while k in tokens:
                tokens.remove(k)
        P=[]
        Tsat=[]
        for x in tokens:
            if x[0]=='(':
                y=x[1:-1]
                Tsat.append(float(y))
            else:
                P.append(float(x))
        while len(Tsat)<len(P):
            Tsat.append(None)
        Pval.extend(P)
        data='\n'.join(lines[plines[-1]+1:])
        ndf=my_split(data,hder,P,Tsat,fixw=(phase=='L'))
        # ndf=ndf.sort_values(by='T')
        # print(ndf.to_string())
        DFS.append(ndf)
        self.data=pd.concat(DFS,ignore_index=True)
        # print(self.data.to_string())
        # self.data.dropna(inplace=True,axis=0)
        # Pressures=self.data['P'].unique()
        # for P in Pressures:
        #     print(f'sorting T at P {P}')
        #     sorted_subframe=self.data[self.data['P']==P].sort_values(by='T').reset_index(drop=True)
        #     print(sorted_subframe.to_string())
        #     self.data[self.data['P']==P]=sorted_subframe
        # self.data.reset_index(drop=True,inplace=True)
        # print(self.data.to_string())
        # return
        # print(f'Superheated steam table: {self.data["T"].min()}<T(C)<{self.data["T"].max()} {self.data["P"].min()}<P(MPa)<{self.data["P"].max()} ')
        dof=self.data.columns
        # print(dof)
        self.uniqs={}
        for d in dof:
            self.uniqs[d]=np.sort(np.array(list(set(self.data[d].to_list()))))
        #     print(f'{d}')
        #     print(f'{self.uniqs[d]}')
        # exit()
        self.interpolators={}
        for i in range(len(dof)):
            X=np.array(self.data[dof[i]].to_list())
            for j in range(i+1,len(dof)):
                Y=np.array(self.data[dof[j]].to_list())
                # print(dof[i],dof[j])
                self.interpolators[f'{dof[i]}{dof[j]}']={}
                for k in range(len(dof)):
                    if k != i and k != j:
                        # print(f'   -> {dof[k]}')
                        Z=np.array(self.data[dof[k]].to_list())
                        # for x,y,z in zip(X,Y,Z):
                        #     print(x,y,z)
                        I=LinearNDInterpolator(list(zip(X,Y)),Z)
                        self.interpolators[f'{dof[i]}{dof[j]}'][dof[k]]=I
        # print(self.interpolators)
    def TPBilinear(self,specdict):
        retdict={}
        xn,yn=specdict.keys()
        assert [xn,yn]==['T','P']
        xi,yi=specdict.values()
        df=self.data
        dof=self.data.columns
        retdict={}
        retdict['T']=xi
        retdict['P']=yi
        # retdict['MATCH']=False
        tdf=df[df['P']==yi]
        if not tdf.empty:
            # print(f'Found SH block for P = {yi}')
            X=np.array(tdf['T'])
            # reidx=np.argsort(X)
            # X=X[reidx]
            # retdict['MATCH']=np.any([xi==x for x in X])
            # print(f'X {X}')
            for d in dof:
                if d not in 'TP':
                    Y=np.array(tdf[d])#[reidx]
                    retdict[d]=np.interp(xi,X,Y)
                    # print(f'interpolating {d} X {X} Y {Y} at x {xi} => {retdict[d]}')
            # print(f'returning {retdict}')
        else:
            for PL,PR in zip(self.uniqs['P'][:-1],self.uniqs['P'][1:]):
                if PL<yi<PR:
                    # print(f'Found {PL} < {yi} < {PR}')
                    break
            else:
                raise Exception(f'P {yi} not between {PL} and {PR}')
            X=np.array([PL,PR])
            LDF=df[df['P']==PL]
            RDF=df[df['P']==PR]
            LT=np.array(LDF['T'])
            RT=np.array(RDF['T'])
            # print(f'LT {LT}')
            # print(f'RT {RT}')
            CT=np.array([T for T in LT if T in RT])
            # print(f'Common T {CT}')
            if xi in CT:
                # print(f'Found line(s) for T {xi}')
                for d in dof:
                    if d not in 'TP':
                        Y=np.array([LDF[LT==xi][d].values[0],RDF[RT==xi][d].values[0]])
                        # print(f'X {X} Y {Y} yi {yi}')
                        retdict[d]=np.interp(yi,X,Y)
            else:
                for TL,TR in zip(CT[:-1],CT[1:]):
                    if TL<xi<TR:
                        # print(f'Found {TL} < {xi} < {TR}')
                        break
                else:
                    raise Exception(f'T {xi} not between {TL} and {TR}')
                LTDF=LDF[(LDF['T']==TL)|(LDF['T']==TR)].sort_values(by='T')
                RTDF=RDF[(RDF['T']==TL)|(RDF['T']==TR)].sort_values(by='T')
                # print(LTDF.to_string())
                # print(RTDF.to_string())
                iv=np.zeros(2)
                for p in dof:
                    if not p in 'TP':
                        for i in range(2):
                            Lp=LTDF[p].values[i]
                            Rp=RTDF[p].values[i]
                            Y=np.array([Lp,Rp])
                            # print(f'{p} {Lp} {Rp}')
                            iv[i]=np.interp(yi,X,Y)
                        retdict[p]=np.interp(xi,np.array([TL,TR]),iv)
        return retdict

    def texprint(self,P):
        block=self.data[self.data['P']==P][['T','V','U','H','S']]
        if not block.empty:
            block_floatsplit=pd.DataFrame()
            for c in ['T','V','U','H','S']:
                w=block[c].astype(int)
                dd=np.round((block[c]-w),10).astype(str)
                d=[]
                block_floatsplit[c+'w']=w
                for a,s in zip(w,dd):
                    if '.' in s:
                        ss=s[1:]
                    if ss=='.0' and c=='T':
                        d.append('')
                    else:
                        if a==0: # this is a fractional number
                            while len(ss)<6:
                                ss=ss+'0'
                            iss=int(ss[1:])
                            if iss>19999:
                                ss=ss[:-1]
                        elif a<10:
                            while len(ss)<5:
                                ss=ss+'0'
                        d.append(ss)
                    
                block_floatsplit[c+'d']=d
                # block_floatsplit[c+'d'][block_floatsplit[c+'d']==0]=None
            title=r'\begin{minipage}{0.6\textwidth}'+'\n'+r'\footnotesize\vspace{5mm}'+'\n'+r'\begin{center}'+'\n'+r'$P$ = '+f'{P}'+r' MPa\\*[1ex]'+'\n'
            fmts =r'>{\raggedleft}p{8mm}@{}p{5mm}' # T
            fmts+=r'>{\raggedleft}p{4mm}@{}p{10mm}' # V
            fmts+=r'>{\raggedleft}p{10mm}@{}p{3mm}' # U
            fmts+=r'>{\raggedleft}p{10mm}@{}p{3mm}' # H
            fmts+=r'>{\raggedleft\arraybackslash}p{3mm}@{}p{8mm}' # S
            tbl=block_floatsplit.to_latex(escape=False,header=False,column_format=fmts,index=False,float_format='%g')
            hdrs=[r'\multicolumn{2}{c}{$T$~($^\circ$C)}',
                  r'\multicolumn{2}{c}{$\hat{V}$}',
                  r'\multicolumn{2}{c}{$\hat{U}$}',
                  r'\multicolumn{2}{c}{$\hat{H}$}',
                  r'\multicolumn{2}{c}{$\hat{S}$}']
            tbl=add_headers(tbl,[hdrs],[''])
            return title+tbl+r'\end{center}'+'\n'+r'\end{minipage}'+'\n'
        else:
            return None

    def TThBilinear(self,specdict):
        xn,yn=specdict.keys()
        assert xn=='T'
        xi,yi=specdict.values()
        df=self.data
        dof=list(df.columns)
        dof.remove(yn)
        dof.remove(xn)
        retdict={}
        retdict['T']=xi
        retdict[yn]=yi
        LLdat={}
        LLdat[yn]=[]
        for d in dof:
            LLdat[d]=[]
        for P in self.uniqs['P'][::-1]:  # VUSH properties decrease with increasing P
            tdf=df[df['P']==P]
            X=np.array(tdf['T'])
            Y=np.array(tdf[yn])
            if Y.min()<yi<Y.max():
                LLdat['P'].append(P)
                for d in 'VUSH':
                    Y=np.array(tdf[d])
                    LLdat[d].append(np.interp(xi,X,Y))
        X=np.array(LLdat[yn])
        retdict['P']=np.interp(yi,X,np.array(LLdat['P']))
        for d in 'VUSH':
            if d!=yn:
                retdict[d]=np.interp(yi,X,LLdat[d])
        return retdict

    def PThBilinear(self,specdict):
        xn,yn=specdict.keys()
        xi,yi=specdict.values()
        assert xn=='P'
        df=self.data
        dof=list(df.columns)
        dof.remove(yn)
        dof.remove(xn)
        retdict={}
        retdict['T']=xi
        retdict[yn]=yi
        if xi in self.uniqs['P']:
            tdf=df[df['P']==xi]
            X=np.array(tdf[yn])
            for pp in dof:
                Y=np.array(tdf[pp])
                retdict[pp]=np.interp(yi,X,Y)
        else:
            for PL,PR in zip(self.uniqs['P'][:-1],self.uniqs['P'][1:]):
                if PL<xi<PR:
                    break
            else:
                raise Exception(f'Error: no two blocks bracket P={xi}')
            ldf,rdf=df[df['P']==PL],df[df['P']==PR]
            ldict,rdict={},{}
            ldict['P'],rdict['P']=PL,PR
            ldict[yn],rdict[yn]=yi,yi
            for d,xdf in zip([ldict,rdict],[ldf,rdf]):
                X=xdf[yn]
                for pp in dof:
                    Y=np.array(xdf[pp])
                    d[pp]=np.interp(yi,X,Y)
            X=np.array([PL,PR])
            for pp in dof:
                Y=np.array([ldict[pp],rdict[pp]])
                retdict[pp]=np.interp(xi,X,Y)
        # print(f'PThBilinear ret {retdict}')
        return retdict
    
    def ThThBilinear(self,specdict):
        xn,yn=specdict.keys()
        assert xn!='T' and yn!='T' and xn!='P' and yn!='P'
        xi,yi=specdict.values()
        df=self.data
        dof=self.data.columns
        LLdat={}
        for d in dof:
            if d not in 'TP':
                LLdat[d]=[]
        LLdat['T']=self.uniqs['T']
        for T in LLdat['T']:
            tdf=df[df['T']==T]
            X=np.array(tdf[xn])
            for d in dof:
                if d!='T' and d!=xn:
                    Y=np.array(tdf[d])
                    if Y.min()<yi<Y.max():
                        LLdat[d]=np.interp(xi,X,Y)
        X=LLdat[yn]
        retdict={}
        retdict[xn]=xi
        retdict[yn]=yi
        for d in dof:
            if d!=xn and d!=yn:
                retdict[d]=np.interp(yi,X,LLdat[d])
        return retdict

    def Bilinear(self,specdict):
        xn,yn=specdict.keys()
        if [xn,yn]==['T','P']:
            return self.TPBilinear(specdict)
        elif xn=='T':
            return self.TThBilinear(specdict)
        elif xn=='P':
            return self.PThBilinear(specdict)
        else:
            return self.ThThBilinear(specdict)

_SATD=SATD()
_SUPH=SUPH('V')
_SUBC=SUPH('L')

class PHASE:
    def __init__(self):
        pass

LARGE=1.e99
NEGLARGE=-LARGE

class SANDLER:
    _p=['T','P','u','v','s','h','x']
    _sp=['T','P','VL','VV','UL','UV','HL','HV','SL','SV']
    def _resolve(self):
        """ Resolve the thermodynamic state of steam/water given specifications
        """
        self.spec=[p for p in self._p if self.__dict__[p]]
        assert len(self.spec)==2,f'Error: must specify two properties (of {self._p}) for steam'
        if self.spec[1]=='x':
            ''' explicitly saturated '''
            self._resolve_satd()
        else:
            if self.spec==['T','P']:
                ''' T and P given explicitly '''
                self._resolve_subsup()
            elif 'T' in self.spec or 'P' in self.spec:
                ''' T OR P given, along with some other property (v,u,s,h) '''
                self._resolve_TPTh()
            else:
                raise Exception('If not explicitly saturated, you must specify either T or P')

    def _resolve_TPTh(self):
        ''' T or P along with one other property (th) are specified '''
        p=self.spec[0]
        cp='P' if p=='T' else 'T'
        th=self.spec[1]
        # print(f'{p} = {self.__dict__[p]} ({self.satd.lim[p][0]}-{self.satd.lim[p][1]}), {th} = {self.__dict__[th]}')
        if self.satd.lim[p][0]<self.__dict__[p]<self.satd.lim[p][1]:
            ''' T or P is between saturation limits; may be a saturated state, so 
                check whether the second property value lies between its liquid
                and vapor phase values at this T or P '''
            thL=self.satd.interpolators[p][f'{th.upper()}L'](self.__dict__[p])
            thV=self.satd.interpolators[p][f'{th.upper()}V'](self.__dict__[p])
            # print(f'{th}L = {thL}, {th}V = {thV}')
            self.__dict__[cp]=self.satd.interpolators[p][cp](self.__dict__[p])
            if thL<self.__dict__[th]<thV:
                ''' This is a saturated state! Use lever rule to get vapor fraction: '''
                self.x=(self.__dict__[th]-thL)/(thV-thL)
                self.Liquid=PHASE()
                self.Vapor=PHASE()
                self.Liquid.__dict__[th]=thL
                self.Vapor.__dict__[th]=thV
                for pp in self._sp:
                    if pp not in ['T','P',f'{th.upper()}V',f'{th.upper()}L']:
                        ppp=self.satd.interpolators[p][pp](self.__dict__[p])
                        if pp[-1]=='V':
                            self.Vapor.__dict__[pp[0].lower()]=ppp
                        elif pp[-1]=='L':
                            self.Liquid.__dict__[pp[0].lower()]=ppp
                for pp in self._p:
                    if pp not in [p,cp,'x']:
                        self.__dict__[pp]=self.x*self.Vapor.__dict__[pp]+(1-self.x)*self.Liquid.__dict__[pp]
            else:
                ''' even though T or P is between saturation limits, the other property is not '''
                specdict={p.upper():self.__dict__[p] for p in self.spec}
                if self.__dict__[th]<thL:
                    ''' Th is below its liquid-state value; assume this is a subcooled state '''
                    # icode=''.join([x.upper() for x in self.spec])
                    # dofv=[self.__dict__[p] for p in self.spec]
                    retdict=self.subc.Bilinear(specdict)
                    for p in self._p: 
                        if p not in self.spec and p!='x':
                            self.__dict__[p]=retdict[p.upper()]
                else:
                    ''' Th is above its vapor-state value; assume this is a superheated state '''
                    retdict=self.suph.Bilinear(specdict)
                    for p in self._p: 
                        if p not in self.spec and p!='x':
                            self.__dict__[p]=retdict[p.upper()]
        elif self.__dict__[p]>self.satd.lim[p][1]:
            ''' Th is above its vapor-state value; assume this is a superheated state '''
            # icode=''.join([x.upper() for x in self.spec])
            # dofv=[self.__dict__[p] for p in self.spec]
            specdict={p.upper():self.__dict__[p] for p in self.spec}
            retdict=self.suph.Bilinear(specdict)
            for p in self._p: 
                if p not in self.spec and p!='x':
                    self.__dict__[p]=retdict[p.upper()]
        else:
            ''' Th is below its liquid-state value; assume this is a subcooled state '''
            specdict={p.upper():self.__dict__[p] for p in self.spec}
            retdict=self.subc.Bilinear(specdict)
            for p in self._p: 
                if p not in self.spec and p!='x':
                    self.__dict__[p]=retdict[p.upper()]

    def _resolve_subsup(self):
        ''' T and P are both given explicitly.  Could be either superheated or subcooled state '''
        assert self.spec==['T','P']
        specdict={'T':self.T,'P':self.P}
        if self.satd.lim['T'][0]<self.T<self.satd.lim['T'][1]:
            Psat=self.satd.interpolators['T']['P'](self.T)
        else:
            Psat=LARGE
        if self.P>Psat:
            ''' P is higher than saturation: this is a subcooled state '''
            retdict=self.subc.Bilinear(specdict)
        else:
            ''' P is lower than saturation: this is a superheated state '''
            retdict=self.suph.Bilinear(specdict)
        for p in self._p: 
            if p not in self.spec and p!='x':
                self.__dict__[p]=retdict[p.upper()]

    def _resolve_satd(self):
        ''' This is an explicitly saturated state with vapor fraction (x) and one 
        other property (p) specified '''
        p=self.spec[0]
        self.Liquid=PHASE()
        self.Vapor=PHASE()
        if p=='T':
            ''' The other property is T; make sure it lies between saturation limits '''
            if self.T<self.satd.lim['T'][0] or self.T>self.satd.lim['T'][1]:
                raise Exception(f'Cannot have a saturated state at T = {self.T} C')
            ''' Assign all other property values by interpolation '''
            for q in self._sp:
                if q!='T':
                    prop=self.satd.interpolators['T'][q](self.T)
                    # print(q,prop)
                    if q=='P': self.__dict__[q]=prop
                    if q[-1]=='V':
                        self.Vapor.__dict__[q[0].lower()]=prop
                    elif q[-1]=='L':
                        self.Liquid.__dict__[q[0].lower()]=prop
            for q in self._p:
                if not q in 'PTx':
                    self.__dict__[q]=self.x*self.Vapor.__dict__[q]+(1-self.x)*self.Liquid.__dict__[q]
        elif p=='P':
            ''' The other property is P; make sure it lies between saturation limits '''
            if self.P<self.satd.lim['P'][0] or self.P>self.satd.lim['P'][1]:
                raise Exception(f'Cannot have a saturated state at P = {self.P} MPa')
            ''' Assign all other property values by interpolation '''
            for q in self._sp:
                if q!='P':
                    prop=self.satd.interpolators['P'][q](self.P)        
                    if q=='T': self.__dict__[q]=prop
                    if q[-1]=='V':
                        self.Vapor.__dict__[q[0].lower()]=prop
                    elif q[-1]=='L':
                        self.Liquid.__dict__[q[0].lower()]=prop
            for q in self._p:
                if not q in 'PTx':
                    self.__dict__[q]=self.x*self.Vapor.__dict__[q]+(1-self.x)*self.Liquid.__dict__[q]
        else:
            ''' The other property is neither T or P; must use a lever-rule-based interpolation '''
            self._resolve_satd_lever()

    def _resolve_satd_lever(self):
        p=self.spec[0]
        assert p!='T' and p!='P'
        ''' Vapor fraction and one other property value (not T or P) is given '''
        th=self.__dict__[p]
        x=self.__dict__['x']
        ''' Build an array of V-L mixed properties based on given value of x '''
        Y=np.array(self.satd.DF['T'][f'{p.upper()}V'])*x+np.array(self.satd.DF['T'][f'{p.upper()}L'])*(1-x)
        X=np.array(self.satd.DF['T']['T'])
        ''' define an interpolator '''
        f=interp1d(X,Y)
        try:
            ''' interpolate the Temperature '''
            self.T=f(x)
            ''' Assign all other property values '''
            for q in self._sp:
                if q!='T':
                    prop=self.satd.interpolators['T'][q](self.T)
                    if q=='P': self.__dict__[q]=prop
                    if q[-1]=='V':
                        self.Vapor.__dict__[q[0].lower()]=prop
                    elif q[-1]=='L':
                        self.Liquid.__dict__[q[0].lower()]=prop
            for q in self._p:
                if not q in 'PTx':
                    self.__dict__[q]=self.x*self.Vapor.__dict__[q]+(1-self.x)*self.Liquid.__dict__[q]
        except:
            raise Exception(f'Could not interpolate {p} = {th} at quality {x} from saturated steam table')
        
    def __init__(self,**kwargs):
        self.satd=_SATD
        self.suph=_SUPH
        self.subc=_SUBC
        for p in self._p:
            self.__dict__[p]=kwargs.get(p,None)
        self._resolve()

if __name__=='__main__':
    state1=SANDLER(T=525,P=10)
    print(state1.T,state1.P,state1.v,state1.u,state1.h,state1.s)   
    state2=SANDLER(T=525,v=0.009)
    print(state2.T,state2.P,state2.v,state2.u,state2.h,state2.s)
    state3=SANDLER(T=525,u=2904.9)
    print(state3.T,state3.P,state3.v,state3.u,state3.h,state3.s)
    state4=SANDLER(T=525,s=6.24)
    print(state4.T,state4.P,state4.v,state4.u,state4.h,state4.s)
    state5=SANDLER(T=525,h=3140.6)
    print(state5.T,state5.P,state5.v,state5.u,state5.h,state5.s)
    