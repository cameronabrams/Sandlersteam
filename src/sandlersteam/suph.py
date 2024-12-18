# Author: Cameron F. Abrams <cfa22@drexel.edu>

from importlib import resources
import os
import numpy as np
import pandas as pd
from .util import add_headers, my_split
from scipy.interpolate import LinearNDInterpolator

class SUPH:
    def __init__(self,phase='V'):
        with resources.path('sandlersteam','__init__.py') as f:
            inst_root=os.path.split(os.path.abspath(f.parent.parent))[0]
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
