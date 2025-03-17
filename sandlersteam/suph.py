# Author: Cameron F. Abrams <cfa22@drexel.edu>

import os
import numpy as np
import pandas as pd
from .util import add_headers, my_split
from scipy.interpolate import RegularGridInterpolator,LinearNDInterpolator
from . import data
__data_path__=os.path.dirname(data.__file__)

class SUPH:
    def __init__(self,phase='V'):
        data=__data_path__
        if phase=='V':
            ff='SandlerSuphSteamTables.txt'
        elif phase=='L':
            ff='SandlerSubcSteamTables.txt'
        fn=os.path.join(data,ff)
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
        DFS.append(ndf)
        self.data=pd.concat(DFS,ignore_index=True)
        dof=self.data.columns
        self.uniqs={}
        for d in dof:
            self.uniqs[d]=np.sort(np.array(list(set(self.data[d].to_list()))))
        
        # specify all interpolators for all possible pairs of state variables
        # these are not used since they don't do bilinear interpolation the way
        # a student would with physical steam tables
        # self.interpolators={}
        # for i in range(len(dof)):
        #     X=np.array(self.data[dof[i]].to_list())
        #     for j in range(i+1,len(dof)):
        #         Y=np.array(self.data[dof[j]].to_list())
        #         mainkey=''.join([dof[i],dof[j]])
        #         self.interpolators[mainkey]={}
        #         for k in range(len(dof)):
        #             if k != i and k != j:
        #                 Z=np.array(self.data[dof[k]].to_list())
        #                 self.interpolators[mainkey][dof[k]]=LinearNDInterpolator(list(zip(X,Y)),Z)

    # def GenBilinear(self,specdict):
    #     retdict={}
    #     # df=self.data
    #     dof=list(self.data.columns)
    #     mainkeyL=[]
    #     deps=[]
    #     for d in dof:
    #         if d in specdict.keys():
    #             mainkeyL.append(d)
    #             retdict[d]=specdict[d]
    #         else:
    #             deps.append(d)
    #     assert len(mainkeyL)==2,f'Must specify two valid degrees of freedom from {dof}; I found {specdict.keys()}'
    #     mainkey=''.join(mainkeyL)
    #     args=[specdict[k] for k in mainkeyL]
    #     for d in deps:
    #         retdict[d]=float(self.interpolators[mainkey][d](*args))
    #     return retdict

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
        tdf=df[df['P']==yi]
        if not tdf.empty:
            X=np.array(tdf['T'])
            for d in dof:
                if d not in 'TP':
                    Y=np.array(tdf[d])
                    retdict[d]=np.interp(xi,X,Y,left=np.nan,right=np.nan)
        else:
            for PL,PR in zip(self.uniqs['P'][:-1],self.uniqs['P'][1:]):
                if PL<yi<PR:
                    break
            else:
                raise Exception(f'P {yi} not between {self.uniqs["P"][0]} and {self.uniqs["P"][-1]} at {xi} C')
            X=np.array([PL,PR])
            LDF=df[df['P']==PL]
            RDF=df[df['P']==PR]
            LT=np.array(LDF['T'])
            RT=np.array(RDF['T'])
            CT=np.array([T for T in LT if T in RT])
            if xi in CT:
                for d in dof:
                    if d not in 'TP':
                        Y=np.array([LDF[LT==xi][d].values[0],RDF[RT==xi][d].values[0]])
                        retdict[d]=np.interp(yi,X,Y,left=np.nan,right=np.nan)
            else:
                for TL,TR in zip(CT[:-1],CT[1:]):
                    if TL<xi<TR:
                        break
                else:
                    raise Exception(f'T {xi} not between {CT[0]} and {CT[-1]} at {yi} MPa')
                LTDF=LDF[(LDF['T']==TL)|(LDF['T']==TR)].sort_values(by='T')
                RTDF=RDF[(RDF['T']==TL)|(RDF['T']==TR)].sort_values(by='T')
                iv=np.zeros(2)
                for p in dof:
                    if not p in 'TP':
                        for i in range(2):
                            Lp=LTDF[p].values[i]
                            Rp=RTDF[p].values[i]
                            Y=np.array([Lp,Rp])
                            iv[i]=np.interp(yi,X,Y,left=np.nan,right=np.nan)
                        retdict[p]=np.interp(xi,np.array([TL,TR]),iv)
        return retdict

    def to_latex(self,P):
        # generates latex version of a P-block of the superheated/subcooled steam table
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
            title=r'\noindent\begin{minipage}{0.6\textwidth}'+'\n'+r'\footnotesize\vspace{5mm}'+'\n'+r'\begin{center}'+'\n'+r'$P$ = '+f'{P}'+r' MPa\\*[1ex]'+'\n'
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
                    LLdat[d].append(np.interp(xi,X,Y,left=np.nan,right=np.nan))
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
                retdict[pp]=np.interp(yi,X,Y,left=np.nan,right=np.nan)
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
                    d[pp]=np.interp(yi,X,Y,left=np.nan,right=np.nan)
            X=np.array([PL,PR])
            for pp in dof:
                Y=np.array([ldict[pp],rdict[pp]])
                retdict[pp]=np.interp(xi,X,Y,left=np.nan,right=np.nan)
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
                        LLdat[d]=np.interp(xi,X,Y,left=np.nan,right=np.nan)
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
