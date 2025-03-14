# Author: Cameron F. Abrams, cfa22@drexel.edu
from .state import SteamTables as st

class Request:
    def __init__(self):
        self.suph=[]
        self.subc=[]
        self.satdP=False
        self.satdT=False
    def register(self,*args,**kwargs):
        if 'satdP' in args:
            self.satdP=True
        if 'satdT' in args:
            self.satdT=True
        if 'suphP' in kwargs:
            P=kwargs['suphP']
            if P in st['suph'].uniqs['P'] and not P in self.suph:
                self.suph.append(P)
        if 'subcP' in kwargs:
            P=kwargs['subcP']
            if P in st['subc'].uniqs['P'] and not P in self.subc:
                self.subc.append(P)

    def to_latex(self):
        unit_string=r"""\noindent $\hat{V}\ [=]\ \mbox{m$^3$/kg}$; $\hat{U}\ [=]\ \mbox{kJ/kg}$; $\hat{H}\ [=]\ \mbox{kJ/kg}$; $\hat{S}\ [=]\ \mbox{kJ/kg-K}$"""
        tables=[]
        if any(self.suph) or self.satdP or self.satdP:
            tables.append(r"""
\clearpage
\noindent THERMODYNAMIC PROPERTIES OF STEAM (Selected)\\*[1cm]""")
        if any(self.suph):
            tables.append(r"""Superheated steam:\\*[0mm]""")
        for p in sorted(self.suph):
            tables.append(st['suph'].to_latex(P=p))
        if any(self.suph):
            tables.append(unit_string+r'\\*[1cm]')

        if any(self.subc):
            tables.append(r"""Subcooled liquid:\\*[1cm]""")
        for p in sorted(self.subc):
            tables.append(st['subc'].to_latex(P=p))
        
        if self.satdP or self.satdT:
            if len(self.suph)+len(self.subc)>1:
                tables.append(r"""\clearpage""")
            tables.append(r"""Saturated steam:\\*[5mm]""")
            if self.satdP:
                tables.append(st['satd'].to_latex(by='P'))
                tables.append(unit_string)
            if self.satdT:
                if self.satdP:
                    tables.append(r"""\clearpage""")
                tables.append(st['satd'].to_latex(by='T'))
                tables.append(unit_string)
        
        return '\n'.join(tables)
        
