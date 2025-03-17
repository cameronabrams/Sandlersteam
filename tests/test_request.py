from sandlersteam.request import Request
from sandlersteam.state import SteamTables as st
import pytest

def test_request_init():
    R=Request()
    assert(R.satdP==False)
    assert(R.satdT==False)
    assert(not any(R.suph))
    assert(len(R.suph)==0)

def test_request_make():
    R=Request()
    R.register('satdP')
    assert(R.satdP)
    R.register('satdT')
    assert(R.satdT)
    R.register(suphP=1.0)
    R.register(suphP=1.0)
    assert(R.suph[0]==1.0)
    R.register(suphP=2.0)
    assert(R.suph[1]==2.0)
    assert(len(R.suph)==2)
    R.register(suphP=-99.999)
    assert(len(R.suph)==2)
    R.register(subcP=50.0)
    assert(len(R.subc)==1)


def test_request_latex():
    R=Request()
    R.register(suphP=0.1)
    R.register('satdP')
    R.register(suphP=0.2)
    R.register(subcP=50.0)
    output=R.to_latex()
    with open('out.tex','w') as f:
        f.write(output)
    
