from sandlersteam.state import RandomSample
from sandlersteam.state import SteamTables
import pytest

def test_random_sample():
    R=RandomSample(phase='suph')
    assert R.x == None
    R=RandomSample(phase='satd')
    assert R.x == 1.0
    assert hasattr(R,'Vapor')
    assert hasattr(R,'Liquid')
    R=RandomSample(phase='satd',satdDOF='P')
    assert R.x == 1.0
    assert hasattr(R,'Vapor')
    assert hasattr(R,'Liquid')
    R1=RandomSample(phase='satd',satdDOF='T',seed=12345)
    R2=RandomSample(phase='satd',satdDOF='T',seed=54321)
    # fingers crossed
    assert R1.T != R2.T
    R3=RandomSample(phase='satd',satdDOF='T',seed=12345)
    assert R1.T == R3.T

def test_random_in_ranges():
    R=RandomSample(phase='suph',Trange=[500,800])
    assert R.T >= 500
    assert R.T <= 800
    R=RandomSample(phase='suph',Prange=[1,10])
    assert R.P >= 1
    assert R.P <= 10
    R=RandomSample(phase='suph',Trange=[400,500],Prange=[1,10])
    assert R.T >= 400
    assert R.T <= 500
    assert R.P >= 1
    assert R.P <= 10
    R=RandomSample(phase='satd',satdDOF='P',Prange=[1,10])
    assert R.P >= 1
    assert R.P <= 10
    R=RandomSample(phase='satd',satdDOF='T',Trange=[100,200])
    assert R.T >= 100
    assert R.T <= 200
    with pytest.raises(ValueError) as excinfo:
        R=RandomSample(phase='satd',satdDOF='T',Trange=[500,800])
    assert excinfo is not None
    with pytest.raises(ValueError) as excinfo:
        R=RandomSample(phase='satd',satdDOF='P',Trange=[500,800])
    assert excinfo is not None
