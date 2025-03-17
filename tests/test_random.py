from sandlersteam.state import RandomSample
from sandlersteam.state import SteamTables

def test_random_sample():
    R=RandomSample(phase='suph')
    assert R.x == None
    R=RandomSample(phase='satd')
    assert R.x == 1.0
    R=RandomSample(phase='satd',satdDOF='P')
    assert R.x == 1.0
    R1=RandomSample(phase='satd',satdDOF='T',seed=12345)
    R2=RandomSample(phase='satd',satdDOF='T',seed=54321)
    # fingers crossed
    assert R1.T != R2.T
    R3=RandomSample(phase='satd',satdDOF='T',seed=12345)
    assert R1.T == R3.T