from sandlersteam.state import RandomSample
from sandlersteam.state import SteamTables

def test_random_sample():
    R=RandomSample(phase='suph')
    assert R.x == None
    R=RandomSample(phase='satd')
    assert R.x == 1.0
    R=RandomSample(phase='satd',satdDOF='P')
    assert R.x == 1.0