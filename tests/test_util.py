import os
from importlib import resources
from sandlersteam.util import data_path

def test_data_location():
    with resources.as_file(resources.files('sandlersteam')) as src_root:
        dp=data_path()
        assert dp==os.path.join(src_root,'data')