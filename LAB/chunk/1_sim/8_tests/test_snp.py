from itertools import count
from pysnptools.snpreader import Bed
from pathlib import Path


snpreader = Bed("all")
print(snpreader)
print(snpreader.iid_count, snpreader.sid_count)
