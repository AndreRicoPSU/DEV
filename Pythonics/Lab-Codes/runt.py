import pandas as pd
import numpy as np
import time
start = time.time()
from pandas_genomics import sim, io, scalars, GenotypeDtype

# Define the variant for two SNPs
from pandas_genomics.scalars import Variant
variant1 = Variant('1', 1, id='rs1', ref='A', alt=['C'])
variant2 = Variant('1', 2, id='rs2', ref='G', alt=['T'])
print(variant1)
