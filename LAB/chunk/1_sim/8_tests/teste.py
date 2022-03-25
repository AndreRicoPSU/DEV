import time
import tables
import numpy as np
import pandas as pd

df = pd.DataFrame(
    {"A": np.random.randn(9)},
    index=pd.MultiIndex.from_product(
        [range(3), list("abc")], names=["first", "second"]
    ),
)

print(df)

hdf = pd.HDFStore("~/public/Genotype.h5")  # Open .h5 files as Append Mode
hdf.put(str("/fam/samples"), df, "t")
