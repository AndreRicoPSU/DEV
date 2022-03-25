import math
import pandas as pd
import numpy as np
import os
import sys
import time
import psutil
import pandas as pd
import numpy as np
import concurrent.futures
import datetime


from pathlib import Path
from typing import Optional, Union
from pandas_genomics.scalars import MISSING_IDX


def split_dataframes_by_positions(df_variant, v_chunk_split, idx):

    hdf_s = 
    idx_s = v_chunk_split * idx
    idx_f = v_chunk_split * (idx + 1)
    hdf.put(str("/bim/variants_" + str(idx)), df_variant.iloc[idx_s:idx_f, :])
    hdf.flush()
    print(f"position {idx} will start in {idx_s} and finish in {idx_f}")
    return idx_s, idx_f


"""
    if v_cicles > 1:
        for i in range(0, v_cicles):
            vs = ve
            ve = ve + v_chunk

            hdf.put(str("/bim/variants_" + str(i)), variant_list.iloc[vs:ve, :])

            print(
                f"\tSaved Genetypes Dataset {i} with {len(variant_list.iloc[vs:ve, :])} variants"
            )
            hdf.flush()

    else:
        hdf.put(str("/bim/variants_0"), variant_list)
        print(f"\tSaved Genetypes single file")
        hdf.flush()

"""


def from_plink(input: Union[str, Path], v_chunk: Optional[int] = None):

    # All three files are used
    input = str(input)  # coerce to string in order to add extensions
    bim_file = Path(input + ".bim")
    fam_file = Path(input + ".fam")

    # Make sure each file exists
    if not Path(bim_file).exists():
        raise ValueError(f"The .bim file was not found\n\t{str(bim_file)}")
    if not Path(fam_file).exists():
        raise ValueError(f"The .fam file was not found\n\t{str(fam_file)}")

    print(f"Loading genetic data from '{fam_file.stem}'")

    # hdf = pd.HDFStore("/storage/home/alr6366/scratch/data.h5")  # Open .h5 files as Append Mode

    hdf = pd.HDFStore("~/public/data.h5")  # Open .h5 files as Append Mode

    # Load fam file (PLINK sample information file) into a list of Samples
    print("----------Started Samples Process---------")
    df = pd.read_table(fam_file, header=None, sep=" ")
    df.columns = ["FID", "IID", "IID_father", "IID_mother", "sex", "phenotype"]
    hdf.put(str("/fam/samples"), df)
    hdf.flush()
    print(f"\tLoaded information for {len(df)} samples from '{fam_file.name}'")
    print(f"\tSaved in /fam/samples table in ....")

    # Load bim file (PLINK extended MAP file) into a list of variants
    print("----------Started Variants Process---------")
    variant_list = pd.read_table(bim_file, header=None, sep="\t")
    v_cicles = math.ceil(len(variant_list) / v_chunk)

    hdf.put("/cicles", pd.Series(v_cicles))
    hdf.close()

    df_results = []

    num_procs = psutil.cpu_count(logical=False)
    if len(sys.argv) > 1:
        num_procs = int(sys.argv[1])

    start = time.time()
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_procs) as executor:
        results = [
            executor.submit(
                split_dataframes_by_positions, variant_list, v_chunk, idx=x
            )
            for x in range(v_cicles)
        ]
        for result in concurrent.futures.as_completed(results):
            try:
                df_results.append(result.result())
            except Exception as ex:
                print(str(ex))
                pass

    end = time.time()

    return "OK"


if __name__ == "__main__":
    DATA_DIR = Path(__file__).parent / "data" / "plink"
    input = DATA_DIR / "plink_test_medium"
    result = from_plink(input, 1000)
    print("Process finished with status: ", result)
