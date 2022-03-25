from gc import callbacks
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
import multiprocessing as mp
import datetime


from pathlib import Path
from typing import Optional, Union
from pandas_genomics.scalars import MISSING_IDX

"""
def split_dataframes_by_positions(df_variant, v_chunk_split, idx):

    idx_s = v_chunk_split * idx
    idx_f = v_chunk_split * (idx + 1)
    result = df_variant.iloc[idx_s:idx_f, :]
    print(f"position {idx} will start in {idx_s} and finish in {idx_f}")
    return result, idx
"""


def handle_output(result, idx):
    dataset = str("/bim/variants_" + str(idx(0)))
    hdf = pd.HDFStore("~/public/data.h5", mode="a")
    hdf.put(str("/bim/variants_" + str(idx(0))), result)

    hdf.flush()
    hdf.close()
    return dataset


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

    variant_splitted = np.array_split(variant_list, v_cicles)
    variant_splitted = np.asarray(variant_splitted, dtype=object)

    num_procs = psutil.cpu_count(logical=False)
    if len(sys.argv) > 1:
        num_procs = int(sys.argv[1])

    num_procs = 1
    df_results = []

    for idx, x in np.ndenumerate(variant_splitted):
        dataset = str("/bim/variants_" + str(idx[0]))
        print(dataset)
        hdf = pd.HDFStore("~/public/data.h5", mode="a")
        hdf.put(dataset, x)

        hdf.flush()
        hdf.close()

    # with concurrent.futures.ProcessPoolExecutor(max_workers=num_procs) as executor:

    #    results = executor.map(handle_output, variant_splitted)

    # for result in results:
    #    print(result)

    # results = [
    #    executor.submit(handle_output, df=df, idx=idx)
    #    for idx, df in np.ndenumerate(variant_splitted)
    # ]
    # for result in concurrent.futures.as_completed(results):
    #    try:
    #        df_results.append(result.result())
    #    except Exception as ex:
    #        print(str(ex))
    #        pass

    """
    #esta acessando apenas a primeira linha da funcao.
    pool = mp.Pool(num_procs)
    for idx, x in np.ndenumerate(variant_splitted):
        pool.apply_async(handle_output, (x, idx))
    pool.close()
    pool.join()
    
    """

    """
    Antiga com chunk manual
    pool = mp.Pool(num_procs)
    for i in range(v_cicles):
        pool.apply_async(
            split_dataframes_by_positions,
            (variant_list, v_chunk, i),
            callback=handle_output(result, i),
        )
    pool.close()
    pool.join()
    """

    end = time.time()

    return "OK"
    """
    with concurrent.futures.ProcessPoolExecutor(max_workers=num_procs) as executor:
        results = [executor.submit(do_something, df=df) for df in splitted_df]
        for result in concurrent.futures.as_completed(results):
            try:
                df_results.append(result.result())
            except Exception as ex:
                print(str(ex))
                pass
    end = time.time()
    print("-------------------------------------------")
    print("PPID %s Completed in %s" % (os.getpid(), round(end - start, 2)))
    df_results = pd.concat(df_results)
    """


if __name__ == "__main__":
    DATA_DIR = Path(__file__).parent / "data" / "plink"
    input = DATA_DIR / "plink_test_medium"
    result = from_plink(input, 1000)
    print("Process finished with status: ", result)
