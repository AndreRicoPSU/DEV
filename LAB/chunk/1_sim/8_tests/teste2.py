import os
import sys
import time
import psutil
import pandas as pd
import numpy as np
import concurrent.futures
import datetime


def do_something(df):
    """
    do stuff
    """
    pid = os.getpid()
    ppid = os.getppid()
    start = time.time()
    print("PPID %s->%s Started" % (ppid, pid))
    df["diff"] = datetime.datetime.now() - pd.to_datetime(df["Order Date"])
    stop = time.time()
    completed_in = round(stop - start, 2)
    return df


if __name__ == "__main__":

    logical = False
    df_results = []
    num_procs = psutil.cpu_count(logical=logical)
    if len(sys.argv) > 1:
        num_procs = int(sys.argv[1])

    big_dataframe = pd.read_csv("5m.csv")
    splitted_df = np.array_split(big_dataframe, num_procs)
    start = time.time()

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
