import concurrent.futures
import os
import pandas as pd
import numpy as np
import clarite
import time
from pandas_genomics import sim, io, scalars, GenotypeDtype
from memory_profiler import profile


# Define the variant for two SNPs
variant1 = scalars.Variant("1", 1, id="rs1", ref="A", alt=["C"])
variant2 = scalars.Variant("1", 2, id="rs2", ref="G", alt=["T"])

# Define Case-Control ratio
num_samples = 200
case_control_ratio = "1:3"
n_cases = int(num_samples / 4)
n_controls = num_samples - n_cases
PEN_BASE = 0.05
PEN_DIFF = 0.25
MAFA = 0.05
MAFB = 0.05
SNR = 0.01

# Interations
n_loops = range(2)

work_group = (
    ["REC", "ADD"],
    ["SUB", "ADD"],
    ["ADD", "ADD"],
    ["SUP", "ADD"],
    ["DOM", "ADD"],
    ["HET", "ADD"])

# "DOMDEV", "Recessive","Dominant", "Codominant", "EDGE"
encoding = ("Additive", "DOMDEV", "Recessive",
            "Dominant", "Codominant", "EDGE")

start = time.perf_counter()
EDGE_alpha_Final = pd.DataFrame()
All_Results_Final = pd.DataFrame()


# Name conversation for parameters
def conv_u(i):
    switcher = {
        "REC": "RECESSIVE",
        "SUB": "SUB_ADDITIVE",
        "ADD": "ADDITIVE",
        "SUP": "SUPER_ADDITIVE",
        "DOM": "DOMINANT",
        "HET": "HET",
    }
    return switcher.get(i, "Invalid Group")


def conv_l(i):
    switcher = {
        "REC": "Recessive",
        "SUB": "Sub-Additive",
        "ADD": "Additive",
        "SUP": "Super_Additive",
        "DOM": "Dominant",
        "HET": "Heterozygous",
    }
    return switcher.get(i, "Invalid Group")


# @profile
def simulations(seed):

    train_seed = seed
    test_seed = seed + 2000

    all_results = pd.DataFrame()
    alpha_Results = pd.DataFrame()

    for ab1, ab2 in work_group:

        ab1u = conv_u(ab1)
        ab2u = conv_u(ab2)
        ab1l = conv_l(ab1)

        vpid = os.getpid()

        result = pd.DataFrame(
            [{"proc_number": vpid, "seed": train_seed, "ab1u": ab1u, "ab2u": ab2u}])

        # Run Regression by using weightes from CLARITE
        for v_enc in encoding:

            encode = str("encode_" + v_enc.lower())

            if v_enc != "DOMDEV":

                test_me_enc = pd.DataFrame([{"proc_number": vpid, "train": train_seed, "effec 2": ab1u, "effec 2": ab2u,
                                             "Encoding": encode, "if_st": "OTHER"}])
                all_results = pd.concat([all_results, test_me_enc], axis=0)
            else:
                test_me_enc = pd.DataFrame([{"proc_number": vpid, "train": train_seed, "effec 2": ab1u, "effec 2": ab2u,
                                             "Encoding": encode, "if_st": "DOM"}])
                all_results = pd.concat([all_results, test_me_enc], axis=0)

        alpha_Results = pd.concat([alpha_Results, result], axis=0)

        # breakpoint

    return alpha_Results, all_results


# instantiating the decorator
# @profile
def mp_treads():
    with concurrent.futures.ProcessPoolExecutor() as executor:

        results = executor.map(simulations, n_loops)

        alpha_Finals = pd.DataFrame()
        all_results_Finals = pd.DataFrame()

        for alpha_Results, All_Results in results:
            alpha_Finals = pd.concat([alpha_Finals, alpha_Results], axis=0)
            all_results_Finals = pd.concat(
                [all_results_Finals, All_Results], axis=0)

        alpha_Finals.to_csv("~/Public/alpha_final.txt", sep=";")
        all_results_Finals.to_csv("~/Public/all_results_final.txt", sep=";")

    print(f"Finished in {round(time.perf_counter()-start,2)} second(s)")


if __name__ == "__main__":
    mp_treads()
