import concurrent.futures
import pandas as pd
import numpy as np
import clarite
import time
from pandas_genomics import sim, io, scalars, GenotypeDtype
from memory_profiler import profile  # Library to monitor memory consumption

# Define the variant for two SNPs
variant1 = scalars.Variant("1", 1, id="rs1", ref="A", alt=["C"])
variant2 = scalars.Variant("1", 2, id="rs2", ref="G", alt=["T"])

# Define Case-Control ratio
num_samples = 5000
n_cases = int(num_samples / 4)
n_controls = num_samples - n_cases

# Interations
n_loops = range(2)
work_group = (["REC", "ADD"], ["ADD", "ADD"], ["DOM", "ADD"], ["HED", "ADD"])

start = time.perf_counter()
EDGE_alpha_Final = pd.DataFrame()
All_Results_Final = pd.DataFrame()


# Name conversation for parameters
def conv_u(i):
    switcher = {"REC": "RECESSIVE", "ADD": "ADDITIVE", "DOM": "DOMINANT", "HED": "HET"}
    return switcher.get(i, "Invalid Group")


def conv_l(i):
    switcher = {
        "REC": "Recessive",
        "ADD": "Additive",
        "DOM": "Dominant",
        "HED": "Heterozygous",
    }
    return switcher.get(i, "Invalid Group")


# @profile
def do_something(seed):

    train_seed = seed
    test_seed = seed + 2000

    All_Results = pd.DataFrame()
    EDGE_alpha_Results = pd.DataFrame()

    for ab1, ab2 in work_group:

        ab1u = conv_u(ab1)
        ab2u = conv_u(ab2)
        ab1l = conv_l(ab1)

        ## Recessive Main Effect for SNP1 without interaction
        # Training data
        train_main = sim.BAMS.from_model(
            eff1=getattr(sim.SNPEffectEncodings, ab1u),
            eff2=getattr(sim.SNPEffectEncodings, ab2u),
            penetrance_base=0.465,
            main1=1,
            main2=0,
            interaction=0,
            snp1=variant1,
            snp2=variant2,
            random_seed=train_seed,
        )
        train_main = train_main.generate_case_control(
            n_cases=n_cases, n_controls=n_controls
        )

        # Calculate weights from the training dataset
        train_main = train_main.genomics.calculate_edge_encoding_values(
            data=train_main["Outcome"], outcome_variable="Outcome"
        )
        edge_weights_me = train_main.copy()
        edge_weights_me.insert(loc=0, column="BioAct", value=ab1l)
        edge_weights_me.insert(loc=0, column="TrainSeed", value=train_seed)
        edge_weights_me.insert(loc=0, column="TestSeed", value=test_seed)

        # Test data
        test_main = sim.BAMS.from_model(
            eff1=getattr(sim.SNPEffectEncodings, ab1u),
            eff2=getattr(sim.SNPEffectEncodings, ab2u),
            penetrance_base=0.465,
            main1=1,
            main2=0,
            interaction=0,
            snp1=variant1,
            snp2=variant2,
            random_seed=test_seed,
        )
        test_main = test_main.generate_case_control(
            n_cases=n_cases, n_controls=n_controls
        )
        test_main["Outcome"].cat.reorder_categories(["Control", "Case"], inplace=True)

        # Run Regression by using weightes from CLARITE
        # EDGE Encoding
        # test_me_CLARITE=test_main.genomics.encode_edge(encoding_info=edge_weights_me_t)
        edge_results_me = clarite.analyze.association_study(
            data=test_main.genomics.encode_edge(encoding_info=train_main),
            outcomes="Outcome",
        )
        edge_results_me["odds ratio"] = np.exp(edge_results_me["Beta"])
        edge_results_me.insert(loc=0, column="Encoding", value="EDGE")
        edge_results_me.insert(loc=0, column="BioAct", value=ab1l)
        edge_results_me.insert(loc=0, column="TrainSeed", value=train_seed)
        edge_results_me.insert(loc=0, column="TestSeed", value=test_seed)

        ##Concat the results
        EDGE_alpha_Results = pd.concat([EDGE_alpha_Results, edge_weights_me], axis=0)
        All_Results = pd.concat([All_Results, edge_results_me], axis=0)

        del train_main, edge_weights_me, edge_results_me, test_main

    return EDGE_alpha_Results, All_Results


# instantiating the decorator
# @profile
def xx():
    with concurrent.futures.ProcessPoolExecutor() as executor:

        # map(func, *iterables, timeout=None, chunksize=1)
        results = executor.map(do_something, n_loops)

        EDGE_alpha_Finals = pd.DataFrame()
        All_Results_Finals = pd.DataFrame()

        for EDGE_alpha_Results, All_Results in results:

            EDGE_alpha_Finals = pd.concat(
                [EDGE_alpha_Finals, EDGE_alpha_Results], axis=0
            )
            EDGE_alpha_Finals.to_csv("~/Public/teste1.txt", sep=";")

            All_Results_Finals = pd.concat([All_Results_Finals, All_Results], axis=0)
            All_Results_Finals.to_csv("~/Public/teste2.txt", sep=";")

    print(f"Finished in {round(time.perf_counter()-start,2)} second(s)")

    del All_Results_Finals, EDGE_alpha_Finals


if __name__ == "__main__":
    xx()
