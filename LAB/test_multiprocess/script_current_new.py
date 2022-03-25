import concurrent.futures
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
num_samples = 5000
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
    ["HET", "ADD"],
)
encoding = ("Additive", "DOMDEV", "Recessive", "Dominant", "Codominant", "EDGE")

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
def do_something(seed):

    train_seed = seed
    test_seed = seed + 2000

    All_Results = pd.DataFrame()
    alpha_Results = pd.DataFrame()

    for ab1, ab2 in work_group:

        ab1u = conv_u(ab1)
        ab2u = conv_u(ab2)
        ab1l = conv_l(ab1)

        ## Recessive Main Effect for SNP1 without interaction
        # Training data
        train_main = sim.BAMS.from_model(
            eff1=getattr(sim.SNPEffectEncodings, ab1u),
            eff2=getattr(sim.SNPEffectEncodings, ab2u),
            penetrance_base=PEN_BASE,
            penetrance_diff=PEN_DIFF,
            main1=1,
            main2=0,
            interaction=0,
            snp1=variant1,
            snp2=variant2,
            random_seed=train_seed,
        )
        train_me = train_main.generate_case_control(
            n_cases=n_cases, n_controls=n_controls, maf1=MAFA, maf2=MAFB, snr=SNR
        )

        # Calculate weights from the training dataset
        weights_me_t = train_me.genomics.calculate_edge_encoding_values(
            data=train_me["Outcome"], outcome_variable="Outcome"
        )
        weights_me = weights_me_t.copy()
        weights_me.insert(loc=0, column="BioAct", value=ab1l)
        weights_me.insert(loc=0, column="TrainSeed", value=train_seed)
        weights_me.insert(loc=0, column="TestSeed", value=test_seed)

        # Test data
        test_main = sim.BAMS.from_model(
            eff1=getattr(sim.SNPEffectEncodings, ab1u),
            eff2=getattr(sim.SNPEffectEncodings, ab2u),
            penetrance_base=PEN_BASE,
            penetrance_diff=PEN_DIFF,
            main1=1,
            main2=0,
            interaction=0,
            snp1=variant1,
            snp2=variant2,
            random_seed=test_seed,
        )
        test_me = test_main.generate_case_control(
            n_cases=n_cases, n_controls=n_controls, maf1=MAFA, maf2=MAFB, snr=SNR
        )
        test_me["Outcome"].cat.reorder_categories(["Control", "Case"], inplace=True)
        # test_me = test_me["Outcome"].cat.reorder_categories(["Control", "Case"])

        # Run Regression by using weightes from CLARITE
        for v_enc in encoding:

            encode = str("encode_" + v_enc.lower())

            if v_enc != "DOMDEV":
                test_me_enc = test_me.genomics.encode_additive()
                # test_me_enc = getattr(test_me.genomics, encode)

                results_me = clarite.analyze.association_study(
                    data=test_me_enc, outcomes="Outcome"
                )
                results_me["odds ratio"] = np.exp(results_me["Beta"])
                results_me.insert(loc=0, column="Encoding", value=v_enc)
                results_me.insert(loc=0, column="BioAct", value=ab1l)
                results_me.insert(loc=0, column="TrainSeed", value=train_seed)
                results_me.insert(loc=0, column="TestSeed", value=test_seed)

                All_Results = pd.concat([All_Results, results_me], axis=0)
                # All_Results = pd.concat([All_Results, test_me_enc], axis=0)
            else:
                # DOMDEV Encoding
                test_me_pb000_DOMDEV = test_me
                test_me_pb000_DOMDEV["COV1"] = test_me_pb000_DOMDEV["SNP1"]
                test_me_pb000_DOMDEV["COV2"] = test_me_pb000_DOMDEV["SNP2"]

                test_me_pb000_DOMDEV["SNP1"] = test_me_pb000_DOMDEV["SNP1"].replace(
                    2, 0
                )
                test_me_pb000_DOMDEV["SNP2"] = test_me_pb000_DOMDEV["SNP2"].replace(
                    2, 0
                )

                test_me_pb000_DOMDEV_SNP1 = test_me_pb000_DOMDEV[
                    ["Outcome", "SNP1", "COV1"]
                ]
                test_me_pb000_DOMDEV_SNP2 = test_me_pb000_DOMDEV[
                    ["Outcome", "SNP2", "COV2"]
                ]
                """
                DOMDEV_results_me_pb000_SNP1 = clarite.analyze.association_study(
                    data=test_me_pb000_DOMDEV_SNP1,
                    outcomes="Outcome",
                    covariates=["COV1"],
                )
                DOMDEV_results_me_pb000_SNP2 = clarite.analyze.association_study(
                    data=test_me_pb000_DOMDEV_SNP2,
                    outcomes="Outcome",
                    covariates=["COV2"],
                )
                DOMDEV_results_me = pd.concat(
                    [DOMDEV_results_me_pb000_SNP1, DOMDEV_results_me_pb000_SNP2]
                )

                DOMDEV_results_me["odds ratio"] = np.exp(DOMDEV_results_me["Beta"])
                DOMDEV_results_me.insert(loc=0, column="Encoding", value="DOMDEV")
                DOMDEV_results_me.insert(loc=0, column="BioAct", value=v_enc)
                DOMDEV_results_me.insert(loc=0, column="TrainSeed", value=train_seed)
                DOMDEV_results_me.insert(loc=0, column="TestSeed", value=test_seed)

                All_Results = pd.concat([All_Results, DOMDEV_results_me], axis=0)"""
                All_Results = pd.concat([All_Results, test_me_pb000_DOMDEV], axis=0)

        ##Concat the results
        alpha_Results = pd.concat([alpha_Results, weights_me], axis=0)

    return alpha_Results, All_Results


# instantiating the decorator
# @profile
def mp_treads():
    with concurrent.futures.ProcessPoolExecutor() as executor:

        results = executor.map(do_something, n_loops)

        alpha_Finals = pd.DataFrame()
        All_Results_Finals = pd.DataFrame()

        for alpha_Results, All_Results in results:

            alpha_Finals = pd.concat([alpha_Finals, alpha_Results], axis=0)
            All_Results_Finals = pd.concat([All_Results_Finals, All_Results], axis=0)

        alpha_Finals.to_csv("~/Public/teste1.txt", sep=";")
        All_Results_Finals.to_csv("~/Public/teste2.txt", sep=";")

    print(f"Finished in {round(time.perf_counter()-start,2)} second(s)")


if __name__ == "__main__":
    mp_treads()
