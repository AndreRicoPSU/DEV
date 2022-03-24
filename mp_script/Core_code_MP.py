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
num_samples = 1000
case_control_ratio = "1:3"
n_cases = int(num_samples / 4)
n_controls = num_samples - n_cases
PEN_BASE = 0.05
PEN_DIFF = 0.25
MAFA = 0.05
MAFB = 0.05
SNR = 0.01

# Interations
n_loops = range(1)

work_group = (
    ["REC", "ADD"],
    ["SUB", "ADD"],
    ["ADD", "ADD"],
    ["SUP", "ADD"],
    ["DOM", "ADD"],
    ["HET", "ADD"],
    ["NUL", "ADD"],)

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
        "NUL": "ADDITIVE",
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
        "NUL": "NULL",
    }
    return switcher.get(i, "Invalid Group")


# @profile
for train_seed, test_seed in zip(range(0, 1), range(2000, 2001)):

    ALL_RESULTS_ENCODING = pd.DataFrame()
    ALL_RESULTS_ENCODING_EDGE = pd.DataFrame()
    ALL_RESULTS_EDGE_ALPHA = pd.DataFrame()

    for ab1, ab2 in work_group:

        ab1u = conv_u(ab1)
        ab2u = conv_u(ab2)
        ab1l = conv_l(ab1)

        # BLOCK 1
        # Recessive Main Effect for SNP1 without interaction
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
            random_seed=train_seed)

        # BLOCK 2
        train_me = train_main.generate_case_control(
            n_cases=n_cases, n_controls=n_controls, maf1=MAFA, maf2=MAFB, snr=SNR)

        # BLOCK 3
        # Calculate weights from the training dataset
        edge_weights_me_t = train_me.genomics.calculate_edge_encoding_values(
            data=train_me["Outcome"], outcome_variable="Outcome")
        edge_weights_me = edge_weights_me_t.copy()
        edge_weights_me.insert(loc=0, column="BioAct", value=ab1l)
        edge_weights_me.insert(loc=0, column="TrainSeed", value=train_seed)
        edge_weights_me.insert(loc=0, column="TestSeed", value=test_seed)

        # BLOCK 4
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
            random_seed=test_seed)
        test_me = test_main.generate_case_control(
            n_cases=n_cases, n_controls=n_controls, maf1=MAFA, maf2=MAFB, snr=SNR)
        # test_me = test_me["Outcome"].cat.reorder_categories(["Control", "Case"]) # didnt works
        test_me["Outcome"].cat.reorder_categories(
            ["Control", "Case"], inplace=True)  # Python will stop to support this inplace

        # Run Regression by using weightes from CLARITE
        for v_enc in encoding:
            #encode = str("encode_" + v_enc.lower())
            # test_me_enc = getattr(test_me.genomics, encode) # More Pythonics but didnt works
            if v_enc != "DOMDEV":
                if v_enc.lower() == "additive":
                    test_me_enc = test_me.genomics.encode_additive()
                elif v_enc.lower() == "recessive":
                    test_me_enc = test_me.genomics.encode_recessive()
                elif v_enc.lower() == "dominant":
                    test_me_enc = test_me.genomics.encode_dominant()
                elif v_enc.lower() == "codominant":
                    test_me_enc = test_me.genomics.encode_codominant()
                else:
                    test_me_enc = test_me.genomics.encode_edge(
                        encoding_info=edge_weights_me_t)

                results_me = clarite.analyze.association_study(
                    data=test_me_enc, outcomes="Outcome"
                )
                results_me["odds ratio"] = np.exp(results_me["Beta"])
                results_me.insert(loc=0, column="Encoding", value=v_enc)
                results_me.insert(loc=0, column="BioAct", value=ab1l)
                results_me.insert(loc=0, column="TrainSeed", value=train_seed)
                results_me.insert(loc=0, column="TestSeed", value=test_seed)

                ALL_RESULTS_ENCODING = pd.concat(
                    [ALL_RESULTS_ENCODING, results_me], axis=0)

                if v_enc.lower() == "edge":
                    ALL_RESULTS_ENCOGIND_EDGE = pd.concat(
                        [ALL_RESULTS_ENCODING_EDGE, results_me], axis=0)

            else:
                # DOMDEV Encoding
                test_me_pb000_DOMDEV = test_me_enc
                test_me_pb000_DOMDEV["COV1"] = test_me_pb000_DOMDEV["SNP1"]
                test_me_pb000_DOMDEV["COV2"] = test_me_pb000_DOMDEV["SNP2"]

                test_me_pb000_DOMDEV["SNP1"] = test_me_pb000_DOMDEV["SNP1"].replace(
                    2, 0)
                test_me_pb000_DOMDEV["SNP2"] = test_me_pb000_DOMDEV["SNP2"].replace(
                    2, 0)

                test_me_pb000_DOMDEV_SNP1_t = test_me_pb000_DOMDEV[[
                    "Outcome", "SNP1", "COV1"]]
                test_me_pb000_DOMDEV_SNP2_t = test_me_pb000_DOMDEV[[
                    "Outcome", "SNP2", "COV2"]]

                DOMDEV_results_me_pb000_SNP1 = clarite.analyze.association_study(
                    data=test_me_pb000_DOMDEV_SNP1_t,
                    outcomes="Outcome",
                    covariates=["COV1"],
                )
                DOMDEV_results_me_pb000_SNP2 = clarite.analyze.association_study(
                    data=test_me_pb000_DOMDEV_SNP2_t,
                    outcomes="Outcome",
                    covariates=["COV2"],
                )

                DOMDEV_results_me = pd.concat(
                    [DOMDEV_results_me_pb000_SNP1, DOMDEV_results_me_pb000_SNP2])

                DOMDEV_results_me["odds ratio"] = np.exp(
                    DOMDEV_results_me["Beta"])
                DOMDEV_results_me.insert(loc=0, column="Encoding", value=v_enc)
                DOMDEV_results_me.insert(loc=0, column="BioAct", value=ab1l)
                DOMDEV_results_me.insert(
                    loc=0, column="TrainSeed", value=train_seed)
                DOMDEV_results_me.insert(
                    loc=0, column="TestSeed", value=test_seed)

                ALL_RESULTS_ENCODING = pd.concat(
                    [ALL_RESULTS_ENCODING, DOMDEV_results_me], axis=0)

        ALL_RESULTS_EDGE_ALPHA = pd.concat(
            [ALL_RESULTS_EDGE_ALPHA, edge_weights_me], axis=0)

    #ALL_ALPHA_RESULTS = pd.concat([ALL_ALPHA_RESULTS, train_me], axis=0)

print("finish")
