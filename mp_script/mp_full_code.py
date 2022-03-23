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
num_samples = 2000
case_control_ratio = "1:3"
n_cases = int(num_samples / 4)
n_controls = num_samples - n_cases
PEN_BASE = 0.05
PEN_DIFF = 0.25
MAFA = 0.05
MAFB = 0.05
SNR = 0.01

# Interations
n_loops = range(10)


start = time.perf_counter()
EDGE_alpha_Final = pd.DataFrame()
All_Results_Final = pd.DataFrame()


# @profile
def simulations(seed):

    train_seed = seed
    test_seed = seed + 2000

    startcycle = time.time()
    # Recessive Main Effect for SNP1 without interaction
    # Training data
    train_rec_main_effect00 = sim.BAMS.from_model(
        eff1=sim.SNPEffectEncodings.RECESSIVE,
        eff2=sim.SNPEffectEncodings.ADDITIVE,
        penetrance_base=PEN_BASE,
        penetrance_diff=PEN_DIFF,
        main1=1,
        main2=0,
        interaction=0,
        snp1=variant1,
        snp2=variant2,
        random_seed=train_seed,
    )
    train_rec_me_pb000 = train_rec_main_effect00.generate_case_control(
        n_cases=n_cases, n_controls=n_controls, maf1=MAFA, maf2=MAFB, snr=SNR
    )

    # Calculate weights from the training dataset
    edge_weights_rec_me_pb000 = (
        train_rec_me_pb000.genomics.calculate_edge_encoding_values(
            data=train_rec_me_pb000["Outcome"], outcome_variable="Outcome"
        )
    )
    edge_weights_rec_me = edge_weights_rec_me_pb000.copy()
    edge_weights_rec_me.insert(loc=0, column="BioAct", value="Recessive")
    edge_weights_rec_me.insert(loc=0, column="TrainSeed", value=train_seed)
    edge_weights_rec_me.insert(loc=0, column="TestSeed", value=test_seed)

    # Test data
    test_rec_main_effect00 = sim.BAMS.from_model(
        eff1=sim.SNPEffectEncodings.RECESSIVE,
        eff2=sim.SNPEffectEncodings.ADDITIVE,
        penetrance_base=PEN_BASE,
        penetrance_diff=PEN_DIFF,
        main1=1,
        main2=0,
        interaction=0,
        snp1=variant1,
        snp2=variant2,
        random_seed=test_seed,
    )
    test_rec_me_pb000 = test_rec_main_effect00.generate_case_control(
        n_cases=n_cases, n_controls=n_controls, maf1=MAFA, maf2=MAFB, snr=SNR
    )
    test_rec_me_pb000["Outcome"].cat.reorder_categories(
        ["Control", "Case"], inplace=True
    )

    # Run Regression by using weightes from CLARITE
    # Addtive Encoding
    test_rec_me_pb000_ADD = test_rec_me_pb000.genomics.encode_additive()
    add_results_rec_me_pb000 = clarite.analyze.association_study(
        data=test_rec_me_pb000_ADD, outcomes="Outcome"
    )
    add_results_rec_me_pb000["odds ratio"] = np.exp(
        add_results_rec_me_pb000["Beta"])
    add_results_rec_me_pb000.insert(loc=0, column="Encoding", value="Additive")
    add_results_rec_me_pb000.insert(loc=0, column="BioAct", value="Recessive")
    add_results_rec_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    add_results_rec_me_pb000.insert(loc=0, column="TestSeed", value=test_seed)

    # DOMDEV Encoding
    test_rec_me_pb000_DOMDEV = test_rec_me_pb000_ADD
    test_rec_me_pb000_DOMDEV["COV1"] = test_rec_me_pb000_DOMDEV["SNP1"]
    test_rec_me_pb000_DOMDEV["COV2"] = test_rec_me_pb000_DOMDEV["SNP2"]

    test_rec_me_pb000_DOMDEV["SNP1"] = test_rec_me_pb000_DOMDEV["SNP1"].replace(
        2, 0)
    test_rec_me_pb000_DOMDEV["SNP2"] = test_rec_me_pb000_DOMDEV["SNP2"].replace(
        2, 0)

    test_rec_me_pb000_DOMDEV_SNP1 = test_rec_me_pb000_DOMDEV[
        ["Outcome", "SNP1", "COV1"]
    ]
    test_rec_me_pb000_DOMDEV_SNP2 = test_rec_me_pb000_DOMDEV[
        ["Outcome", "SNP2", "COV2"]
    ]
    DOMDEV_results_rec_me_pb000_SNP1 = clarite.analyze.association_study(
        data=test_rec_me_pb000_DOMDEV_SNP1, outcomes="Outcome", covariates=["COV1"]
    )
    DOMDEV_results_rec_me_pb000_SNP2 = clarite.analyze.association_study(
        data=test_rec_me_pb000_DOMDEV_SNP2, outcomes="Outcome", covariates=["COV2"]
    )
    DOMDEV_results_rec_me_pb000 = pd.concat(
        [DOMDEV_results_rec_me_pb000_SNP1, DOMDEV_results_rec_me_pb000_SNP2]
    )

    DOMDEV_results_rec_me_pb000["odds ratio"] = np.exp(
        DOMDEV_results_rec_me_pb000["Beta"]
    )
    DOMDEV_results_rec_me_pb000.insert(
        loc=0, column="Encoding", value="DOMDEV")
    DOMDEV_results_rec_me_pb000.insert(
        loc=0, column="BioAct", value="Recessive")
    DOMDEV_results_rec_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    DOMDEV_results_rec_me_pb000.insert(
        loc=0, column="TestSeed", value=test_seed)

    # Recessive Encoding
    test_rec_me_pb000_REC = test_rec_me_pb000.genomics.encode_recessive()
    rec_results_rec_me_pb000 = clarite.analyze.association_study(
        data=test_rec_me_pb000_REC, outcomes="Outcome"
    )
    rec_results_rec_me_pb000["odds ratio"] = np.exp(
        rec_results_rec_me_pb000["Beta"])
    rec_results_rec_me_pb000.insert(
        loc=0, column="Encoding", value="Recessive")
    rec_results_rec_me_pb000.insert(loc=0, column="BioAct", value="Recessive")
    rec_results_rec_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    rec_results_rec_me_pb000.insert(loc=0, column="TestSeed", value=test_seed)

    # Dominant Encoding
    test_rec_me_pb000_DOM = test_rec_me_pb000.genomics.encode_dominant()
    dom_results_rec_me_pb000 = clarite.analyze.association_study(
        data=test_rec_me_pb000_DOM, outcomes="Outcome"
    )
    dom_results_rec_me_pb000["odds ratio"] = np.exp(
        dom_results_rec_me_pb000["Beta"])
    dom_results_rec_me_pb000.insert(loc=0, column="Encoding", value="Dominant")
    dom_results_rec_me_pb000.insert(loc=0, column="BioAct", value="Recessive")
    dom_results_rec_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    dom_results_rec_me_pb000.insert(loc=0, column="TestSeed", value=test_seed)

    # Codominant Encoding
    test_rec_me_pb000_CODOM = test_rec_me_pb000.genomics.encode_codominant()
    codom_results_rec_me_pb000 = clarite.analyze.association_study(
        data=test_rec_me_pb000_CODOM, outcomes="Outcome"
    )
    codom_results_rec_me_pb000["odds ratio"] = np.exp(
        codom_results_rec_me_pb000["Beta"]
    )
    codom_results_rec_me_pb000.insert(
        loc=0, column="Encoding", value="Codominant")
    codom_results_rec_me_pb000.insert(
        loc=0, column="BioAct", value="Recessive")
    codom_results_rec_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    codom_results_rec_me_pb000.insert(
        loc=0, column="TestSeed", value=test_seed)

    # EDGE Encoding
    test_rec_me_pb000_CLARITE = test_rec_me_pb000.genomics.encode_edge(
        encoding_info=edge_weights_rec_me_pb000
    )

    edge_results_rec_me_pb000 = clarite.analyze.association_study(
        data=test_rec_me_pb000_CLARITE, outcomes="Outcome"
    )
    edge_results_rec_me_pb000["odds ratio"] = np.exp(
        edge_results_rec_me_pb000["Beta"])
    edge_results_rec_me_pb000.insert(loc=0, column="Encoding", value="EDGE")
    edge_results_rec_me_pb000.insert(loc=0, column="BioAct", value="Recessive")
    edge_results_rec_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    edge_results_rec_me_pb000.insert(loc=0, column="TestSeed", value=test_seed)
    Recessive_Results = pd.concat(
        [
            add_results_rec_me_pb000,
            DOMDEV_results_rec_me_pb000,
            rec_results_rec_me_pb000,
            dom_results_rec_me_pb000,
            codom_results_rec_me_pb000,
            edge_results_rec_me_pb000,
        ]
    )
    Recessive_Results_Final = pd.concat(
        [Recessive_Results_Final, Recessive_Results], axis=0
    )

    # Sub-Additive Main Effect for SNP1 without interaction
    # Training data
    train_sub_add_main_effect00 = sim.BAMS.from_model(
        eff1=sim.SNPEffectEncodings.SUB_ADDITIVE,
        eff2=sim.SNPEffectEncodings.ADDITIVE,
        penetrance_base=PEN_BASE,
        penetrance_diff=PEN_DIFF,
        main1=1,
        main2=0,
        interaction=0,
        snp1=variant1,
        snp2=variant2,
        random_seed=train_seed,
    )
    train_sub_add_me_pb000 = train_sub_add_main_effect00.generate_case_control(
        n_cases=n_cases, n_controls=n_controls, maf1=MAFA, maf2=MAFB, snr=SNR
    )

    # Calculate weights from the training dataset
    edge_weights_sub_add_me_pb000 = (
        train_sub_add_me_pb000.genomics.calculate_edge_encoding_values(
            data=train_sub_add_me_pb000["Outcome"], outcome_variable="Outcome"
        )
    )
    edge_weights_sub_add_me = edge_weights_sub_add_me_pb000.copy()
    edge_weights_sub_add_me.insert(
        loc=0, column="BioAct", value="Sub-Additive")
    edge_weights_sub_add_me.insert(loc=0, column="TrainSeed", value=train_seed)
    edge_weights_sub_add_me.insert(loc=0, column="TestSeed", value=test_seed)

    # Test data
    test_sub_add_main_effect00 = sim.BAMS.from_model(
        eff1=sim.SNPEffectEncodings.SUB_ADDITIVE,
        eff2=sim.SNPEffectEncodings.ADDITIVE,
        penetrance_base=PEN_BASE,
        penetrance_diff=PEN_DIFF,
        main1=1,
        main2=0,
        interaction=0,
        snp1=variant1,
        snp2=variant2,
        random_seed=test_seed,
    )
    test_sub_add_me_pb000 = test_sub_add_main_effect00.generate_case_control(
        n_cases=n_cases, n_controls=n_controls, maf1=MAFA, maf2=MAFB, snr=SNR
    )
    test_sub_add_me_pb000["Outcome"].cat.reorder_categories(
        ["Control", "Case"], inplace=True
    )

    # Run Regression by using weightes from CLARITE
    # Addtive Encoding
    test_sub_add_me_pb000_ADD = test_sub_add_me_pb000.genomics.encode_additive()
    add_results_sub_add_me_pb000 = clarite.analyze.association_study(
        data=test_sub_add_me_pb000_ADD, outcomes="Outcome"
    )
    add_results_sub_add_me_pb000["odds ratio"] = np.exp(
        add_results_sub_add_me_pb000["Beta"]
    )
    add_results_sub_add_me_pb000.insert(
        loc=0, column="Encoding", value="Additive")
    add_results_sub_add_me_pb000.insert(
        loc=0, column="BioAct", value="Sub-Additive")
    add_results_sub_add_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    add_results_sub_add_me_pb000.insert(
        loc=0, column="TestSeed", value=test_seed)

    # DOMDEV Encoding
    test_sub_add_me_pb000_DOMDEV = test_sub_add_me_pb000_ADD
    test_sub_add_me_pb000_DOMDEV["COV1"] = test_sub_add_me_pb000_DOMDEV["SNP1"]
    test_sub_add_me_pb000_DOMDEV["COV2"] = test_sub_add_me_pb000_DOMDEV["SNP2"]

    test_sub_add_me_pb000_DOMDEV["SNP1"] = test_sub_add_me_pb000_DOMDEV["SNP1"].replace(
        2, 0
    )
    test_sub_add_me_pb000_DOMDEV["SNP2"] = test_sub_add_me_pb000_DOMDEV["SNP2"].replace(
        2, 0
    )

    test_sub_add_me_pb000_DOMDEV_SNP1 = test_sub_add_me_pb000_DOMDEV[
        ["Outcome", "SNP1", "COV1"]
    ]
    test_sub_add_me_pb000_DOMDEV_SNP2 = test_sub_add_me_pb000_DOMDEV[
        ["Outcome", "SNP2", "COV2"]
    ]
    DOMDEV_results_sub_add_me_pb000_SNP1 = clarite.analyze.association_study(
        data=test_sub_add_me_pb000_DOMDEV_SNP1, outcomes="Outcome", covariates=["COV1"]
    )
    DOMDEV_results_sub_add_me_pb000_SNP2 = clarite.analyze.association_study(
        data=test_sub_add_me_pb000_DOMDEV_SNP2, outcomes="Outcome", covariates=["COV2"]
    )
    DOMDEV_results_sub_add_me_pb000 = pd.concat(
        [DOMDEV_results_sub_add_me_pb000_SNP1,
            DOMDEV_results_sub_add_me_pb000_SNP2]
    )

    DOMDEV_results_sub_add_me_pb000["odds ratio"] = np.exp(
        DOMDEV_results_sub_add_me_pb000["Beta"]
    )
    DOMDEV_results_sub_add_me_pb000.insert(
        loc=0, column="Encoding", value="DOMDEV")
    DOMDEV_results_sub_add_me_pb000.insert(
        loc=0, column="BioAct", value="Sub-Additive")
    DOMDEV_results_sub_add_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    DOMDEV_results_sub_add_me_pb000.insert(
        loc=0, column="TestSeed", value=test_seed)

    # Recessive Encoding
    test_sub_add_me_pb000_REC = test_sub_add_me_pb000.genomics.encode_recessive()
    rec_results_sub_add_me_pb000 = clarite.analyze.association_study(
        data=test_sub_add_me_pb000_REC, outcomes="Outcome"
    )
    rec_results_sub_add_me_pb000["odds ratio"] = np.exp(
        rec_results_sub_add_me_pb000["Beta"]
    )
    rec_results_sub_add_me_pb000.insert(
        loc=0, column="Encoding", value="Recessive")
    rec_results_sub_add_me_pb000.insert(
        loc=0, column="BioAct", value="Sub-Additive")
    rec_results_sub_add_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    rec_results_sub_add_me_pb000.insert(
        loc=0, column="TestSeed", value=test_seed)

    # Dominant Encoding
    test_sub_add_me_pb000_DOM = test_sub_add_me_pb000.genomics.encode_dominant()
    dom_results_sub_add_me_pb000 = clarite.analyze.association_study(
        data=test_sub_add_me_pb000_DOM, outcomes="Outcome"
    )
    dom_results_sub_add_me_pb000["odds ratio"] = np.exp(
        dom_results_sub_add_me_pb000["Beta"]
    )
    dom_results_sub_add_me_pb000.insert(
        loc=0, column="Encoding", value="Dominant")
    dom_results_sub_add_me_pb000.insert(
        loc=0, column="BioAct", value="Sub-Additive")
    dom_results_sub_add_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    dom_results_sub_add_me_pb000.insert(
        loc=0, column="TestSeed", value=test_seed)

    # Codominant Encoding
    test_sub_add_me_pb000_CODOM = test_sub_add_me_pb000.genomics.encode_codominant()
    codom_results_sub_add_me_pb000 = clarite.analyze.association_study(
        data=test_sub_add_me_pb000_CODOM, outcomes="Outcome"
    )
    codom_results_sub_add_me_pb000["odds ratio"] = np.exp(
        codom_results_sub_add_me_pb000["Beta"]
    )
    codom_results_sub_add_me_pb000.insert(
        loc=0, column="Encoding", value="Codominant")
    codom_results_sub_add_me_pb000.insert(
        loc=0, column="BioAct", value="Sub-Additive")
    codom_results_sub_add_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    codom_results_sub_add_me_pb000.insert(
        loc=0, column="TestSeed", value=test_seed)

    # EDGE Encoding
    test_sub_add_me_pb000_CLARITE = test_sub_add_me_pb000.genomics.encode_edge(
        encoding_info=edge_weights_sub_add_me_pb000
    )

    edge_results_sub_add_me_pb000 = clarite.analyze.association_study(
        data=test_sub_add_me_pb000_CLARITE, outcomes="Outcome"
    )
    edge_results_sub_add_me_pb000["odds ratio"] = np.exp(
        edge_results_sub_add_me_pb000["Beta"]
    )
    edge_results_sub_add_me_pb000.insert(
        loc=0, column="Encoding", value="EDGE")
    edge_results_sub_add_me_pb000.insert(
        loc=0, column="BioAct", value="Sub-Additive")
    edge_results_sub_add_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    edge_results_sub_add_me_pb000.insert(
        loc=0, column="TestSeed", value=test_seed)

    Sub_Additive_Results = pd.concat(
        [
            add_results_sub_add_me_pb000,
            DOMDEV_results_sub_add_me_pb000,
            rec_results_sub_add_me_pb000,
            dom_results_sub_add_me_pb000,
            codom_results_sub_add_me_pb000,
            edge_results_sub_add_me_pb000,
        ]
    )
    Sub_Additive_Results_Final = pd.concat(
        [Sub_Additive_Results_Final, Sub_Additive_Results], axis=0
    )

    # Additive Main Effect for SNP1 without interaction
    # Training data
    train_add_main_effect00 = sim.BAMS.from_model(
        eff1=sim.SNPEffectEncodings.ADDITIVE,
        eff2=sim.SNPEffectEncodings.ADDITIVE,
        penetrance_base=PEN_BASE,
        penetrance_diff=PEN_DIFF,
        main1=1,
        main2=0,
        interaction=0,
        snp1=variant1,
        snp2=variant2,
        random_seed=train_seed,
    )
    train_add_me_pb000 = train_add_main_effect00.generate_case_control(
        n_cases=n_cases, n_controls=n_controls, maf1=MAFA, maf2=MAFB, snr=SNR
    )

    # Calculate weights from the training dataset
    edge_weights_add_me_pb000 = (
        train_add_me_pb000.genomics.calculate_edge_encoding_values(
            data=train_add_me_pb000["Outcome"], outcome_variable="Outcome"
        )
    )
    edge_weights_add_me = edge_weights_add_me_pb000.copy()
    edge_weights_add_me.insert(loc=0, column="BioAct", value="Additive")
    edge_weights_add_me.insert(loc=0, column="TrainSeed", value=train_seed)
    edge_weights_add_me.insert(loc=0, column="TestSeed", value=test_seed)

    # Test data
    test_add_main_effect00 = sim.BAMS.from_model(
        eff1=sim.SNPEffectEncodings.ADDITIVE,
        eff2=sim.SNPEffectEncodings.ADDITIVE,
        penetrance_base=PEN_BASE,
        penetrance_diff=PEN_DIFF,
        main1=1,
        main2=0,
        interaction=0,
        snp1=variant1,
        snp2=variant2,
        random_seed=test_seed,
    )
    test_add_me_pb000 = test_add_main_effect00.generate_case_control(
        n_cases=n_cases, n_controls=n_controls, maf1=MAFA, maf2=MAFB, snr=SNR
    )
    test_add_me_pb000["Outcome"].cat.reorder_categories(
        ["Control", "Case"], inplace=True
    )

    # Run Regression by using weightes from CLARITE
    # Addtive Encoding
    test_add_me_pb000_ADD = test_add_me_pb000.genomics.encode_additive()
    add_results_add_me_pb000 = clarite.analyze.association_study(
        data=test_add_me_pb000_ADD, outcomes="Outcome"
    )
    add_results_add_me_pb000["odds ratio"] = np.exp(
        add_results_add_me_pb000["Beta"])
    add_results_add_me_pb000.insert(loc=0, column="Encoding", value="Additive")
    add_results_add_me_pb000.insert(loc=0, column="BioAct", value="Additive")
    add_results_add_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    add_results_add_me_pb000.insert(loc=0, column="TestSeed", value=test_seed)

    # DOMDEV Encoding
    test_add_me_pb000_DOMDEV = test_add_me_pb000_ADD
    test_add_me_pb000_DOMDEV["COV1"] = test_add_me_pb000_DOMDEV["SNP1"]
    test_add_me_pb000_DOMDEV["COV2"] = test_add_me_pb000_DOMDEV["SNP2"]

    test_add_me_pb000_DOMDEV["SNP1"] = test_add_me_pb000_DOMDEV["SNP1"].replace(
        2, 0)
    test_add_me_pb000_DOMDEV["SNP2"] = test_add_me_pb000_DOMDEV["SNP2"].replace(
        2, 0)

    test_add_me_pb000_DOMDEV_SNP1 = test_add_me_pb000_DOMDEV[
        ["Outcome", "SNP1", "COV1"]
    ]
    test_add_me_pb000_DOMDEV_SNP2 = test_add_me_pb000_DOMDEV[
        ["Outcome", "SNP2", "COV2"]
    ]
    DOMDEV_results_add_me_pb000_SNP1 = clarite.analyze.association_study(
        data=test_add_me_pb000_DOMDEV_SNP1, outcomes="Outcome", covariates=["COV1"]
    )
    DOMDEV_results_add_me_pb000_SNP2 = clarite.analyze.association_study(
        data=test_add_me_pb000_DOMDEV_SNP2, outcomes="Outcome", covariates=["COV2"]
    )
    DOMDEV_results_add_me_pb000 = pd.concat(
        [DOMDEV_results_add_me_pb000_SNP1, DOMDEV_results_add_me_pb000_SNP2]
    )

    DOMDEV_results_add_me_pb000["odds ratio"] = np.exp(
        DOMDEV_results_add_me_pb000["Beta"]
    )
    DOMDEV_results_add_me_pb000.insert(
        loc=0, column="Encoding", value="DOMDEV")
    DOMDEV_results_add_me_pb000.insert(
        loc=0, column="BioAct", value="Additive")
    DOMDEV_results_add_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    DOMDEV_results_add_me_pb000.insert(
        loc=0, column="TestSeed", value=test_seed)

    # Recessive Encoding
    test_add_me_pb000_REC = test_add_me_pb000.genomics.encode_recessive()
    rec_results_add_me_pb000 = clarite.analyze.association_study(
        data=test_add_me_pb000_REC, outcomes="Outcome"
    )
    rec_results_add_me_pb000["odds ratio"] = np.exp(
        rec_results_add_me_pb000["Beta"])
    rec_results_add_me_pb000.insert(
        loc=0, column="Encoding", value="Recessive")
    rec_results_add_me_pb000.insert(loc=0, column="BioAct", value="Additive")
    rec_results_add_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    rec_results_add_me_pb000.insert(loc=0, column="TestSeed", value=test_seed)

    # Dominant Encoding
    test_add_me_pb000_DOM = test_add_me_pb000.genomics.encode_dominant()
    dom_results_add_me_pb000 = clarite.analyze.association_study(
        data=test_add_me_pb000_DOM, outcomes="Outcome"
    )
    dom_results_add_me_pb000["odds ratio"] = np.exp(
        dom_results_add_me_pb000["Beta"])
    dom_results_add_me_pb000.insert(loc=0, column="Encoding", value="Dominant")
    dom_results_add_me_pb000.insert(loc=0, column="BioAct", value="Additive")
    dom_results_add_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    dom_results_add_me_pb000.insert(loc=0, column="TestSeed", value=test_seed)

    # Codominant Encoding
    test_add_me_pb000_CODOM = test_add_me_pb000.genomics.encode_codominant()
    codom_results_add_me_pb000 = clarite.analyze.association_study(
        data=test_add_me_pb000_CODOM, outcomes="Outcome"
    )
    codom_results_add_me_pb000["odds ratio"] = np.exp(
        codom_results_add_me_pb000["Beta"]
    )
    codom_results_add_me_pb000.insert(
        loc=0, column="Encoding", value="Codominant")
    codom_results_add_me_pb000.insert(loc=0, column="BioAct", value="Additive")
    codom_results_add_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    codom_results_add_me_pb000.insert(
        loc=0, column="TestSeed", value=test_seed)

    # EDGE Encoding
    test_add_me_pb000_CLARITE = test_add_me_pb000.genomics.encode_edge(
        encoding_info=edge_weights_add_me_pb000
    )

    edge_results_add_me_pb000 = clarite.analyze.association_study(
        data=test_add_me_pb000_CLARITE, outcomes="Outcome"
    )
    edge_results_add_me_pb000["odds ratio"] = np.exp(
        edge_results_add_me_pb000["Beta"])
    edge_results_add_me_pb000.insert(loc=0, column="Encoding", value="EDGE")
    edge_results_add_me_pb000.insert(loc=0, column="BioAct", value="Additive")
    edge_results_add_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    edge_results_add_me_pb000.insert(loc=0, column="TestSeed", value=test_seed)

    Additive_Results = pd.concat(
        [
            add_results_add_me_pb000,
            DOMDEV_results_add_me_pb000,
            rec_results_add_me_pb000,
            dom_results_add_me_pb000,
            codom_results_add_me_pb000,
            edge_results_add_me_pb000,
        ]
    )
    Additive_Results_Final = pd.concat(
        [Additive_Results_Final, Additive_Results], axis=0
    )

    # Super-Additive Main Effect for SNP1 without interaction
    # Training data
    train_sup_add_main_effect00 = sim.BAMS.from_model(
        eff1=sim.SNPEffectEncodings.SUPER_ADDITIVE,
        eff2=sim.SNPEffectEncodings.ADDITIVE,
        penetrance_base=PEN_BASE,
        penetrance_diff=PEN_DIFF,
        main1=1,
        main2=0,
        interaction=0,
        snp1=variant1,
        snp2=variant2,
        random_seed=train_seed,
    )
    train_sup_add_me_pb000 = train_sup_add_main_effect00.generate_case_control(
        n_cases=n_cases, n_controls=n_controls, maf1=MAFA, maf2=MAFB, snr=SNR
    )

    # Calculate weights from the training dataset
    edge_weights_sup_add_me_pb000 = (
        train_sup_add_me_pb000.genomics.calculate_edge_encoding_values(
            data=train_sup_add_me_pb000["Outcome"], outcome_variable="Outcome"
        )
    )
    edge_weights_sup_add_me = edge_weights_sup_add_me_pb000.copy()
    edge_weights_sup_add_me.insert(
        loc=0, column="BioAct", value="Super-Additive")
    edge_weights_sup_add_me.insert(loc=0, column="TrainSeed", value=train_seed)
    edge_weights_sup_add_me.insert(loc=0, column="TestSeed", value=test_seed)

    # Test data
    test_sup_add_main_effect00 = sim.BAMS.from_model(
        eff1=sim.SNPEffectEncodings.SUPER_ADDITIVE,
        eff2=sim.SNPEffectEncodings.ADDITIVE,
        penetrance_base=PEN_BASE,
        penetrance_diff=PEN_DIFF,
        main1=1,
        main2=0,
        interaction=0,
        snp1=variant1,
        snp2=variant2,
        random_seed=test_seed,
    )
    test_sup_add_me_pb000 = test_sup_add_main_effect00.generate_case_control(
        n_cases=n_cases, n_controls=n_controls, maf1=MAFA, maf2=MAFB, snr=SNR
    )
    test_sup_add_me_pb000["Outcome"].cat.reorder_categories(
        ["Control", "Case"], inplace=True
    )

    # Run Regression by using weightes from CLARITE
    # Addtive Encoding
    test_sup_add_me_pb000_ADD = test_sup_add_me_pb000.genomics.encode_additive()
    add_results_sup_add_me_pb000 = clarite.analyze.association_study(
        data=test_sup_add_me_pb000_ADD, outcomes="Outcome"
    )
    add_results_sup_add_me_pb000["odds ratio"] = np.exp(
        add_results_sup_add_me_pb000["Beta"]
    )
    add_results_sup_add_me_pb000.insert(
        loc=0, column="Encoding", value="Additive")
    add_results_sup_add_me_pb000.insert(
        loc=0, column="BioAct", value="Super-Additive")
    add_results_sup_add_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    add_results_sup_add_me_pb000.insert(
        loc=0, column="TestSeed", value=test_seed)

    # DOMDEV Encoding
    test_sup_add_me_pb000_DOMDEV = test_sup_add_me_pb000_ADD
    test_sup_add_me_pb000_DOMDEV["COV1"] = test_sup_add_me_pb000_DOMDEV["SNP1"]
    test_sup_add_me_pb000_DOMDEV["COV2"] = test_sup_add_me_pb000_DOMDEV["SNP2"]

    test_sup_add_me_pb000_DOMDEV["SNP1"] = test_sup_add_me_pb000_DOMDEV["SNP1"].replace(
        2, 0
    )
    test_sup_add_me_pb000_DOMDEV["SNP2"] = test_sup_add_me_pb000_DOMDEV["SNP2"].replace(
        2, 0
    )

    test_sup_add_me_pb000_DOMDEV_SNP1 = test_sup_add_me_pb000_DOMDEV[
        ["Outcome", "SNP1", "COV1"]
    ]
    test_sup_add_me_pb000_DOMDEV_SNP2 = test_sup_add_me_pb000_DOMDEV[
        ["Outcome", "SNP2", "COV2"]
    ]
    DOMDEV_results_sup_add_me_pb000_SNP1 = clarite.analyze.association_study(
        data=test_sup_add_me_pb000_DOMDEV_SNP1, outcomes="Outcome", covariates=["COV1"]
    )
    DOMDEV_results_sup_add_me_pb000_SNP2 = clarite.analyze.association_study(
        data=test_sup_add_me_pb000_DOMDEV_SNP2, outcomes="Outcome", covariates=["COV2"]
    )
    DOMDEV_results_sup_add_me_pb000 = pd.concat(
        [DOMDEV_results_sup_add_me_pb000_SNP1,
            DOMDEV_results_sup_add_me_pb000_SNP2]
    )

    DOMDEV_results_sup_add_me_pb000["odds ratio"] = np.exp(
        DOMDEV_results_sup_add_me_pb000["Beta"]
    )
    DOMDEV_results_sup_add_me_pb000.insert(
        loc=0, column="Encoding", value="DOMDEV")
    DOMDEV_results_sup_add_me_pb000.insert(
        loc=0, column="BioAct", value="Super-Additive"
    )
    DOMDEV_results_sup_add_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    DOMDEV_results_sup_add_me_pb000.insert(
        loc=0, column="TestSeed", value=test_seed)

    # Recessive Encoding
    test_sup_add_me_pb000_REC = test_sup_add_me_pb000.genomics.encode_recessive()
    rec_results_sup_add_me_pb000 = clarite.analyze.association_study(
        data=test_sup_add_me_pb000_REC, outcomes="Outcome"
    )
    rec_results_sup_add_me_pb000["odds ratio"] = np.exp(
        rec_results_sup_add_me_pb000["Beta"]
    )
    rec_results_sup_add_me_pb000.insert(
        loc=0, column="Encoding", value="Recessive")
    rec_results_sup_add_me_pb000.insert(
        loc=0, column="BioAct", value="Super-Additive")
    rec_results_sup_add_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    rec_results_sup_add_me_pb000.insert(
        loc=0, column="TestSeed", value=test_seed)

    # Dominant Encoding
    test_sup_add_me_pb000_DOM = test_sup_add_me_pb000.genomics.encode_dominant()
    dom_results_sup_add_me_pb000 = clarite.analyze.association_study(
        data=test_sup_add_me_pb000_DOM, outcomes="Outcome"
    )
    dom_results_sup_add_me_pb000["odds ratio"] = np.exp(
        dom_results_sup_add_me_pb000["Beta"]
    )
    dom_results_sup_add_me_pb000.insert(
        loc=0, column="Encoding", value="Dominant")
    dom_results_sup_add_me_pb000.insert(
        loc=0, column="BioAct", value="Super-Additive")
    dom_results_sup_add_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    dom_results_sup_add_me_pb000.insert(
        loc=0, column="TestSeed", value=test_seed)

    # Codominant Encoding
    test_sup_add_me_pb000_CODOM = test_sup_add_me_pb000.genomics.encode_codominant()
    codom_results_sup_add_me_pb000 = clarite.analyze.association_study(
        data=test_sup_add_me_pb000_CODOM, outcomes="Outcome"
    )
    codom_results_sup_add_me_pb000["odds ratio"] = np.exp(
        codom_results_sup_add_me_pb000["Beta"]
    )
    codom_results_sup_add_me_pb000.insert(
        loc=0, column="Encoding", value="Codominant")
    codom_results_sup_add_me_pb000.insert(
        loc=0, column="BioAct", value="Super-Additive"
    )
    codom_results_sup_add_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    codom_results_sup_add_me_pb000.insert(
        loc=0, column="TestSeed", value=test_seed)

    # EDGE Encoding
    test_sup_add_me_pb000_CLARITE = test_sup_add_me_pb000.genomics.encode_edge(
        encoding_info=edge_weights_sup_add_me_pb000
    )

    edge_results_sup_add_me_pb000 = clarite.analyze.association_study(
        data=test_sup_add_me_pb000_CLARITE, outcomes="Outcome"
    )
    edge_results_sup_add_me_pb000["odds ratio"] = np.exp(
        edge_results_sup_add_me_pb000["Beta"]
    )
    edge_results_sup_add_me_pb000.insert(
        loc=0, column="Encoding", value="EDGE")
    edge_results_sup_add_me_pb000.insert(
        loc=0, column="BioAct", value="Super-Additive")
    edge_results_sup_add_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    edge_results_sup_add_me_pb000.insert(
        loc=0, column="TestSeed", value=test_seed)

    Super_Additive_Results = pd.concat(
        [
            add_results_sup_add_me_pb000,
            DOMDEV_results_sup_add_me_pb000,
            rec_results_sup_add_me_pb000,
            dom_results_sup_add_me_pb000,
            codom_results_sup_add_me_pb000,
            edge_results_sup_add_me_pb000,
        ]
    )
    Super_Additive_Results_Final = pd.concat(
        [Super_Additive_Results_Final, Super_Additive_Results], axis=0
    )

    # Dominant Main Effect for SNP1 without interaction
    # Training data
    train_dom_main_effect00 = sim.BAMS.from_model(
        eff1=sim.SNPEffectEncodings.DOMINANT,
        eff2=sim.SNPEffectEncodings.ADDITIVE,
        penetrance_base=PEN_BASE,
        penetrance_diff=PEN_DIFF,
        main1=1,
        main2=0,
        interaction=0,
        snp1=variant1,
        snp2=variant2,
        random_seed=train_seed,
    )
    train_dom_me_pb000 = train_dom_main_effect00.generate_case_control(
        n_cases=n_cases, n_controls=n_controls, maf1=MAFA, maf2=MAFB, snr=SNR
    )

    # Calculate weights from the training dataset
    edge_weights_dom_me_pb000 = (
        train_dom_me_pb000.genomics.calculate_edge_encoding_values(
            data=train_dom_me_pb000["Outcome"], outcome_variable="Outcome"
        )
    )
    edge_weights_dom_me = edge_weights_dom_me_pb000.copy()
    edge_weights_dom_me.insert(loc=0, column="BioAct", value="Dominant")
    edge_weights_dom_me.insert(loc=0, column="TrainSeed", value=train_seed)
    edge_weights_dom_me.insert(loc=0, column="TestSeed", value=test_seed)

    # Test data
    test_dom_main_effect00 = sim.BAMS.from_model(
        eff1=sim.SNPEffectEncodings.DOMINANT,
        eff2=sim.SNPEffectEncodings.ADDITIVE,
        penetrance_base=PEN_BASE,
        penetrance_diff=PEN_DIFF,
        main1=1,
        main2=0,
        interaction=0,
        snp1=variant1,
        snp2=variant2,
        random_seed=test_seed,
    )
    test_dom_me_pb000 = test_dom_main_effect00.generate_case_control(
        n_cases=n_cases, n_controls=n_controls, maf1=MAFA, maf2=MAFB, snr=SNR
    )
    test_dom_me_pb000["Outcome"].cat.reorder_categories(
        ["Control", "Case"], inplace=True
    )

    # Run Regression by using weightes from CLARITE
    # Addtive Encoding
    test_dom_me_pb000_ADD = test_dom_me_pb000.genomics.encode_additive()
    add_results_dom_me_pb000 = clarite.analyze.association_study(
        data=test_dom_me_pb000_ADD, outcomes="Outcome"
    )
    add_results_dom_me_pb000["odds ratio"] = np.exp(
        add_results_dom_me_pb000["Beta"])
    add_results_dom_me_pb000.insert(loc=0, column="Encoding", value="Additive")
    add_results_dom_me_pb000.insert(loc=0, column="BioAct", value="Dominant")
    add_results_dom_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    add_results_dom_me_pb000.insert(loc=0, column="TestSeed", value=test_seed)

    # DOMDEV Encoding
    test_dom_me_pb000_DOMDEV = test_dom_me_pb000_ADD
    test_dom_me_pb000_DOMDEV["COV1"] = test_dom_me_pb000_DOMDEV["SNP1"]
    test_dom_me_pb000_DOMDEV["COV2"] = test_dom_me_pb000_DOMDEV["SNP2"]

    test_dom_me_pb000_DOMDEV["SNP1"] = test_dom_me_pb000_DOMDEV["SNP1"].replace(
        2, 0)
    test_dom_me_pb000_DOMDEV["SNP2"] = test_dom_me_pb000_DOMDEV["SNP2"].replace(
        2, 0)

    test_dom_me_pb000_DOMDEV_SNP1 = test_dom_me_pb000_DOMDEV[
        ["Outcome", "SNP1", "COV1"]
    ]
    test_dom_me_pb000_DOMDEV_SNP2 = test_dom_me_pb000_DOMDEV[
        ["Outcome", "SNP2", "COV2"]
    ]
    DOMDEV_results_dom_me_pb000_SNP1 = clarite.analyze.association_study(
        data=test_dom_me_pb000_DOMDEV_SNP1, outcomes="Outcome", covariates=["COV1"]
    )
    DOMDEV_results_dom_me_pb000_SNP2 = clarite.analyze.association_study(
        data=test_dom_me_pb000_DOMDEV_SNP2, outcomes="Outcome", covariates=["COV2"]
    )
    DOMDEV_results_dom_me_pb000 = pd.concat(
        [DOMDEV_results_dom_me_pb000_SNP1, DOMDEV_results_dom_me_pb000_SNP2]
    )

    DOMDEV_results_dom_me_pb000["odds ratio"] = np.exp(
        DOMDEV_results_dom_me_pb000["Beta"]
    )
    DOMDEV_results_dom_me_pb000.insert(
        loc=0, column="Encoding", value="DOMDEV")
    DOMDEV_results_dom_me_pb000.insert(
        loc=0, column="BioAct", value="Dominant")
    DOMDEV_results_dom_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    DOMDEV_results_dom_me_pb000.insert(
        loc=0, column="TestSeed", value=test_seed)

    # Recessive Encoding
    test_dom_me_pb000_REC = test_dom_me_pb000.genomics.encode_recessive()
    rec_results_dom_me_pb000 = clarite.analyze.association_study(
        data=test_dom_me_pb000_REC, outcomes="Outcome"
    )
    rec_results_dom_me_pb000["odds ratio"] = np.exp(
        rec_results_dom_me_pb000["Beta"])
    rec_results_dom_me_pb000.insert(
        loc=0, column="Encoding", value="Recessive")
    rec_results_dom_me_pb000.insert(loc=0, column="BioAct", value="Dominant")
    rec_results_dom_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    rec_results_dom_me_pb000.insert(loc=0, column="TestSeed", value=test_seed)

    # Dominant Encoding
    test_dom_me_pb000_DOM = test_dom_me_pb000.genomics.encode_dominant()
    dom_results_dom_me_pb000 = clarite.analyze.association_study(
        data=test_dom_me_pb000_DOM, outcomes="Outcome"
    )
    dom_results_dom_me_pb000["odds ratio"] = np.exp(
        dom_results_dom_me_pb000["Beta"])
    dom_results_dom_me_pb000.insert(loc=0, column="Encoding", value="Dominant")
    dom_results_dom_me_pb000.insert(loc=0, column="BioAct", value="Dominant")
    dom_results_dom_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    dom_results_dom_me_pb000.insert(loc=0, column="TestSeed", value=test_seed)

    # Codominant Encoding
    test_dom_me_pb000_CODOM = test_dom_me_pb000.genomics.encode_codominant()
    codom_results_dom_me_pb000 = clarite.analyze.association_study(
        data=test_dom_me_pb000_CODOM, outcomes="Outcome"
    )
    codom_results_dom_me_pb000["odds ratio"] = np.exp(
        codom_results_dom_me_pb000["Beta"]
    )
    codom_results_dom_me_pb000.insert(
        loc=0, column="Encoding", value="Codominant")
    codom_results_dom_me_pb000.insert(loc=0, column="BioAct", value="Dominant")
    codom_results_dom_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    codom_results_dom_me_pb000.insert(
        loc=0, column="TestSeed", value=test_seed)

    # EDGE Encoding
    test_dom_me_pb000_CLARITE = test_dom_me_pb000.genomics.encode_edge(
        encoding_info=edge_weights_dom_me_pb000
    )

    edge_results_dom_me_pb000 = clarite.analyze.association_study(
        data=test_dom_me_pb000_CLARITE, outcomes="Outcome"
    )
    edge_results_dom_me_pb000["odds ratio"] = np.exp(
        edge_results_dom_me_pb000["Beta"])
    edge_results_dom_me_pb000.insert(loc=0, column="Encoding", value="EDGE")
    edge_results_dom_me_pb000.insert(loc=0, column="BioAct", value="Dominant")
    edge_results_dom_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    edge_results_dom_me_pb000.insert(loc=0, column="TestSeed", value=test_seed)

    Dominant_Results = pd.concat(
        [
            add_results_dom_me_pb000,
            DOMDEV_results_dom_me_pb000,
            rec_results_dom_me_pb000,
            dom_results_dom_me_pb000,
            codom_results_dom_me_pb000,
            edge_results_dom_me_pb000,
        ]
    )
    Dominant_Results_Final = pd.concat(
        [Dominant_Results_Final, Dominant_Results], axis=0
    )

    # Heterozygous Main Effect for SNP1 without interaction
    # Training data
    train_het_main_effect00 = sim.BAMS.from_model(
        eff1=sim.SNPEffectEncodings.HET,
        eff2=sim.SNPEffectEncodings.ADDITIVE,
        penetrance_base=PEN_BASE,
        penetrance_diff=PEN_DIFF,
        main1=1,
        main2=0,
        interaction=0,
        snp1=variant1,
        snp2=variant2,
        random_seed=train_seed,
    )
    train_het_me_pb000 = train_het_main_effect00.generate_case_control(
        n_cases=n_cases, n_controls=n_controls, maf1=MAFA, maf2=MAFB, snr=SNR
    )

    # Calculate weights from the training dataset
    edge_weights_het_me_pb000 = (
        train_het_me_pb000.genomics.calculate_edge_encoding_values(
            data=train_het_me_pb000["Outcome"], outcome_variable="Outcome"
        )
    )
    edge_weights_het_me = edge_weights_het_me_pb000.copy()
    edge_weights_het_me.insert(loc=0, column="BioAct", value="Heterozygous")
    edge_weights_het_me.insert(loc=0, column="TrainSeed", value=train_seed)
    edge_weights_het_me.insert(loc=0, column="TestSeed", value=test_seed)

    # Test data
    test_het_main_effect00 = sim.BAMS.from_model(
        eff1=sim.SNPEffectEncodings.HET,
        eff2=sim.SNPEffectEncodings.ADDITIVE,
        penetrance_base=PEN_BASE,
        penetrance_diff=PEN_DIFF,
        main1=1,
        main2=0,
        interaction=0,
        snp1=variant1,
        snp2=variant2,
        random_seed=test_seed,
    )
    test_het_me_pb000 = test_het_main_effect00.generate_case_control(
        n_cases=n_cases, n_controls=n_controls, maf1=MAFA, maf2=MAFB, snr=SNR
    )
    test_het_me_pb000["Outcome"].cat.reorder_categories(
        ["Control", "Case"], inplace=True
    )

    # Run Regression by using weightes from CLARITE
    # Addtive Encoding
    test_het_me_pb000_ADD = test_het_me_pb000.genomics.encode_additive()
    add_results_het_me_pb000 = clarite.analyze.association_study(
        data=test_het_me_pb000_ADD, outcomes="Outcome"
    )
    add_results_het_me_pb000["odds ratio"] = np.exp(
        add_results_het_me_pb000["Beta"])
    add_results_het_me_pb000.insert(loc=0, column="Encoding", value="Additive")
    add_results_het_me_pb000.insert(
        loc=0, column="BioAct", value="Heterozygous")
    add_results_het_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    add_results_het_me_pb000.insert(loc=0, column="TestSeed", value=test_seed)

    # DOMDEV Encoding
    test_het_me_pb000_DOMDEV = test_het_me_pb000_ADD
    test_het_me_pb000_DOMDEV["COV1"] = test_het_me_pb000_DOMDEV["SNP1"]
    test_het_me_pb000_DOMDEV["COV2"] = test_het_me_pb000_DOMDEV["SNP2"]

    test_het_me_pb000_DOMDEV["SNP1"] = test_het_me_pb000_DOMDEV["SNP1"].replace(
        2, 0)
    test_het_me_pb000_DOMDEV["SNP2"] = test_het_me_pb000_DOMDEV["SNP2"].replace(
        2, 0)

    test_het_me_pb000_DOMDEV_SNP1 = test_het_me_pb000_DOMDEV[
        ["Outcome", "SNP1", "COV1"]
    ]
    test_het_me_pb000_DOMDEV_SNP2 = test_het_me_pb000_DOMDEV[
        ["Outcome", "SNP2", "COV2"]
    ]
    DOMDEV_results_het_me_pb000_SNP1 = clarite.analyze.association_study(
        data=test_het_me_pb000_DOMDEV_SNP1, outcomes="Outcome", covariates=["COV1"]
    )
    DOMDEV_results_het_me_pb000_SNP2 = clarite.analyze.association_study(
        data=test_het_me_pb000_DOMDEV_SNP2, outcomes="Outcome", covariates=["COV2"]
    )
    DOMDEV_results_het_me_pb000 = pd.concat(
        [DOMDEV_results_het_me_pb000_SNP1, DOMDEV_results_het_me_pb000_SNP2]
    )

    DOMDEV_results_het_me_pb000["odds ratio"] = np.exp(
        DOMDEV_results_het_me_pb000["Beta"]
    )
    DOMDEV_results_het_me_pb000.insert(
        loc=0, column="Encoding", value="DOMDEV")
    DOMDEV_results_het_me_pb000.insert(
        loc=0, column="BioAct", value="Heterozygous")
    DOMDEV_results_het_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    DOMDEV_results_het_me_pb000.insert(
        loc=0, column="TestSeed", value=test_seed)

    # Recessive Encoding
    test_het_me_pb000_REC = test_het_me_pb000.genomics.encode_recessive()
    rec_results_het_me_pb000 = clarite.analyze.association_study(
        data=test_het_me_pb000_REC, outcomes="Outcome"
    )
    rec_results_het_me_pb000["odds ratio"] = np.exp(
        rec_results_het_me_pb000["Beta"])
    rec_results_het_me_pb000.insert(
        loc=0, column="Encoding", value="Recessive")
    rec_results_het_me_pb000.insert(
        loc=0, column="BioAct", value="Heterozygous")
    rec_results_het_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    rec_results_het_me_pb000.insert(loc=0, column="TestSeed", value=test_seed)

    # Dominant Encoding
    test_het_me_pb000_DOM = test_het_me_pb000.genomics.encode_dominant()
    dom_results_het_me_pb000 = clarite.analyze.association_study(
        data=test_het_me_pb000_DOM, outcomes="Outcome"
    )
    dom_results_het_me_pb000["odds ratio"] = np.exp(
        dom_results_het_me_pb000["Beta"])
    dom_results_het_me_pb000.insert(loc=0, column="Encoding", value="Dominant")
    dom_results_het_me_pb000.insert(
        loc=0, column="BioAct", value="Heterozygous")
    dom_results_het_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    dom_results_het_me_pb000.insert(loc=0, column="TestSeed", value=test_seed)

    # Codominant Encoding
    test_het_me_pb000_CODOM = test_het_me_pb000.genomics.encode_codominant()
    codom_results_het_me_pb000 = clarite.analyze.association_study(
        data=test_het_me_pb000_CODOM, outcomes="Outcome"
    )
    codom_results_het_me_pb000["odds ratio"] = np.exp(
        codom_results_het_me_pb000["Beta"]
    )
    codom_results_het_me_pb000.insert(
        loc=0, column="Encoding", value="Codominant")
    codom_results_het_me_pb000.insert(
        loc=0, column="BioAct", value="Heterozygous")
    codom_results_het_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    codom_results_het_me_pb000.insert(
        loc=0, column="TestSeed", value=test_seed)

    # EDGE Encoding
    test_het_me_pb000_CLARITE = test_het_me_pb000.genomics.encode_edge(
        encoding_info=edge_weights_het_me_pb000
    )

    edge_results_het_me_pb000 = clarite.analyze.association_study(
        data=test_het_me_pb000_CLARITE, outcomes="Outcome"
    )
    edge_results_het_me_pb000["odds ratio"] = np.exp(
        edge_results_het_me_pb000["Beta"])
    edge_results_het_me_pb000.insert(loc=0, column="Encoding", value="EDGE")
    edge_results_het_me_pb000.insert(
        loc=0, column="BioAct", value="Heterozygous")
    edge_results_het_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    edge_results_het_me_pb000.insert(loc=0, column="TestSeed", value=test_seed)

    Het_Results = pd.concat(
        [
            add_results_het_me_pb000,
            DOMDEV_results_het_me_pb000,
            rec_results_het_me_pb000,
            dom_results_het_me_pb000,
            codom_results_het_me_pb000,
            edge_results_het_me_pb000,
        ]
    )
    Het_Results_Final = pd.concat([Het_Results_Final, Het_Results], axis=0)

    # NULL Main Effect for SNP1 without interaction
    # Training data
    train_null_main_effect00 = sim.BAMS.from_model(
        eff1=sim.SNPEffectEncodings.ADDITIVE,
        eff2=sim.SNPEffectEncodings.ADDITIVE,
        penetrance_base=PEN_BASE,
        penetrance_diff=PEN_DIFF,
        main1=0,
        main2=0,
        interaction=0,
        snp1=variant1,
        snp2=variant2,
        random_seed=train_seed,
    )
    train_null_me_pb000 = train_null_main_effect00.generate_case_control(
        n_cases=n_cases, n_controls=n_controls, maf1=MAFA, maf2=MAFB
    )

    # Calculate weights from the training dataset
    edge_weights_null_me_pb000 = (
        train_null_me_pb000.genomics.calculate_edge_encoding_values(
            data=train_null_me_pb000["Outcome"], outcome_variable="Outcome"
        )
    )
    edge_weights_null_me = edge_weights_null_me_pb000.copy()
    edge_weights_null_me.insert(loc=0, column="BioAct", value="NULL")
    edge_weights_null_me.insert(loc=0, column="TrainSeed", value=train_seed)
    edge_weights_null_me.insert(loc=0, column="TestSeed", value=test_seed)

    # Test data
    test_null_main_effect00 = sim.BAMS.from_model(
        eff1=sim.SNPEffectEncodings.ADDITIVE,
        eff2=sim.SNPEffectEncodings.ADDITIVE,
        penetrance_base=PEN_BASE,
        penetrance_diff=PEN_DIFF,
        main1=0,
        main2=0,
        interaction=0,
        snp1=variant1,
        snp2=variant2,
        random_seed=test_seed,
    )
    test_null_me_pb000 = test_null_main_effect00.generate_case_control(
        n_cases=n_cases, n_controls=n_controls, maf1=MAFA, maf2=MAFB
    )
    test_null_me_pb000["Outcome"].cat.reorder_categories(
        ["Control", "Case"], inplace=True
    )

    # Run Regression by using weightes from CLARITE
    # Addtive Encoding
    test_null_me_pb000_ADD = test_null_me_pb000.genomics.encode_additive()
    add_results_null_me_pb000 = clarite.analyze.association_study(
        data=test_null_me_pb000_ADD, outcomes="Outcome"
    )
    add_results_null_me_pb000["odds ratio"] = np.exp(
        add_results_null_me_pb000["Beta"])
    add_results_null_me_pb000.insert(
        loc=0, column="Encoding", value="Additive")
    add_results_null_me_pb000.insert(loc=0, column="BioAct", value="NULL")
    add_results_null_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    add_results_null_me_pb000.insert(loc=0, column="TestSeed", value=test_seed)

    # DOMDEV Encoding
    test_null_me_pb000_DOMDEV = test_null_me_pb000_ADD
    test_null_me_pb000_DOMDEV["COV1"] = test_null_me_pb000_DOMDEV["SNP1"]
    test_null_me_pb000_DOMDEV["COV2"] = test_null_me_pb000_DOMDEV["SNP2"]

    test_null_me_pb000_DOMDEV["SNP1"] = test_null_me_pb000_DOMDEV["SNP1"].replace(
        2, 0)
    test_null_me_pb000_DOMDEV["SNP2"] = test_null_me_pb000_DOMDEV["SNP2"].replace(
        2, 0)

    test_null_me_pb000_DOMDEV_SNP1 = test_null_me_pb000_DOMDEV[
        ["Outcome", "SNP1", "COV1"]
    ]
    test_null_me_pb000_DOMDEV_SNP2 = test_null_me_pb000_DOMDEV[
        ["Outcome", "SNP2", "COV2"]
    ]
    DOMDEV_results_null_me_pb000_SNP1 = clarite.analyze.association_study(
        data=test_null_me_pb000_DOMDEV_SNP1, outcomes="Outcome", covariates=["COV1"]
    )
    DOMDEV_results_null_me_pb000_SNP2 = clarite.analyze.association_study(
        data=test_null_me_pb000_DOMDEV_SNP2, outcomes="Outcome", covariates=["COV2"]
    )
    DOMDEV_results_null_me_pb000 = pd.concat(
        [DOMDEV_results_null_me_pb000_SNP1, DOMDEV_results_null_me_pb000_SNP2]
    )

    DOMDEV_results_null_me_pb000["odds ratio"] = np.exp(
        DOMDEV_results_null_me_pb000["Beta"]
    )
    DOMDEV_results_null_me_pb000.insert(
        loc=0, column="Encoding", value="DOMDEV")
    DOMDEV_results_null_me_pb000.insert(loc=0, column="BioAct", value="NULL")
    DOMDEV_results_null_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    DOMDEV_results_null_me_pb000.insert(
        loc=0, column="TestSeed", value=test_seed)

    # Recessive Encoding
    test_null_me_pb000_REC = test_null_me_pb000.genomics.encode_recessive()
    rec_results_null_me_pb000 = clarite.analyze.association_study(
        data=test_null_me_pb000_REC, outcomes="Outcome"
    )
    rec_results_null_me_pb000["odds ratio"] = np.exp(
        rec_results_null_me_pb000["Beta"])
    rec_results_null_me_pb000.insert(
        loc=0, column="Encoding", value="Recessive")
    rec_results_null_me_pb000.insert(loc=0, column="BioAct", value="NULL")
    rec_results_null_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    rec_results_null_me_pb000.insert(loc=0, column="TestSeed", value=test_seed)

    # Dominant Encoding
    test_null_me_pb000_DOM = test_null_me_pb000.genomics.encode_dominant()
    dom_results_null_me_pb000 = clarite.analyze.association_study(
        data=test_null_me_pb000_DOM, outcomes="Outcome"
    )
    dom_results_null_me_pb000["odds ratio"] = np.exp(
        dom_results_null_me_pb000["Beta"])
    dom_results_null_me_pb000.insert(
        loc=0, column="Encoding", value="Dominant")
    dom_results_null_me_pb000.insert(loc=0, column="BioAct", value="NULL")
    dom_results_null_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    dom_results_null_me_pb000.insert(loc=0, column="TestSeed", value=test_seed)

    # Codominant Encoding
    test_null_me_pb000_CODOM = test_null_me_pb000.genomics.encode_codominant()
    codom_results_null_me_pb000 = clarite.analyze.association_study(
        data=test_null_me_pb000_CODOM, outcomes="Outcome"
    )
    codom_results_null_me_pb000["odds ratio"] = np.exp(
        codom_results_null_me_pb000["Beta"]
    )
    codom_results_null_me_pb000.insert(
        loc=0, column="Encoding", value="Codominant")
    codom_results_null_me_pb000.insert(loc=0, column="BioAct", value="NULL")
    codom_results_null_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    codom_results_null_me_pb000.insert(
        loc=0, column="TestSeed", value=test_seed)

    # EDGE Encoding
    test_null_me_pb000_CLARITE = test_null_me_pb000.genomics.encode_edge(
        encoding_info=edge_weights_null_me_pb000
    )

    edge_results_null_me_pb000 = clarite.analyze.association_study(
        data=test_null_me_pb000_CLARITE, outcomes="Outcome"
    )
    edge_results_null_me_pb000["odds ratio"] = np.exp(
        edge_results_null_me_pb000["Beta"]
    )
    edge_results_null_me_pb000.insert(loc=0, column="Encoding", value="EDGE")
    edge_results_null_me_pb000.insert(loc=0, column="BioAct", value="NULL")
    edge_results_null_me_pb000.insert(
        loc=0, column="TrainSeed", value=train_seed)
    edge_results_null_me_pb000.insert(
        loc=0, column="TestSeed", value=test_seed)

    NULL_Results = pd.concat(
        [
            add_results_null_me_pb000,
            DOMDEV_results_null_me_pb000,
            rec_results_null_me_pb000,
            dom_results_null_me_pb000,
            codom_results_null_me_pb000,
            edge_results_null_me_pb000,
        ]
    )
    NULL_Results_Final = pd.concat([NULL_Results_Final, NULL_Results], axis=0)
    # NULL_Results_Final.to_csv(
    #    f"/storage/home/jpz5091/work/bams/Sim1/NULL_Results_{num_samples}_{case_control_ratio}_pb{PEN_BASE}_pd{PEN_DIFF}_maf{MAFA}_snr{SNR}.txt",
    #    sep=";",
    # )

    EDGE_Results = pd.concat(
        [
            edge_results_rec_me_pb000,
            edge_results_sub_add_me_pb000,
            edge_results_add_me_pb000,
            edge_results_sup_add_me_pb000,
            edge_results_dom_me_pb000,
            edge_results_het_me_pb000,
        ]
    )
    EDGE_Results_Final = pd.concat([EDGE_Results_Final, EDGE_Results], axis=0)
    # EDGE_Results_Final.to_csv(
    #    f"/storage/home/jpz5091/work/bams/Sim1/EDGE_Results_{num_samples}_{case_control_ratio}_pb{PEN_BASE}_pd{PEN_DIFF}_maf{MAFA}_snr{SNR}.txt",
    #    sep=";",
    # )

    EDGE_alpha_Results = pd.concat(
        [
            edge_weights_rec_me,
            edge_weights_sub_add_me,
            edge_weights_add_me,
            edge_weights_sup_add_me,
            edge_weights_dom_me,
            edge_weights_het_me,
            edge_weights_null_me,
        ]
    )
    EDGE_alpha_Final = pd.concat(
        [EDGE_alpha_Final, EDGE_alpha_Results], axis=0)
    # EDGE_alpha_Final.to_csv(
    #    f"/storage/home/jpz5091/work/bams/Sim1/EDGE_alpha_Results_{num_samples}_{case_control_ratio}_pb{PEN_BASE}_pd{PEN_DIFF}_maf{MAFA}_snr{SNR}.txt",
    #    sep=";",
    # )

    All_Results = pd.concat(
        [
            Recessive_Results_Final,
            Sub_Additive_Results_Final,
            Additive_Results_Final,
            Super_Additive_Results_Final,
            Dominant_Results_Final,
            Het_Results_Final,
            NULL_Results_Final,
        ]
    )
    # All_Results.to_csv(
    #   f"/storage/home/jpz5091/work/bams/Sim1/All_Results_{num_samples}_{case_control_ratio}_pb{PEN_BASE}_pd{PEN_DIFF}_maf{MAFA}_snr{SNR}.txt",
    #    sep=";",
    # )

    endcycle = time.time()
    print("The time of execution of one cycle is :", endcycle - startcycle)


# instantiating the decorator
# @profile


def mp_treads():
    with concurrent.futures.ProcessPoolExecutor() as executor:

        results = executor.map(simulations, n_loops)

        #All_Results_Finals = pd.DataFrame()

        # for vpid in results:
        # print(vpid)
        # All_Results_Finals = pd.concat(
        #    [All_Results_Finals, All_Results], axis=0)

        # alpha_Finals.to_csv("~/Public/teste1.txt", sep=";")
        # All_Results_Finals.to_csv("~/Public/teste2.txt", sep=";")

    print(f"Finished in {round(time.perf_counter()-start,2)} second(s)")


if __name__ == "__main__":
    mp_treads()
