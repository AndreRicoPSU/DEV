
    ## Recessive Main Effect for SNP1 without interaction
    # Training data
    train_rec_main_effect00 = sim.BAMS.from_model(
        eff1=sim.SNPEffectEncodings.RECESSIVE,
        eff2=sim.SNPEffectEncodings.ADDITIVE,
    edge_weights_rec_me.insert(loc=0, column="BioAct", value="Recessive")
    test_rec_main_effect00 = sim.BAMS.from_model(
        eff1=sim.SNPEffectEncodings.RECESSIVE,
        eff2=sim.SNPEffectEncodings.ADDITIVE,    
    add_results_rec_me_pb000.insert(loc=0, column="Encoding", value="Additive")
    add_results_rec_me_pb000.insert(loc=0, column="BioAct", value="Recessive")    
    DOMDEV_results_rec_me_pb000.insert(loc=0, column="Encoding", value="DOMDEV")
    DOMDEV_results_rec_me_pb000.insert(loc=0, column="BioAct", value="Recessive")
    rec_results_rec_me_pb000.insert(loc=0, column="BioAct", value="Recessive")    
    dom_results_rec_me_pb000.insert(loc=0, column="Encoding", value="Dominant")
    dom_results_rec_me_pb000.insert(loc=0, column="BioAct", value="Recessive")
    codom_results_rec_me_pb000.insert(loc=0, column="Encoding", value="Codominant")
    codom_results_rec_me_pb000.insert(loc=0, column="BioAct", value="Recessive")
    edge_results_rec_me_pb000.insert(loc=0, column="Encoding", value="EDGE")
    edge_results_rec_me_pb000.insert(loc=0, column="BioAct", value="Recessive")


    ## Sub-Additive Main Effect for SNP1 without interaction
    # Training data
    train_sub_add_main_effect00 = sim.BAMS.from_model(
        eff1=sim.SNPEffectEncodings.SUB_ADDITIVE,
        eff2=sim.SNPEffectEncodings.ADDITIVE,
    edge_weights_sub_add_me.insert(loc=0, column="BioAct", value="Sub-Additive")
      test_sub_add_main_effect00 = sim.BAMS.from_model(
        eff1=sim.SNPEffectEncodings.SUB_ADDITIVE,
        eff2=sim.SNPEffectEncodings.ADDITIVE,
    add_results_sub_add_me_pb000.insert(loc=0, column="Encoding", value="Additive")
    add_results_sub_add_me_pb000.insert(loc=0, column="BioAct", value="Sub-Additive")    )
    DOMDEV_results_sub_add_me_pb000.insert(loc=0, column="Encoding", value="DOMDEV")
    DOMDEV_results_sub_add_me_pb000.insert(loc=0, column="BioAct", value="Sub-Additive")
    rec_results_sub_add_me_pb000.insert(loc=0, column="Encoding", value="Recessive")
    rec_results_sub_add_me_pb000.insert(loc=0, column="BioAct", value="Sub-Additive")
    dom_results_sub_add_me_pb000.insert(loc=0, column="Encoding", value="Dominant")
    dom_results_sub_add_me_pb000.insert(loc=0, column="BioAct", value="Sub-Additive")
    codom_results_sub_add_me_pb000.insert(loc=0, column="Encoding", value="Codominant")
    codom_results_sub_add_me_pb000.insert(loc=0, column="BioAct", value="Sub-Additive")
    edge_results_sub_add_me_pb000.insert(loc=0, column="Encoding", value="EDGE")
    edge_results_sub_add_me_pb000.insert(loc=0, column="BioAct", value="Sub-Additive")


    ## Additive Main Effect for SNP1 without interaction
    # Training data
    train_add_main_effect00 = sim.BAMS.from_model(
        eff1=sim.SNPEffectEncodings.ADDITIVE,
        eff2=sim.SNPEffectEncodings.ADDITIVE,
     edge_weights_add_me.insert(loc=0, column="BioAct", value="Additive")
    test_add_main_effect00 = sim.BAMS.from_model(
        eff1=sim.SNPEffectEncodings.ADDITIVE,
        eff2=sim.SNPEffectEncodings.ADDITIVE,
    add_results_add_me_pb000.insert(loc=0, column="Encoding", value="Additive")
    add_results_add_me_pb000.insert(loc=0, column="BioAct", value="Additive")
    DOMDEV_results_add_me_pb000.insert(loc=0, column="Encoding", value="DOMDEV")
    DOMDEV_results_add_me_pb000.insert(loc=0, column="BioAct", value="Additive")
    rec_results_add_me_pb000.insert(loc=0, column="Encoding", value="Recessive")
    rec_results_add_me_pb000.insert(loc=0, column="BioAct", value="Additive")
    dom_results_add_me_pb000.insert(loc=0, column="Encoding", value="Dominant")
    dom_results_add_me_pb000.insert(loc=0, column="BioAct", value="Additive")
    codom_results_add_me_pb000.insert(loc=0, column="Encoding", value="Codominant")
    codom_results_add_me_pb000.insert(loc=0, column="BioAct", value="Additive")
    edge_results_add_me_pb000.insert(loc=0, column="Encoding", value="EDGE")
    edge_results_add_me_pb000.insert(loc=0, column="BioAct", value="Additive")

    ## Super-Additive Main Effect for SNP1 without interaction
    # Training data
    train_sup_add_main_effect00 = sim.BAMS.from_model(
        eff1=sim.SNPEffectEncodings.SUPER_ADDITIVE,
        eff2=sim.SNPEffectEncodings.ADDITIVE,
    edge_weights_sup_add_me.insert(loc=0, column="BioAct", value="Super-Additive")
    test_sup_add_main_effect00 = sim.BAMS.from_model(
        eff1=sim.SNPEffectEncodings.SUPER_ADDITIVE,
        eff2=sim.SNPEffectEncodings.ADDITIVE,
    add_results_sup_add_me_pb000.insert(loc=0, column="Encoding", value="Additive")
    add_results_sup_add_me_pb000.insert(loc=0, column="BioAct", value="Super-Additive")
    DOMDEV_results_sup_add_me_pb000.insert(loc=0, column="Encoding", value="DOMDEV")
    DOMDEV_results_sup_add_me_pb000.insert(loc=0, column="BioAct", value="Super-Additive")
    rec_results_sup_add_me_pb000.insert(loc=0, column="Encoding", value="Recessive")
    rec_results_sup_add_me_pb000.insert(loc=0, column="BioAct", value="Super-Additive")
    dom_results_sup_add_me_pb000.insert(loc=0, column="Encoding", value="Dominant")
    dom_results_sup_add_me_pb000.insert(loc=0, column="BioAct", value="Super-Additive")
    codom_results_sup_add_me_pb000.insert(loc=0, column="Encoding", value="Codominant")
    codom_results_sup_add_me_pb000.insert(loc=0, column="BioAct", value="Super-Additive")
    edge_results_sup_add_me_pb000.insert(loc=0, column="Encoding", value="EDGE")
    edge_results_sup_add_me_pb000.insert(loc=0, column="BioAct", value="Super-Additive")

    ## Dominant Main Effect for SNP1 without interaction
    # Training data
    train_dom_main_effect00 = sim.BAMS.from_model(
        eff1=sim.SNPEffectEncodings.DOMINANT,
        eff2=sim.SNPEffectEncodings.ADDITIVE,
    edge_weights_dom_me.insert(loc=0, column="BioAct", value="Dominant")
    test_dom_main_effect00 = sim.BAMS.from_model(
        eff1=sim.SNPEffectEncodings.DOMINANT,
        eff2=sim.SNPEffectEncodings.ADDITIVE,
    add_results_dom_me_pb000.insert(loc=0, column="Encoding", value="Additive")
    add_results_dom_me_pb000.insert(loc=0, column="BioAct", value="Dominant")
    DOMDEV_results_dom_me_pb000.insert(loc=0, column="Encoding", value="DOMDEV")
    DOMDEV_results_dom_me_pb000.insert(loc=0, column="BioAct", value="Dominant")
    rec_results_dom_me_pb000.insert(loc=0, column="Encoding", value="Recessive")
    rec_results_dom_me_pb000.insert(loc=0, column="BioAct", value="Dominant")
    dom_results_dom_me_pb000.insert(loc=0, column="Encoding", value="Dominant")
    dom_results_dom_me_pb000.insert(loc=0, column="BioAct", value="Dominant")
    codom_results_dom_me_pb000.insert(loc=0, column="Encoding", value="Codominant")
    codom_results_dom_me_pb000.insert(loc=0, column="BioAct", value="Dominant")
    edge_results_dom_me_pb000.insert(loc=0, column="Encoding", value="EDGE")
    edge_results_dom_me_pb000.insert(loc=0, column="BioAct", value="Dominant")

    ## Heterozygous Main Effect for SNP1 without interaction
    # Training data
    train_het_main_effect00 = sim.BAMS.from_model(
        eff1=sim.SNPEffectEncodings.HET,
        eff2=sim.SNPEffectEncodings.ADDITIVE,
    edge_weights_het_me.insert(loc=0, column="BioAct", value="Heterozygous")
    test_het_main_effect00 = sim.BAMS.from_model(
        eff1=sim.SNPEffectEncodings.HET,
        eff2=sim.SNPEffectEncodings.ADDITIVE,
    add_results_het_me_pb000.insert(loc=0, column="Encoding", value="Additive")
    add_results_het_me_pb000.insert(loc=0, column="BioAct", value="Heterozygous")
    DOMDEV_results_het_me_pb000.insert(loc=0, column="Encoding", value="DOMDEV")
    DOMDEV_results_het_me_pb000.insert(loc=0, column="BioAct", value="Heterozygous")
    rec_results_het_me_pb000.insert(loc=0, column="Encoding", value="Recessive")
    rec_results_het_me_pb000.insert(loc=0, column="BioAct", value="Heterozygous")
    dom_results_het_me_pb000.insert(loc=0, column="Encoding", value="Dominant")
    dom_results_het_me_pb000.insert(loc=0, column="BioAct", value="Heterozygous")
    codom_results_het_me_pb000.insert(loc=0, column="Encoding", value="Codominant")
    codom_results_het_me_pb000.insert(loc=0, column="BioAct", value="Heterozygous")
    edge_results_het_me_pb000.insert(loc=0, column="Encoding", value="EDGE")
    edge_results_het_me_pb000.insert(loc=0, column="BioAct", value="Heterozygous")


    ## NULL Main Effect for SNP1 without interaction
    # Training data
    train_null_main_effect00 = sim.BAMS.from_model(
        eff1=sim.SNPEffectEncodings.ADDITIVE,
        eff2=sim.SNPEffectEncodings.ADDITIVE,
    edge_weights_null_me = edge_weights_null_me_pb000.copy()
    edge_weights_null_me.insert(loc=0, column="BioAct", value="NULL")

    test_null_main_effect00 = sim.BAMS.from_model(
        eff1=sim.SNPEffectEncodings.ADDITIVE,
        eff2=sim.SNPEffectEncodings.ADDITIVE,
    add_results_null_me_pb000.insert(loc=0, column="Encoding", value="Additive")
    add_results_null_me_pb000.insert(loc=0, column="BioAct", value="NULL")
    DOMDEV_results_null_me_pb000.insert(loc=0, column="Encoding", value="DOMDEV")
    DOMDEV_results_null_me_pb000.insert(loc=0, column="BioAct", value="NULL")
    rec_results_null_me_pb000.insert(loc=0, column="Encoding", value="Recessive")
    rec_results_null_me_pb000.insert(loc=0, column="BioAct", value="NULL")
    dom_results_null_me_pb000.insert(loc=0, column="Encoding", value="Dominant")
    dom_results_null_me_pb000.insert(loc=0, column="BioAct", value="NULL")
    codom_results_null_me_pb000.insert(loc=0, column="Encoding", value="Codominant")
    codom_results_null_me_pb000.insert(loc=0, column="BioAct", value="NULL")
    edge_results_null_me_pb000.insert(loc=0, column="Encoding", value="EDGE")
    edge_results_null_me_pb000.insert(loc=0, column="BioAct", value="NULL")