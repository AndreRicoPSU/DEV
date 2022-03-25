import pandas as pd
import numpy as np
import clarite
import time
start = time.time()
from pandas_genomics import sim, io, scalars, GenotypeDtype
from memory_profiler import profile  

# Define the variant for two SNPs
from pandas_genomics.scalars import Variant
variant1 = Variant('1', 1, id='rs1', ref='A', alt=['C'])
variant2 = Variant('1', 2, id='rs2', ref='G', alt=['T'])

# Define Case-Control ratio
num_samples = 5000
n_cases = int(num_samples/4)
n_controls = num_samples - n_cases

EDGE_alpha_Final = pd.DataFrame()
All_Results_Final = pd.DataFrame()

@profile
def xx():
    for train_seed, test_seed in zip(range(0,2),range(2000,2)): #1000
        startcycle = time.time()
        ## Recessive Main Effect for SNP1 without interaction 
        # Training data
        train_rec_main_effect465 = sim.BAMS.from_model(eff1=sim.SNPEffectEncodings.RECESSIVE, eff2=sim.SNPEffectEncodings.ADDITIVE, 
            penetrance_base=0.465, main1=1, main2=0, interaction=0, snp1=variant1, snp2=variant2, random_seed=train_seed)
        train_rec_me_pb0465 = train_rec_main_effect465.generate_case_control(n_cases=n_cases, n_controls=n_controls)
        
        # Calculate weights from the training dataset 
        edge_weights_rec_me_pb0465 = train_rec_me_pb0465.genomics.calculate_edge_encoding_values(data=train_rec_me_pb0465["Outcome"], outcome_variable="Outcome")
        edge_weights_rec_me = edge_weights_rec_me_pb0465.copy()
        edge_weights_rec_me.insert(loc=0, column='BioAct', value="Recessive")
        edge_weights_rec_me.insert(loc=0, column='TrainSeed', value=train_seed)
        edge_weights_rec_me.insert(loc=0, column='TestSeed', value=test_seed)
        
        # Test data
        test_rec_main_effect465 = sim.BAMS.from_model(eff1=sim.SNPEffectEncodings.RECESSIVE, eff2=sim.SNPEffectEncodings.ADDITIVE, penetrance_base=0.465, main1=1, main2=0, interaction=0, snp1=variant1, snp2=variant2, random_seed=test_seed)
        test_rec_me_pb0465 = test_rec_main_effect465.generate_case_control(n_cases=n_cases, n_controls=n_controls)
        test_rec_me_pb0465['Outcome'].cat.reorder_categories(['Control', 'Case'], inplace=True)
        
        # Run Regression by using weightes from CLARITE
        # EDGE Encoding
        test_rec_me_pb0465_CLARITE=test_rec_me_pb0465.genomics.encode_edge(encoding_info=edge_weights_rec_me_pb0465)
        
        edge_results_rec_me_pb0465=clarite.analyze.association_study(data=test_rec_me_pb0465_CLARITE, outcomes="Outcome")
        edge_results_rec_me_pb0465['odds ratio'] = np.exp(edge_results_rec_me_pb0465['Beta'])
        edge_results_rec_me_pb0465.insert(loc=0, column='Encoding', value="EDGE")
        edge_results_rec_me_pb0465.insert(loc=0, column='BioAct', value="Recessive")
        edge_results_rec_me_pb0465.insert(loc=0, column='TrainSeed', value=train_seed)
        edge_results_rec_me_pb0465.insert(loc=0, column='TestSeed', value=test_seed)
        




        ## Additive Main Effect for SNP1 without interaction 
        # Training data
        train_add_main_effect465 = sim.BAMS.from_model(eff1=sim.SNPEffectEncodings.ADDITIVE, eff2=sim.SNPEffectEncodings.ADDITIVE, penetrance_base=0.465, main1=1, main2=0, interaction=0, snp1=variant1, snp2=variant2, random_seed=train_seed)
        train_add_me_pb0465 = train_add_main_effect465.generate_case_control(n_cases=n_cases, n_controls=n_controls)
        
        # Calculate weights from the training dataset 
        edge_weights_add_me_pb0465 = train_add_me_pb0465.genomics.calculate_edge_encoding_values(data=train_add_me_pb0465["Outcome"], outcome_variable="Outcome")
        edge_weights_add_me = edge_weights_add_me_pb0465.copy()
        edge_weights_add_me.insert(loc=0, column='BioAct', value="Additive")
        edge_weights_add_me.insert(loc=0, column='TrainSeed', value=train_seed)
        edge_weights_add_me.insert(loc=0, column='TestSeed', value=test_seed)
        
        # Test data
        test_add_main_effect465 = sim.BAMS.from_model(eff1=sim.SNPEffectEncodings.ADDITIVE, eff2=sim.SNPEffectEncodings.ADDITIVE, penetrance_base=0.465, main1=1, main2=0, interaction=0, snp1=variant1, snp2=variant2, random_seed=test_seed)
        test_add_me_pb0465 = test_add_main_effect465.generate_case_control(n_cases=n_cases, n_controls=n_controls)
        test_add_me_pb0465['Outcome'].cat.reorder_categories(['Control', 'Case'], inplace=True)
        
        # Run Regression by using weightes from CLARITE
        # EDGE Encoding
        test_add_me_pb0465_CLARITE=test_add_me_pb0465.genomics.encode_edge(encoding_info=edge_weights_add_me_pb0465)
        
        edge_results_add_me_pb0465=clarite.analyze.association_study(data=test_add_me_pb0465_CLARITE, outcomes="Outcome")
        edge_results_add_me_pb0465['odds ratio'] = np.exp(edge_results_add_me_pb0465['Beta'])
        edge_results_add_me_pb0465.insert(loc=0, column='Encoding', value="EDGE")
        edge_results_add_me_pb0465.insert(loc=0, column='BioAct', value="Additive")
        edge_results_add_me_pb0465.insert(loc=0, column='TrainSeed', value=train_seed)
        edge_results_add_me_pb0465.insert(loc=0, column='TestSeed', value=test_seed)
        


        ## Dominant Main Effect for SNP1 without interaction 
        # Training data
        train_dom_main_effect465 = sim.BAMS.from_model(eff1=sim.SNPEffectEncodings.DOMINANT, eff2=sim.SNPEffectEncodings.ADDITIVE, penetrance_base=0.465, main1=1, main2=0, interaction=0, snp1=variant1, snp2=variant2, random_seed=train_seed)
        train_dom_me_pb0465 = train_dom_main_effect465.generate_case_control(n_cases=n_cases, n_controls=n_controls)
        
        # Calculate weights from the training dataset 
        edge_weights_dom_me_pb0465 = train_dom_me_pb0465.genomics.calculate_edge_encoding_values(data=train_dom_me_pb0465["Outcome"], outcome_variable="Outcome")
        edge_weights_dom_me = edge_weights_dom_me_pb0465.copy()
        edge_weights_dom_me.insert(loc=0, column='BioAct', value="Dominant")
        edge_weights_dom_me.insert(loc=0, column='TrainSeed', value=train_seed)
        edge_weights_dom_me.insert(loc=0, column='TestSeed', value=test_seed)
        
        # Test data
        test_dom_main_effect465 = sim.BAMS.from_model(eff1=sim.SNPEffectEncodings.DOMINANT, eff2=sim.SNPEffectEncodings.ADDITIVE, penetrance_base=0.465, main1=1, main2=0, interaction=0, snp1=variant1, snp2=variant2, random_seed=test_seed)
        test_dom_me_pb0465 = test_dom_main_effect465.generate_case_control(n_cases=n_cases, n_controls=n_controls)
        test_dom_me_pb0465['Outcome'].cat.reorder_categories(['Control', 'Case'], inplace=True)
        
        # Run Regression by using weightes from CLARITE
        # EDGE Encoding
        test_dom_me_pb0465_CLARITE=test_dom_me_pb0465.genomics.encode_edge(encoding_info=edge_weights_dom_me_pb0465)
        
        edge_results_dom_me_pb0465=clarite.analyze.association_study(data=test_dom_me_pb0465_CLARITE, outcomes="Outcome")
        edge_results_dom_me_pb0465['odds ratio'] = np.exp(edge_results_dom_me_pb0465['Beta'])
        edge_results_dom_me_pb0465.insert(loc=0, column='Encoding', value="EDGE")
        edge_results_dom_me_pb0465.insert(loc=0, column='BioAct', value="Dominant")
        edge_results_dom_me_pb0465.insert(loc=0, column='TrainSeed', value=train_seed)
        edge_results_dom_me_pb0465.insert(loc=0, column='TestSeed', value=test_seed)
        


        ## Heterozygous Main Effect for SNP1 without interaction 
        # Training data
        train_het_main_effect465 = sim.BAMS.from_model(eff1=sim.SNPEffectEncodings.HET, eff2=sim.SNPEffectEncodings.ADDITIVE, penetrance_base=0.465, main1=1, main2=0, interaction=0, snp1=variant1, snp2=variant2, random_seed=train_seed)
        train_het_me_pb0465 = train_het_main_effect465.generate_case_control(n_cases=n_cases, n_controls=n_controls)
        
        # Calculate weights from the training dataset 
        edge_weights_het_me_pb0465 = train_het_me_pb0465.genomics.calculate_edge_encoding_values(data=train_het_me_pb0465["Outcome"], outcome_variable="Outcome")
        edge_weights_het_me = edge_weights_het_me_pb0465.copy()
        edge_weights_het_me.insert(loc=0, column='BioAct', value="Heterozygous")
        edge_weights_het_me.insert(loc=0, column='TrainSeed', value=train_seed)
        edge_weights_het_me.insert(loc=0, column='TestSeed', value=test_seed)
        
        # Test data
        test_het_main_effect465 = sim.BAMS.from_model(eff1=sim.SNPEffectEncodings.HET, eff2=sim.SNPEffectEncodings.ADDITIVE, penetrance_base=0.465, main1=1, main2=0, interaction=0, snp1=variant1, snp2=variant2, random_seed=test_seed)
        test_het_me_pb0465 = test_het_main_effect465.generate_case_control(n_cases=n_cases, n_controls=n_controls)
        test_het_me_pb0465['Outcome'].cat.reorder_categories(['Control', 'Case'], inplace=True)
        
        # Run Regression by using weightes from CLARITE
        # EDGE Encoding
        test_het_me_pb0465_CLARITE=test_het_me_pb0465.genomics.encode_edge(encoding_info=edge_weights_het_me_pb0465)
        
        edge_results_het_me_pb0465=clarite.analyze.association_study(data=test_het_me_pb0465_CLARITE, outcomes="Outcome")
        edge_results_het_me_pb0465['odds ratio'] = np.exp(edge_results_het_me_pb0465['Beta'])
        edge_results_het_me_pb0465.insert(loc=0, column='Encoding', value="EDGE")
        edge_results_het_me_pb0465.insert(loc=0, column='BioAct', value="Heterozygous")
        edge_results_het_me_pb0465.insert(loc=0, column='TrainSeed', value=train_seed)
        edge_results_het_me_pb0465.insert(loc=0, column='TestSeed', value=test_seed)
        



        ##Concat the results
        EDGE_alpha_Results=pd.concat([edge_weights_rec_me,edge_weights_add_me,edge_weights_dom_me,edge_weights_het_me])
        EDGE_alpha_Final=pd.concat([EDGE_alpha_Final,EDGE_alpha_Results], axis = 0)
        #EDGE_alpha_Final.to_csv('/storage/home/jpz5091/work/bams/SimPower50000_PB0465/EDGE_alpha_Results_5000_pb0.465_case1control3.txt', sep=";")
        
        All_Results=pd.concat([edge_results_rec_me_pb0465,edge_results_add_me_pb0465,edge_results_dom_me_pb0465,edge_results_het_me_pb0465])
        All_Results_Final=pd.concat([All_Results_Final,All_Results], axis = 0)
        #All_Results_Final.to_csv('/storage/home/jpz5091/work/bams/SimPower50000_PB0465/All_Results_5000_pb0.465_case1control3.txt', sep=";")
        
        endcycle = time.time()
        print("The time of execution of one cycle is :", endcycle-startcycle)


    end = time.time()
    print("The time of execution of above program is :", end-start)

if __name__ == '__main__':
    xx()