import pandas as pd
import numpy as np
import clarite
import time


start = time.time()
from pandas_genomics import sim, io, scalars, GenotypeDtype
from multiprocessing import Pool
from pandas_genomics.scalars import Variant


# Name conversation for class and methods parameters
def conv_group_eff(i):
    switcher={
        'REC':'RECESSIVE',
        'ADD':'ADDITIVE',
        'DOM':'DOMINANT',
        'HED':'HET'
        }
    return switcher.get(i,"Invalid Group")
def conv_group_edge(i):
    switcher={
        'REC':'Recessive',
        'ADD':'Additive',
        'DOM':'Dominant',
        'HED':'Heterozygous'
        }
    return switcher.get(i,"Invalid Group")

def worker_process(work_group_data):
    #Normalize the nomenclature
    eff1 = conv_group_eff(work_group_data[0])
    eff2 = conv_group_eff(work_group_data[1])
    edge = conv_group_edge(work_group_data[0])

    #print(f"inside worker: {eff1} and {eff2} and {edge}")
  
    # Training data
    train_main = sim.BAMS.from_model(eff1=getattr(sim.SNPEffectEncodings, eff1), eff2=getattr(sim.SNPEffectEncodings, eff2),penetrance_base=0.465, main1=1, main2=0, interaction=0, snp1=variant1, snp2=variant2, random_seed=train_seed)
    train_me = train_main.generate_case_control(n_cases=n_cases, n_controls=n_controls)

    # Calculate weights from the training dataset 
    edge_weights_me = train_me.genomics.calculate_edge_encoding_values(data=train_me["Outcome"], outcome_variable="Outcome")
    edge_weights_me = edge_weights_me.copy()
    edge_weights_me.insert(loc=0, column='BioAct', value=edge)
    edge_weights_me.insert(loc=0, column='TrainSeed', value=train_seed)
    edge_weights_me.insert(loc=0, column='TestSeed', value=test_seed)

    # Test data
    test_main = sim.BAMS.from_model(eff1=getattr(sim.SNPEffectEncodings, eff1), eff2=getattr(sim.SNPEffectEncodings, eff2), penetrance_base=0.465, main1=1, main2=0, interaction=0, snp1=variant1, snp2=variant2, random_seed=test_seed)
    test_me = test_main.generate_case_control(n_cases=n_cases, n_controls=n_controls)
    test_me['Outcome'].cat.reorder_categories(['Control', 'Case'], inplace=True)


    # Run Regression by using weightes from CLARITE
    # EDGE Encoding
    test_me_CLARITE=test_me.genomics.encode_edge(encoding_info=edge_weights_me)
   
    edge_results_me=clarite.analyze.association_study(data=test_me_CLARITE, outcomes="Outcome")
    edge_results_me['odds ratio'] = np.exp(edge_results_me['Beta'])
    edge_results_me.insert(loc=0, column='Encoding', value="EDGE")
    edge_results_me.insert(loc=0, column='BioAct', value="Additive")
    edge_results_me.insert(loc=0, column='TrainSeed', value=train_seed)
    edge_results_me.insert(loc=0, column='TestSeed', value=test_seed)

    ##Concat the results
    EDGE_alpha_Final=pd.concat([EDGE_alpha_Final,edge_weights_me], axis = 0)
    All_Results_Final=pd.concat([All_Results_Final,edge_results_me], axis = 0)

def poolHandler():
    p = Pool(1)
    p.map(worker_process, work_group)


if __name__ == '__main__':
  
    ## ----------------------------------------------------- ##
    ## START: SESSION TO DEFINE EXECUTION PARAMETERS         ##
    
    # Define the variant for two SNPs
    variant1 = Variant('1', 1, id='rs1', ref='A', alt=['C'])
    variant2 = Variant('1', 2, id='rs2', ref='G', alt=['T'])
    
    # Define Case-Control ratio
    num_samples = 200 #5000
    n_cases = int(num_samples/4)
    n_controls = num_samples - n_cases

    EDGE_alpha_Final = pd.DataFrame()
    All_Results_Final = pd.DataFrame()

    # Groups of Eff1 and Eff2
    work_group = (['REC','ADD'],['ADD','ADD'],['DOM','ADD'],['HED','ADD'])

    ## END: SESSION TO DEFINE EXECUTION PARAMETERS           ##


    for train_seed, test_seed in zip(range(0,1),range(2000,2001)):
        startcycle = time.time()

        poolHandler()

        ##Concat the results
        #EDGE_alpha_Final.to_csv('/storage/home/jpz5091/work/bams/SimPower50000_PB0465/EDGE_alpha_Results_5000_pb0.465_case1control3.txt', sep=";")
        #All_Results_Final.to_csv('/storage/home/jpz5091/work/bams/SimPower50000_PB0465/All_Results_5000_pb0.465_case1control3.txt', sep=";")
        
        endcycle = time.time()
        print("The time of execution of one cycle is :", endcycle-startcycle)


    end = time.time()
    print("The time of execution of above program is :", end-start)