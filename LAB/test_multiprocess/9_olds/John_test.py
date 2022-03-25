from cmd import IDENTCHARS
from operator import ge
import pandas_genomics as pg
import pandas as pd
import clarite
#from memory_profiler import profile

# Define size
num_samples = 1000
num_snps = 10


n_cases = int(num_samples/2)
n_controls = num_samples - n_cases
maf1=0.3
maf2=0.3

# instantiating the decorator
#@profile

def xx():
  # Generate model SNPs
  bams = pg.sim.BAMS.from_model(eff1=pg.sim.SNPEffectEncodings.ADDITIVE, eff2=pg.sim.SNPEffectEncodings.ADDITIVE)
  #print('bams is; ', bams)
  genotypes = bams.generate_case_control(n_cases=n_cases, n_controls=n_controls, maf1=maf1, maf2=maf2)
  genotypes.to_csv('~/Public/gen1.txt', sep=";")

  # Add random SNPs
  random_genotypes = pd.DataFrame({f"var{i}": pg.sim.generate_random_gt(pg.scalars.Variant(ref="A", alt="a"), alt_allele_freq=0.3, n=num_samples, random_seed=i) for i in range(2, num_snps)})
  random_genotypes.to_csv('~/Public/gen2.txt', sep=";")
  genotypes = pd.concat([genotypes, random_genotypes], axis=1)

  # Clarite expects index to be named ID, this is just to avoid a warning
  genotypes.index.name = "ID"
  genotypes.to_csv('~/Public/gen3.txt', sep=";")

  import time
  print("--INICIO DO CLARITE ALL--")
  start = time.process_time()
  results = clarite.analyze.association_study(data=genotypes, outcomes="Outcome", encoding="additive")
  results.to_csv('~/Public/res4.txt', sep=";")
  end = time.process_time()
  print("--FIM DO CLARITE ALL--")


  #Split the independent variables
  """
  print("--INICIO DO CLARITE Split--")
  genotypes1 = genotypes[['Outcome','SNP1','SNP2']]
  genotypes2 = genotypes[['Outcome','var2','var3']]
  results1 = clarite.analyze.association_study(data=genotypes1, outcomes="Outcome", encoding="additive")
  results2 = clarite.analyze.association_study(data=genotypes2, outcomes="Outcome", encoding="additive")
  genotypes1.to_csv('~/Public/gen3_1.txt', sep=";")
  genotypes2.to_csv('~/Public/gen3_2.txt', sep=";")
  results1.to_csv('~/Public/res4_1.txt', sep=";")
  results2.to_csv('~/Public/res4_2.txt', sep=";")
  print("--FIM DO CLARITE Split--")
  

  #Drop Duplicates Rows Samples
  genotypes.drop_duplicates(inplace=True)
  genotypes.to_csv('~/Public/gen3_drop.txt', sep=";")
  results = clarite.analyze.association_study(data=genotypes, outcomes="Outcome", encoding="additive")
  results.to_csv('~/Public/res4_drop.txt', sep=";")
  """

  # Log Elapsed time
  seconds = (end-start)%60
  minutes = (end-start)//60
  if minutes > 60:
      hours = minutes//60
      minutes -= (hours*60)
  else:
      hours = 0
      print(f"Processing {int(num_snps)} SNPs for {int(num_samples)} samples" 
            f" took {int(hours)}:{int(minutes)}:{seconds:.3f}")


if __name__ == '__main__':
	xx()