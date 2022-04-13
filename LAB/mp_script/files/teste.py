import pandas as pd

df = pd.read_csv("teste.csv")
df.to_parquet("teste.parquet")
print(df)
