from matplotlib.image import BboxImage
import pandas as pd
import numpy as np
from pandas_genomics.scalars import Variant, MISSING_IDX
from pandas_genomics.arrays import GenotypeDtype, GenotypeArray


def create_variant(variant_info_row):
    variant_id = str(variant_info_row["variant_id"])
    a1 = str(variant_info_row["allele1"])
    a2 = str(variant_info_row["allele2"])
    # 0 indicates a missing allele
    if a2 == "0":
        a2 = None
    if a1 == "0":
        a1 = None
    else:
        a1 = [a1]  # pass as list
    # Ensure chromosome is None instead of nan
    if np.isnan(variant_info_row["chromosome"]):
        chromosome = None
    else:
        chromosome = str(variant_info_row["chromosome"])
    variant = Variant(
        chromosome=chromosome,
        position=int(variant_info_row["coordinate"]),
        id=variant_id,
        ref=a2,
        alt=a1,
        ploidy=2,
    )
    return variant


hdf = pd.HDFStore("~/Public/data.h5", "r")

print(hdf.info())  # Print structure of the file

# for key in hdf.keys():
#    print("the key in PD is: ", key)
#    v_bed = hdf.get(key)

# Tratar a Dataset de samples
df_samples = hdf.get("/fam/samples")
num_samples = len(df_samples)

# 1. Create parameter to choise change sex and phenotype

if 1 == 1:  # Update 'sex'
    df_samples["sex"] = df_samples["sex"].astype("category")
    df_samples["sex"] = df_samples["sex"].cat.rename_categories(
        {1: "male", 2: "female", 0: "unknown"}
    )

if 1 == 1:  # categorical_phenotype:
    DEFAULT_CAT_MAP = {1: "Control", 2: "Case"}
    df_samples["phenotype"] = df_samples["phenotype"].astype("category")
    df_samples["phenotype"].cat.rename_categories(DEFAULT_CAT_MAP, inplace=True)
    df_samples.loc[
        ~df_samples["phenotype"].isin(DEFAULT_CAT_MAP.values()), "phenotype"
    ] = None


v_cicles = 0


# Tratar a Dataset de Variante List

df_variant = hdf.get(str("/bim/variants_" + str(v_cicles)))

df_variant["chromosome"] = df_variant["chromosome"].astype("category")

variant_list = [create_variant(row) for idx, row in df_variant.iterrows()]
print(f"\tLoaded information for {len(variant_list)} variants from '{hdf}'")


# trtar o bed
genotypes = hdf.get(str("/bed/dataset_" + str(v_cicles)))

# Process each variant
gt_array_dict = {}

for x, variant in enumerate(variant_list):

    # transform int to [1, 1]
    genotypes[genotypes.columns[x]] = np.where(
        (genotypes[genotypes.columns[x]] == 1)[:, None], [0, 1], [1, 1]
    ).tolist()

    dtype = GenotypeDtype(variant)

    scores = np.ones(num_samples) * MISSING_IDX  # Missing Scores

    data = np.array(
        list(zip(genotypes[genotypes.columns[x]], scores)), dtype=dtype._record_type
    )

    gt_array = GenotypeArray(values=data, dtype=dtype)

    gt_array_dict[f"{x}_{gt_array.variant.id}"] = gt_array

df = pd.concat([df_samples, pd.DataFrame.from_dict(gt_array_dict)], axis=1)
df = df.set_index(["FID", "IID", "IID_father", "IID_mother", "sex", "phenotype"])

print(df)

df.to_csv("~/Public/teste.csv")
