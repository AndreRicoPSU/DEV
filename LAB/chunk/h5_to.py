from matplotlib.image import BboxImage
import pandas as pd
import numpy as np
from pandas_genomics.scalars import Variant, MISSING_IDX
from pandas_genomics.arrays import GenotypeDtype, GenotypeArray
from pathlib import Path
from typing import Optional, Union

from psutil import swap_memory


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


def xxx(
    input: Union[str, Path],
    output: Union[str, Path],
    db_cicle,
    v_category,
    swap_alleles,
):
    input = str(input)
    h5_file = Path(input + ".h5")
    hdf = pd.HDFStore(h5_file, "r")
    print(hdf.info())  # Print structure of the file

    # SAMPLES - FROM FAM FOLDER/
    df_samples = hdf.get("/fam/samples")
    # num_samples = len(df_samples)

    if v_category:  # categorical for sex
        df_samples["sex"] = df_samples["sex"].astype("category")
        df_samples["sex"] = df_samples["sex"].cat.rename_categories(
            {1: "male", 2: "female", 0: "unknown"}
        )

    if v_category:  # categorical for phenotype
        DEFAULT_CAT_MAP = {1: "Control", 2: "Case"}
        df_samples["phenotype"] = df_samples["phenotype"].astype("category")
        df_samples["phenotype"].cat.rename_categories(DEFAULT_CAT_MAP, inplace=True)
        df_samples.loc[
            ~df_samples["phenotype"].isin(DEFAULT_CAT_MAP.values()), "phenotype"
        ] = None

    # VARIANTS - FROM BIM [HDF/BIM/VARIANTS_X]

    df_variant = hdf.get(str("/bim/variants_" + str(db_cicle)))

    df_variant.columns = [
        "chromosome",
        "variant_id",
        "position",
        "coordinate",
        "allele1",
        "allele2",
    ]

    df_variant["chromosome"] = df_variant["chromosome"].astype("category")

    variant_list = [create_variant(row) for idx, row in df_variant.iterrows()]

    print(f"\tLoaded information for {len(variant_list)} variants from '{hdf}'")

    # GENOTYPES - FROM BED [HDF/BED/GENOTYPE_X]

    gt_bytes_df = hdf.get(str("/bed/genotype_" + str(db_cicle)))
    gt_bytes = gt_bytes_df.to_numpy().astype(
        "uint8"
    )  # improve this point with pytables
    # gt_bytes = gt_bytes.astype("uint8")

    gt_array_dict = {}

    # for v_idx, variant in variant_list["variant_id"].items():
    for v_idx, variant in enumerate(variant_list):
        test = np.unpackbits(gt_bytes[v_idx])
        genotypes = np.flip(np.unpackbits(gt_bytes[v_idx]).reshape(-1, 4, 2), axis=1)
        genotypes = genotypes.reshape(-1, 2)[: len(df_samples)]
        missing_gt = (genotypes == (0, 1)).all(axis=1)
        genotypes[missing_gt] = (MISSING_IDX, MISSING_IDX)
        het_gt = (genotypes == (1, 0)).all(axis=1)
        genotypes[het_gt] = (0, 1)
        dtype = GenotypeDtype(variant)
        scores = np.ones(len(df_samples)) * MISSING_IDX  # Missing Scores
        data = np.array(list(zip(genotypes, scores)), dtype=dtype._record_type)
        gt_array = GenotypeArray(values=data, dtype=dtype)

        if swap_alleles:
            gt_array.set_reference(1)
        gt_array_dict[
            f"{v_idx}_{gt_array.variant.id}"
        ] = gt_array  # maybe alter to variant_id

    df_samples = pd.concat(
        [df_samples, pd.DataFrame.from_dict(gt_array_dict)], axis=1
    )  # maybe isn't necessary

    df_samples = df_samples.set_index(
        ["FID", "IID", "IID_father", "IID_mother", "sex", "phenotype"]
    )

    return df_samples


# print(f"\tLoaded information for {len(variant_list)} variants from '{hdf}'")


if __name__ == "__main__":
    DATA_DIR_IN = Path(__file__).parent / "0_data" / "output"
    DATA_DIR_OUT = Path(__file__).parent / "0_data" / "output"
    input = DATA_DIR_IN / "plink_test_small"
    output = DATA_DIR_OUT / "plink_test_small"
    db_cicle = 0
    v_category = True
    swap_alleles = False

    v_temp = xxx(input, output, db_cicle, v_category, swap_alleles)
    print("done")
