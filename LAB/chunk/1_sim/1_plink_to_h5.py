import math
import pandas as pd
import numpy as np

from pathlib import Path
from typing import Optional, Union
from pandas_genomics.scalars import MISSING_IDX


def from_plink(input: Union[str, Path], v_chunk: Optional[int] = None):

    # All three files are used
    input = str(input)  # coerce to string in order to add extensions
    bed_file = Path(input + ".bed")
    bim_file = Path(input + ".bim")
    fam_file = Path(input + ".fam")

    # Make sure each file exists
    if not Path(bed_file).exists():
        raise ValueError(f"The .bed file was not found\n\t{str(bed_file)}")
    if not Path(bim_file).exists():
        raise ValueError(f"The .bim file was not found\n\t{str(bim_file)}")
    if not Path(fam_file).exists():
        raise ValueError(f"The .fam file was not found\n\t{str(fam_file)}")

    print(f"Loading genetic data from '{bed_file.stem}'")

    hdf = pd.HDFStore("~/public/data.h5")  # Open .h5 files as Append Mode

    # Load fam file (PLINK sample information file) into a list of Samples
    df = pd.read_table(fam_file, header=None, sep=" ")
    df.columns = ["FID", "IID", "IID_father", "IID_mother", "sex", "phenotype"]
    hdf.put(str("/fam/samples"), df)
    print(f"\tLoaded information for {len(df)} samples from '{fam_file.name}'")

    # Load bim file (PLINK extended MAP file) into a list of variants
    variant_list = pd.read_table(bim_file, header=None, sep="\t")
    # Note 'position' is in centimorgans, 'coordinate' is what pandas-genomics refers to as 'position' (in base-pairs)
    variant_list.columns = [
        "chromosome",
        "variant_id",
        "position",
        "coordinate",
        "allele1",
        "allele2",
    ]
    # hdf.put(str("/bim/variant_list"), variant_list)
    # print(
    #    f"\tLoaded information for {len(variant_list)} variants from '{bim_file.name}'"
    # )

    # Load bed file
    gt_bytes = np.fromfile(bed_file, dtype="uint8")

    # Ensure the file is valid
    CORRECT_FIRST_BYTES = np.array([108, 27, 1], dtype="uint8")
    if not (gt_bytes[:3] == CORRECT_FIRST_BYTES).all():
        raise ValueError(
            f"The first 3 bytes {bed_file.name} were not correct.  The file may be corrupted."
        )
    gt_bytes = gt_bytes[3:]
    # Divide array into one row per variant
    chunk_size = len(df) // 4
    if len(df) % 4 > 0:
        chunk_size += 1
    gt_bytes = gt_bytes.reshape(-1, chunk_size)
    df2 = pd.DataFrame()
    gt_array_dict = {}

    for v_idx, variant in variant_list["variant_id"].items():
        genotypes = np.flip(np.unpackbits(gt_bytes[v_idx]).reshape(-1, 4, 2), axis=1)
        genotypes = genotypes.reshape(-1, 2)[: len(df)]
        missing_gt = (genotypes == (0, 1)).all(axis=1)
        genotypes[missing_gt] = (MISSING_IDX, MISSING_IDX)

        data = np.sum(genotypes, axis=1, dtype=np.uint8)

        gt_array_dict[f"{v_idx}_{variant}"] = data

    df2 = pd.concat([df2, pd.DataFrame.from_dict(gt_array_dict)], axis=1)

    print(f"\tLoading genotypes from '{bed_file.name}' ... ")
    # Quebrar o arquivo
    v_cicles = math.ceil(len(variant_list) / v_chunk)
    vs = 0
    ve = 0
    if v_cicles > 1:
        for i in range(0, v_cicles):
            vs = ve
            ve = ve + v_chunk
            # df2.iloc[:, vs:ve].to_csv("~/Public/fromPlink_result_ciclo_{}.txt".format(i), sep=";")

            hdf.put(str("/bim/variants_" + str(i)), variant_list.iloc[vs:ve, :])

            hdf.put(str("/bed/dataset_" + str(i)), df2.iloc[:, vs:ve])

            print(
                f"\tSaved Genetypes Dataset {i} with {len(df2.iloc[:, vs:ve].columns)} variants"
            )
            hdf.flush()
    else:
        hdf.put(str("/bed/dataset_0"), df2)
        hdf.put(str("/bim/variants_0"), variant_list)
        print(f"\tSaved Genetypes single file")
        hdf.flush()

    hdf.put("/cicles", pd.Series(v_cicles))

    hdf.close()

    return "OK"


if __name__ == "__main__":
    DATA_DIR = Path(__file__).parent / "data" / "plink"
    input = DATA_DIR / "plink_test_medium"
    result = from_plink(input, 5000)
    print("Process finished with status: ", result)
