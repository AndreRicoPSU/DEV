import math
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional, Union


def fam(fam_file, output_file):

    print("----------Started Samples Process---------")
    hdf = pd.HDFStore(output_file)
    df = pd.read_table(fam_file, header=None, sep=" ")
    df.columns = ["FID", "IID", "IID_father", "IID_mother", "sex", "phenotype"]
    v_samples = len(df)
    hdf.put(str("/fam/samples"), df)
    hdf.flush()
    hdf.close()
    print(
        f"\tLoaded information for {v_samples} samples from '{fam_file.name}'")
    return v_samples


def bim(bim_file, output_file, v_chunk):

    print("----------Started Variants Process---------")
    hdf = pd.HDFStore(output_file)
    variant_list = pd.read_table(bim_file, header=None, sep="\t")

    v_cicles = math.ceil(len(variant_list) / v_chunk)

    hdf.put("/cicles", pd.Series(v_cicles))

    variant_splitted = np.array_split(variant_list, v_cicles)
    # variant_splitted = np.asarray(variant_splitted, dtype=object)

    # for idx, x in np.ndenumerate(variant_splitted):
    for x in range(v_cicles):
        base = pd.DataFrame(variant_splitted[x], dtype=object)
        basen = pd.DataFrame({
            "chromosome": base[0].astype(np.int16),
            "variant_id": base[1].astype(np.str0),
            "position": base[2].astype(np.int16),
            "coordinate": base[3].astype(np.int16),
            "allele1": base[4].astype(np.str0),
            "allele2": base[5].astype(np.str0)
        })  # Create this way to avoid warning performance mensage
        dataset = str("/bim/variants_" + str(x))
        hdf.put(dataset, basen)
        hdf.flush()
        print(f"\tAdd datas to dataset {dataset}")
    hdf.close()
    return v_cicles


def bed(bed_file, samples, v_cicles, output_file):
    print("----------Started Genotypes Process---------")
    hdf = pd.HDFStore(output_file)

    gt_bytes = np.fromfile(bed_file, dtype="uint8")

    # Ensure the file is valid
    CORRECT_FIRST_BYTES = np.array([108, 27, 1], dtype="uint8")
    if not (gt_bytes[: 3] == CORRECT_FIRST_BYTES).all():
        raise ValueError(
            f"The first 3 bytes {bed_file.name} were not correct.  The file may be corrupted."
        )
    gt_bytes = gt_bytes[3:]
    # Divide array into one row per variant
    chunk_size = samples // 4
    if samples % 4 > 0:
        chunk_size += 1
    gt_bytes = gt_bytes.reshape(-1, chunk_size)

    gt_bytes_splitted = np.array_split(gt_bytes, v_cicles)

    # for idx, x in np.ndenumerate(variant_splitted):
    for x in range(v_cicles):
        base = pd.DataFrame(gt_bytes_splitted[x], dtype=np.int16)
        dataset = str("/bed/genotype_" + str(x))
        hdf.put(dataset, base)
        hdf.flush()
        print(f"\tAdd datas to dataset {dataset}")

    hdf.close()

    return "OK"


def from_plink(
    input: Union[str, Path], output: Union[str, Path], v_chunk: Optional[int] = None
):

    # All three files are used
    input = str(input)
    bim_file = Path(input + ".bim")
    fam_file = Path(input + ".fam")
    bed_file = Path(input + ".bed")

    output_file = Path(str(output) + ".h5")

    # Make sure each file exists
    if not Path(bim_file).exists():
        raise ValueError(f"The .bim file was not found\n\t{str(bim_file)}")
    if not Path(fam_file).exists():
        raise ValueError(f"The .fam file was not found\n\t{str(fam_file)}")

    print(f"Loading genetic data from '{fam_file.stem}'")

    samples = fam(fam_file, output_file)
    v_cicles = bim(bim_file, output_file, v_chunk)
    genotypes = bed(bed_file, samples, v_cicles, output_file)

    return v_cicles
