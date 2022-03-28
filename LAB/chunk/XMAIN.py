
from plink_h5 import *
from h5_df import *

# Enter with the user parameters:
usr_file_name = "plink_test_small"
usr_qtd_variants = 2000
usr_category = True
usr_swap_alleles = False
usr_extract_plink = True
usr_extract_hdf5 = False
usr_data_dir_in = Path(__file__).parent.parent.parent / \
    "Data" / "chunck" / "0_data" / "plink"
usr_data_dir_out = Path(__file__).parent.parent.parent / \
    "Data" / "chunck" / "0_data" / "output"


# Internal Parameters
int_input = usr_data_dir_in / usr_file_name
int_output = usr_data_dir_out / usr_file_name


if __name__ == "__main__":

    # Convert PLINK 1.9 to H5
    if usr_extract_plink:
        result = from_plink(int_input, int_output, usr_qtd_variants)
        print("--->>>  HDF5 file created success with", result, "datasets")

    # Convert H5 to DataFrame
    if usr_extract_hdf5:
        v_db_cicle = get_cicles(int_output)
        print(v_db_cicle)

        for v_cicle in range(v_db_cicle):
            v_db_genomic = get_db(int_output, v_cicle,
                                  usr_category, usr_swap_alleles)

            # Here we have a DataFrame with n Variants with all samples
            # Add the sprint code

        print("loop done")

    print("done")
