
from plink_h5 import *
from h5_df import *
import clarite

# Enter with the user parameters:
usr_file_name = "sim_plink_big"
usr_qtd_variants = 1000
usr_category = True
usr_swap_alleles = False
usr_extract_plink = False
usr_extract_hdf5 = True
usr_data_dir_in = Path(__file__).parent.parent / "plink"
usr_data_dir_out = Path(__file__).parent.parent / "output"


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

        # for v_cicle in range(v_db_cicle):
        for v_cicle in range(1):
            v_db_genomic = get_db(int_output, v_cicle,
                                  usr_category, usr_swap_alleles)

            print("Cicle", v_cicle)
            v_db_genomic.to_clipboard

            # Here we have a DataFrame with n Variants with all samples
            # v_db_genomic.reset_index
            # v_db_genomic.droplevel(level=0)
            print(v_db_genomic)

            #v_db_genomic.drop("FID", axis=1, inplace=True)

            #df = v_db_genomic.drop(v_db_genomic.columns[0], axis=1)
            #df.drop('column_name', axis=1, inplace=True)
            result_clarite = clarite.analyze.association_study(
                data=v_db_genomic, outcomes="phenotype")
            # Add the sprint code
            # Test with CLARITE call
            # Add the outcome results file or return to .h5

            result_clarite.to_csv(
                f"{usr_data_dir_out}_clarite_{v_cicle}.txt", sep=";")

            #v_db_genomic["Outcome"].cat.reorder_categories(["Control", "Case"], inplace=True)

            # Run Regression by using weightes from CLARITE
            # Addtive Encoding
            #test_rec_me_pb000_ADD = test_rec_me_pb000.genomics.encode_additive()

        print("--->>>  HDF5 file was read success with",
              v_db_cicle, "datasets reads")

    print("done")
