
from transf_plink_h5 import *
from transf_h5_pg import *
import clarite


# Enter the steps will be performed:
usr_extract_plink = False  # Read the plink files and create the .h5 file
# Read the .h5, convert in panda-genomics and run CUSTOM SCRIPT
usr_extract_hdf5 = True

# Enter with the user parameters:
usr_qtd_variants = 1000  # Qtd of variants in each dataset to save in .h5
usr_cicle_run = 40  # If None will extract all dataset in .h5 to panda-genomics
usr_category = True
usr_swap_alleles = False

# Enter with the file name and location:
usr_file_name = "plink_test_medium"
int_input = Path(__file__).parent.parent / "plink" / \
    usr_file_name  # Folder Path to plink
int_output = Path(__file__).parent.parent / "output" / \
    usr_file_name  # Folder Path to .h5


def cust_script(df_genomics, v_cicle):
    output_script = Path(__file__).parent.parent / "output" / "file"
    result_clarite = clarite.analyze.association_study(
        data=df_genomics, outcomes="phenotype")
    result_clarite.to_csv(f"{output_script}_clarite_{v_cicle}.txt", sep=";")

    return "ok"


if __name__ == "__main__":

    # Convert plink 1.9 to .h5
    if usr_extract_plink:
        print("> star process of converting plink to .h5")
        result = from_plink(int_input, int_output, usr_qtd_variants)
        print("> end process of conversion plink to .h5 in", result, "datasets")

    # Convert .h5 to pandas-genomics
    if usr_extract_hdf5:

        if usr_cicle_run == None:
            v_db_cicle = get_cicles(int_output)
            print(">  start process to converting .h5 to pandas-genomics with:",
                  v_db_cicle, "datasets")

            for v_cicle in range(v_db_cicle):
                print(">   start extract dataset",
                      v_cicle, "to pandas-genomics")
                df_genomics = get_db(int_output, v_cicle,
                                     usr_category, usr_swap_alleles)
                print(">   end extract dataset",
                      v_cicle, "to pandas-genomics \n")

                # CUSTOM SCRIPT
                cust_script(df_genomics, v_cicle)

        else:
            v_db_cicle = usr_cicle_run
            print(">  start extract dataset",
                  v_db_cicle, "to pandas-genomics")
            df_genomics = get_db(int_output, v_db_cicle,
                                 usr_category, usr_swap_alleles)
            print(">  end extract dataset",
                  v_db_cicle, "to pandas-genomics")

            # CUSTOM SCRIPT
            cust_script(df_genomics, v_db_cicle)


print(">END OF PROCESS")
