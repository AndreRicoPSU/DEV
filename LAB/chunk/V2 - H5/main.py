
from transf_plink_h5 import *
from transf_h5_pg import *
import clarite
import argparse


# Enter with the user parameters:
usr_qtd_variants = 1000  # Qtd of variants in each dataset to save in .h5
usr_category = True
usr_swap_alleles = False

# Enter with the file name and location:
usr_file_name = "plink_test_medium"
int_input = "/Users/andrerico/DEV/LAB/chunk/plink/" + usr_file_name
int_output = "/Users/andrerico/DEV/LAB/chunk/output/" + usr_file_name
#int_input = "/gpfs/group/mah546/default/datasets/simulated_genotypes/" + usr_file_name
#int_output = "/gpfs/group/mah546/default/datasets/simulated_genotypes/" + usr_file_name


def cust_script(df_genomics, v_cicle):
    #output_script = "/Users/andrerico/DEV/LAB/chunk/output/file"
    output_script = "/gpfs/group/mah546/default/projects/edge/HDF5_convert/output/result"
    result_clarite = clarite.analyze.association_study(
        data=df_genomics, outcomes="phenotype")
    result_clarite.to_csv(f"{output_script}_clarite_{v_cicle}.txt", sep=";")

    return "ok"


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--start", type=int,
                        metavar="dataset",
                        action="store",
                        default=None,
                        help="Inform with dataset will start",
                        )
    parser.add_argument("-e", "--end", type=int,
                        metavar="dataset",
                        action="store",
                        default=None,
                        help="Inform with dataset will end",
                        )

    parser.add_argument("-p", "--extract_plink", action="store_true", help="")
    parser.add_argument("-r", "--extract_hdf5", action="store_true", help="")

    args = parser.parse_args()

    #args.extract_plink = True
    #args.extract_hdf5 = True
    #args.start = 10
    #args.end = 19

    # Convert plink 1.9 to .h5
    # if usr_extract_plink:
    if args.extract_plink:
        print("> star process of converting plink to .h5")
        result = from_plink(int_input, int_output, usr_qtd_variants)
        print("> end process of conversion plink to .h5 in", result, "datasets")

    # Convert .h5 to pandas-genomics
    # if usr_extract_hdf5:
    if args.extract_hdf5:
        if args.start == None and args.end == None:
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

        elif args.start != None and args.end != None:
            print(">  start process to converting .h5 to pandas-genomics start:",
                  args.start, "dataset and end:", args.end, "dataset")

            args.end += 1
            for v_cicle in range(args.start, args.end):

                print(">  start extract dataset",
                      v_cicle, "to pandas-genomics")
                df_genomics = get_db(int_output, v_cicle,
                                     usr_category, usr_swap_alleles)
                print(">  end extract dataset",
                      v_cicle, "to pandas-genomics")

                # CUSTOM SCRIPT
                cust_script(df_genomics, v_cicle)

        else:
            print("> Informa start and end dataset to process")

print(">END OF PROCESS")
