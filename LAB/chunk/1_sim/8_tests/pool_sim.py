import multiprocessing as mp
import tables as pt


num_arrays = 100
num_processes = mp.cpu_count()
num_simulations = 1000


def Simulation(ii):
    result = []
    result.append(("createGroup", ("/", "A%s" % ii)))
    for i in range(num_arrays):
        result.append(("createArray", ("/A%s" % ii, "B%s" % i, [ii, i])))
    return result


def handle_output(result):
    hdf = pt.open_file("simulation.h5", mode="a")
    for args in result:
        method, args = args
        getattr(hdf, method)(*args)
    hdf.close()


if __name__ == "__main__":
    # clear the file
    hdf = pt.open_file("simulation.h5", mode="w")
    hdf.close()
    pool = mp.Pool(num_processes)
    for i in range(num_simulations):
        pool.apply_async(Simulation, (i,), callback=handle_output)
    pool.close()
    pool.join()
