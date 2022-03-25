#import multiprocessing
import time
from tracemalloc import start
import concurrent.futures
# importing the library
from memory_profiler import profile

# instantiating the decorator


start = time.perf_counter()
@profile
def do_something(second):
    print(f'Sleeping {second} second...')
    time.sleep(second)
    x = f'Done Sleeping...{second}'
    y = f'teste...{second}'
    return x, y

if __name__ == '__main__':

    with concurrent.futures.ProcessPoolExecutor() as executor:
        #f1 = executor.submit(do_something, 1)
        #print(f1.result())
        secs = [5,4,3,2,1]
        #results = [executor.submit(do_something, sec) for sec in secs]
        #for f in concurrent.futures.as_completed(results):
        #    print(f.result())

        results = executor.map(do_something, secs)
        for x, y in results:
            print(x)
            print(y)

    #processes = []

    # with Process
    #for _ in range(10):
    #    p = multiprocessing.Process(target=do_something)
    #    p.start()
    #    processes.append(p)
    #
    #for process in processes:
    #    process.join()
    #
    finish = time.perf_counter()

    print(f"Finished in {round(finish-start,2)} second(s)")