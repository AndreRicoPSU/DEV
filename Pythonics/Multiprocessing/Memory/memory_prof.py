"""Memory profiler from PyPI is a python library module used for monitoring 
process memory. It uses psutil code to create a decorator and then uses it to 
get the memory distribution. With this pypi module by importing one can save 
lines and directly call the decorator. To install use the following-

pip install -U memory_profiler
Lets see this through a code-"""

# importing the library
from memory_profiler import profile
import sys

# instantiating the decorator
@profile
# code for which memory has to
# be monitored
def my_func():
	x = [x for x in range(0, 10000)]
	y = [y*100 for y in range(0, 15000)]
	del x
	return y

if __name__ == '__main__':
	my_func()



print(sys.getsizeof(df))
print(df.info(memory_usage="deep"), " S")
print(df.memory_usage(deep=True))

    for name, size in sorted(
        ((name, sys.getsizeof(value)) for name, value in locals().items()),
        key=lambda x: -x[1],
    )[:10]:
        print("{:>30}: {:>8}".format(name, sizeof_fmt(size)))