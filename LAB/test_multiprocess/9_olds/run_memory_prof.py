
# importing the library
from memory_profiler import profile

# instantiating the decorator
@profile
# code for which memory has to
# be monitored
def my_func():
	x = [x for x in range(0, 100000)]
	y = [y*100 for y in range(0, 15000)]
	del x
	return y

if __name__ == '__main__':
	my_func()
