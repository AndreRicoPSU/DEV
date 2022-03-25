"""
ARGPARSE
--------
Help to input arguments when we call a program
ex:
$ Python3 main.py teste
Link Doc: https://docs.python.org/3/howto/argparse.html#combining-positional-and-optional-arguments

"""


"""
#exemplo 1
import argparse
parser = argparse.ArgumentParser() # Cria objeto
parser.add_argument("echo", help="esse e um teste") # adiciona 
args = parser.parse_args() # recebe
print(args.echo) 


#Exemplo 2
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("square", type=int,
                    help="display a square of a given number")
parser.add_argument("-v", "--verbose", action="store_true",
                    help="increase output verbosity")
args = parser.parse_args()
answer = args.square**2
if args.verbose:
    print(f"the square of {args.square} equals {answer}")
else:
    print(answer)
"""

#Exemplo 3
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("square", type=int,
                    help="display a square of a given number")
parser.add_argument("-v", "--verbosity", type=int,
                    help="increase output verbosity")
args = parser.parse_args()
answer = args.square**2
if args.verbosity == 2:
    print(f"the square of {args.square} equals {answer}")
elif args.verbosity == 1:
    print(f"{args.square}^2 == {answer}")
else:
    print(answer)