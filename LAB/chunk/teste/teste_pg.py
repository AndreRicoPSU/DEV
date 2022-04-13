import pandas_genomics

# teste = pandas_genomics.io.from_plink(
#    "plink_test_small")
# print(teste)

my_1 = pandas_genomics.scalars.Variant(
    '22', 434434343, 'rs007', 'D', 'T',)
my_2 = pandas_genomics.scalars.Variant(
    '22', 434434343, 'rs007', 'D', 'A', ploidy=3)

if my_1 == my_2:
    print("same")
else:
    print("not the same")

my_1.add_allele("C")
my_2.add_allele("R")

print(my_1)


my_g = my_1.make_genotype("D", "T")
print(my_g)

my_f = pandas_genomics.scalars.Genotype(my_2, (0, 1, 2))
print(my_f)

print("-----")
v = pandas_genomics.arrays.GenotypeDtype(my_2)

print(v)

a = v.
print(a)
