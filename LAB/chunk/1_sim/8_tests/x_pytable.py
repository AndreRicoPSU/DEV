from tables import *
import numpy as np

# http://mtweb.cs.ucl.ac.uk/mus/martha/PythonPackages/tables-3.0.0rc2/doc/html/usersguide/tutorials.html#browsing-the-object-tree

# PANDAS ==> https://www.numpyninja.com/post/hdf5-file-format-with-pandas


class Particle(IsDescription):
    name = StringCol(16)  # 16-character String
    idnumber = Int64Col()  # Signed 64-bit integer
    ADCcount = UInt16Col()  # Unsigned short integer
    TDCcount = UInt8Col()  # unsigned byte
    grid_i = Int32Col()  # 32-bit integer
    grid_j = Int32Col()  # 32-bit integer
    pressure = Float32Col()  # float  (single-precision)
    energy = Float64Col()  # double (double-precision)


# Creating a PyTables file from scratch
h5file = open_file("tutorial1.h5", mode="w", title="Test file")

# Creating a new group
group = h5file.create_group("/", "detector", "Detector information")

# Creating a new table
table = h5file.create_table(group, "readout", Particle, "Readout example")

# The time has come to fill this table with some values. First we will get a pointer to the Row (see The Row class) instance of this table instance:
particle = table.row

for i in range(10):
    particle["name"] = "Particle: %6d" % (i)
    particle["TDCcount"] = i % 256
    particle["ADCcount"] = (i * 256) % (1 << 16)
    particle["grid_i"] = i
    particle["grid_j"] = 10 - i
    particle["pressure"] = float(i * i)
    particle["energy"] = float(particle["pressure"] ** 4)
    particle["idnumber"] = i * (2**34)
    # Insert a new particle record
    # A call to its append() method writes this information to the table I/O buffer.
    particle.append()

# After we have processed all our data, we should flush the table’s I/O buffer if we want to write all this data to disk. We achieve that by calling the table.flush() method:
table.flush()

print(h5file)


# READING AND SELECTION DATA IN A TABLE

table = h5file.root.detector.readout
pressure = [
    x["pressure"]
    for x in table.iterrows()
    if x["TDCcount"] > 3 and 20 <= x["pressure"] < 50
]
print(pressure)

# talbe.where() is a in-Kernel selection for very large tables and high query speed
names = [
    x["name"]
    for x in table.where("""(TDCcount > 3) & (20 <= pressure) & (pressure < 50)""")
]
print(names)
print(table[6])


# CREATING NEW ARRAY OBJECTS
print("---------------------- \n")
gcolumns = h5file.create_group(h5file.root, "columns", "Pressure and Name")

h5file.create_array(
    gcolumns, "pressure", np.array(pressure), "Pressure column selection"
)
h5file.create_array(gcolumns, "name", names, "Name column selection")

# GUI by vitables ====> GREAT INTERFACE


############################################################################

# BROWSING THE OBJECT TREE


############################################################################

# NESTED


class Info(IsDescription):
    """A sub-structure of Test"""

    _v_pos = 2  # The position in the whole structure
    name = StringCol(10)
    value = Float64Col(pos=0)


colors = Enum(["red", "green", "blue"])


class NestedDescr(IsDescription):
    """A description that has several nested columns"""

    color = EnumCol(colors, "red", base="uint32")
    info1 = Info()

    class info2(IsDescription):
        _v_pos = 1
        name = StringCol(10)
        value = Float64Col(pos=0)

        class info3(IsDescription):
            x = Float64Col(dflt=1)
            y = UInt8Col(dflt=1)


fileh = open_file("nested-tut.h5", "w")
table = fileh.create_table(fileh.root, "table", NestedDescr)

row = table.row
for i in range(10):
    row["color"] = colors[["red", "green", "blue"][i % 3]]
    row["info1/name"] = "name1-%s" % i
    row["info2/name"] = "name2-%s" % i
    row["info2/info3/y"] = i
    # All the rest will be filled with defaults
    row.append()
table.flush()
table.nrows
nra = table[::4]
print(nra)
table.append(nra)
print(table.nrows)
