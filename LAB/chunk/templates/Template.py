from datetime import timedelta
import pandas as pd
import numpy as np
from scipy import randn

# Source: https://pandas.pydata.org/pandas-docs/dev/user_guide/io.html#hdf5-pytables


store = pd.HDFStore("templates/store.h5")
print(store)

index = pd.date_range("1/1/2000", periods=8)
s = pd.Series(np.random.randn(5), index=["a", "b", "c", "d", "e"])
df = pd.DataFrame(np.random.randn(8, 3), index=index, columns=["A", "B", "C"])

# store.put('s',s) is an equivalent method
store["s"] = s
store["df"] = df

# store.get('df') is an equivalent method
print(store["df"])
print(store.df)

# store.remove('df') is an equivalent method
del store["df"]

store.close()
print(store.is_open)

# Read/Write API
df_tl = pd.DataFrame({"A": list(range(5)), "B": list(range(5))})
df_tl.to_hdf("templates/store_tl.h5", "table", append=True)
print(pd.read_hdf("templates/store_tl.h5", "table", where=["index>2"]))

df_with_missing = pd.DataFrame(
    {
        "col1": [0, np.nan, 2],
        "col2": [1, np.nan, np.nan],
    }
)
df_with_missing.to_hdf("templates/file.h5",
                       "df_with_missing", format="table", mode="w")
df_with_missing.to_hdf("file.h5", "df_with_missing",
                       format="table", mode="w", dropna=True)


"""The examples above show storing using put, which write the HDF5 to PyTables in a fixed array format, 
called the fixed format. These types of stores are not appendable once written (though you can simply remove them and rewrite).
 Nor are they queryable; they must be retrieved in their entirety. 
 They also do not support dataframes with non-unique column names. 
 The fixed format stores offer very fast writing and slightly faster reading than table stores. 
 This format is specified by default when using put or to_hdf or by format='fixed' or format='f'.

pd.DataFrame(np.random.randn(10, 2)).to_hdf("test_fixed.h5", "df")
pd.read_hdf("test_fixed.h5", "df", where="index>5")
"""

store2 = pd.HDFStore("templates/store.h5")
df1 = df[0:4]
df2 = df[4:]
# append data (creates a table automatically)
store2.append("df", df1)
store2.append("df", df2)
print(store2)
print(store2.select("df"))
print(store2.root.df._v_attrs.pandas_type)


# KEYS
store3 = pd.HDFStore("templates/store.h5")
store3.put("foo/bar/bah", df)
store3.append("food/orange", df)
store3.append("food/apple", df)
store3.keys()
store3.remove("food")

# WALK
for (path, subgroups, subkeys) in store3.walk():
    for subgroup in subgroups:
        print("GROUP: {}/{}".format(path, subgroup))
    for subkey in subkeys:
        key = "/".join([path, subkey])
        print("KEY: {}".format(key))
        print(store3.get(key))

# tips
print(store3.root.foo.bar.bah)
print(store3["foo/bar/bah"])


# Storing mixed types in a table
df_mixed = pd.DataFrame(
    {
        "A": np.random.randn(8),
        "B": np.random.randn(8),
        "C": np.array(np.random.randn(8), dtype="float32"),
        "string": "string",
        "int": 1,
        "bool": True,
        "datetime64": pd.Timestamp("20010102"),
    },
    index=list(range(8)),
)

df_mixed.loc[df_mixed.index[3:5], ["A", "B", "string", "datetime64"]] = np.nan
store3.append("df_mixed", df_mixed, min_itemsize={"values": 50})
df_mixed1 = store3.select("df_mixed")
print(df_mixed1)
print(df_mixed1.dtypes.value_counts())
print(store3.root.df_mixed.table)


# Storing Multindex DataFrames
index = pd.MultiIndex(
    levels=[["foo", "bar", "baz", "qux"], ["one", "two", "three"]],
    codes=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3], [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
    names=["foo", "bar"],
)
df_mi = pd.DataFrame(np.random.randn(
    10, 3), index=index, columns=["A", "B", "C"])
print(df_mi)
store3.append("df_mi", df_mi)
print(store3.select("df_mi"))
print(store3.select("df_mi", "foo=bar"))


# Querying
dfq = pd.DataFrame(
    np.random.randn(10, 4),
    columns=list("ABCD"),
    index=pd.date_range("20130101", periods=10),
)
store = pd.HDFStore("templates/store.h5")
store.append("dfq", dfq, format="table", data_columns=True)
print(store.select(
    "dfq", "index>pd.Timestamp('20130104') & columns=['A', 'B']"))
print(store.select("dfq", where="A>0 or C>0"))

# Query timedelta64(ns)

dftd = pd.DataFrame({"A": pd.Timestamp("20130101"), "B": [pd.Timestamp(
    "20130101") + timedelta(days=i, seconds=10)for i in range(10)], })
dftd["C"] = dftd["A"] - dftd["B"]
store.append("dftd", dftd, data_columns=True)
print(store.select("dftd", "C<'-3.5D'"))

# Query Multindex
print(df_mi.index.names)
print(store.select("df_mi", "foo=baz and bar=two"))
index = pd.MultiIndex(
    levels=[["foo", "bar", "baz", "qux"], ["one", "two", "three"]],
    codes=[[0, 0, 0, 1, 1, 2, 2, 3, 3, 3], [0, 1, 2, 0, 1, 1, 2, 0, 1, 2]],
)
df_mi_2 = pd.DataFrame(np.random.randn(
    10, 3), index=index, columns=["A", "B", "C"])
store.append("df_mi_2", df_mi_2)
store.select("df_mi_2", "level_0=foo and level_1=two")

# indexing
i = store.root.df.table.cols.index.index
print(i.optlevel, i.kind)
store.create_table_index("df", optlevel=9, kind="full")
i = store.root.df.table.cols.index.index
print(i.optlevel, i.kind)

df_1 = pd.DataFrame(np.random.randn(10, 2), columns=list("AB"))
df_2 = pd.DataFrame(np.random.randn(10, 2), columns=list("AB"))
st = pd.HDFStore("appends.h5", mode="w")
st.append("df", df_1, data_columns=["B"], index=False)
st.append("df", df_2, data_columns=["B"], index=False)
print(st.get_storer("df").table)
print(st.df)

# Query via data columns
df_dc = df.copy()
df_dc["string"] = "foo"
df_dc.loc[df_dc.index[4:6], "string"] = np.nan
df_dc.loc[df_dc.index[7:9], "string"] = "bar"
df_dc["string"] = "cool"
df_dc.loc[df_dc.index[1:3], ["B", "C"]] = 1.0
print(df_dc)

store.append("df_dc", df_dc, data_columns=["B", "C", "string", "string2"])
print(store.select("df_dc", where="B > 0"))
print(store.select("df_dc", "B > 0 & C > 0 & string == foo"))
print(df_dc[(df_dc.B > 0) & (df_dc.C > 0) &
      (df_dc.string == "foo")])  # in memory

# Interator
for df in store.select("df", chunksize=3):
    print(df)

dfeq = pd.DataFrame({"number": np.arange(1, 11)})
store.append("dfeq", dfeq, data_columns=["number"])


def chunck(l, n):
    return[l[i: i + n] for i in range(0, len(l), n)]


evens = [2, 4, 6, 8, 10]
coordinates = store.select_as_coordinates("dfeq", "number=evens")
for c in chunck(coordinates, 2):
    print(store.select("dfeq", where=c))


# Advanced queries
print(store.select_column("df_dc", "index"))
print(store.select_column("df_dc", "string"))

# select coordinates
df_coord = pd.DataFrame(np.random.randn(
    1000, 2), index=pd.date_range("20000101", periods=1000))
store.append("df_coord", df_coord)
c = store.select_as_coordinates("df_coord", "index > 20020101")
print(c)
print(store.select("df_coord", where=c))

# Selecting using a where mask
df_mask = pd.DataFrame(np.random.randn(
    1000, 2), index=pd.date_range("20000101", periods=1000))
store.append("df_mask", df_mask)
c = store.select_column("df_mask", "index")
where = c[pd.DatetimeIndex(c).month == 5].index
store.select("df_mask", where=where)

# Storer Object
print(store.get_storer("df_dc").nrows)

# Multiple Table queries
df_mt = pd.DataFrame(np.random.randn(8, 6), index=pd.date_range(
    "1/1/2000", periods=8), columns=["A", "B", "C", "D", "E", "F"],)
df_mt["foo"] = "bar"
df_mt.loc[df_mt.index[1], ("A", "B")] = np.nan

store.append_to_multiple(
    {"df1_mt": ["A", "B"], "df2_mt": None}, df_mt, selector="df1_mt")
print(store.select("df1_mt"))
print(store.select("df2_mt"))
print(store.select_as_multiple(
    ["df1_mt", "df2_mt"], where=["A>0", "B>0"], selector="df1_mt",))

# Delete from a table
"""After delete the menory still allocated, to return the memory, use ptrepack"""

# Compression
store_compressed = pd.HDFStore(
    "store_compressed.h5", complevel=9, complib="blosc:blosclz"
)

store.append("df", df, complib="zlib", complevel=5)  # or on the fly

# ptrepack
# ptrepack --chunkshape=auto --propindexes --complevel=9 --complib=blosc in.h5 out.h5

# DataTypes
# Categorical data
dfcat = pd.DataFrame({"A": pd.Series(list("aabbcdba")).astype(
    "category"), "B": np.random.randn(8)})
print(dfcat)
cstore = pd.HDFStore("cats.h5", mode='w')
cstore.append("dfcat", dfcat, format="table", data_columns=["A"])
result = cstore.select("dfcat", where="A in ['b','c']")
print(result)
print(result.dtypes)

# String Columns
