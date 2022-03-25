import sqlite3
from pathlib import Path
import time
from datetime import datetime

ts = int("1284101485")
print(datetime.utcfromtimestamp(ts).strftime("%Y-%m-%d %H:%M:%S"))


DATA_DIR = Path(__file__).parent
db_GE = DATA_DIR / "GE.db"

conn = sqlite3.connect(db_GE)
print("Opened database successfully")

"""
        ID
         ABBREVIATION       VARCHAR(20)     NOT NULL,         
         DATABASE           VARCHAR(60)    NOT NULL,
         CATEGORY           VARCHAR(60),
         DATASET            VARCHAR(60),
         SOURCE_PATH        VARCHAR(120),
         DESTINATION_PATH   VARCHAR(120),
         FILE_FORMAT        CHAR(5),
         FILE_NAME          VARCHAR(30),
         VERSION            CHAR(5),
         UPD_FILE           BLOB,
         SIZE               INT,
         LAST_UPDATE       INT  
"""

# Conectar um arquivo de carga
# Adicionar os registros no v_records
# drop table Conections

v_records = [
    (
        2,
        "HMDB",
        "Human Metabolome Database",
        "Protein/Gene Sequences",
        "All Metabolite Metabolizing Enzymes",
        "https://hmdb.ca/system/downloads/current/sequences/protein.fasta.zip",
        "default",
        "fasta",
        "protein_fasta",
        "v1.1",
        False,
        333333,
        1521462189,
    ),
    (
        1,
        "CTD",
        "Comparative Toxicogenomics Database",
        "Gene–pathway associations",
        "Gene–pathway associations",
        "http://ctdbase.org/reports/CTD_genes_pathways.csv.gz",
        "default",
        "csv",
        "CTD_genes_pathways.csv.gz",
        "v1.2",
        True,
        333333,
        1521462189,
    ),
]

conn.executemany(
    "INSERT INTO CONNECTIONS (ID,ABBREVIATION,DATABASE,CATEGORY,DATASET,SOURCE_PATH,DESTINATION_PATH,FILE_FORMAT,FILE_NAME,VERSION,UPD_FILE,SIZE,LAST_UPDATE) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?);",
    v_records,
)


conn.commit()

print("Records created successfully")
print("We have inserted", conn.rowcount, "records to the table.")
conn.close()
