import sqlite3
from pathlib import Path
import time
from datetime import datetime
import pandas as pd

ts = int("1284101485")
print(datetime.utcfromtimestamp(ts).strftime("%Y-%m-%d %H:%M:%S"))

DATA_DIR = Path(__file__).parent
db_GE = DATA_DIR / "GE.db"
v_file = DATA_DIR / "GE_Conn_load.csv"

df = pd.read_csv(v_file)

conn = sqlite3.connect(db_GE)
print("Opened database successfully")
df.to_sql("CONNECTIONS", conn, if_exists="append", index=False)

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

conn.commit()

print("Records created successfully")
# print("We have inserted", conn.rowcount, "records to the table.")
conn.close()
