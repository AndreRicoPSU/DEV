import sqlite3
from sqlite3 import Error
from pathlib import Path


DATA_DIR = Path(__file__).parent / "database"
db_GE = DATA_DIR / "GE.db"


conn = sqlite3.connect(db_GE)
print("Opened database successfully")

conn.execute("""DROP table CONNECTIONS""")

conn.execute(
    """
    CREATE TABLE CONNECTIONS
         (ID INT      NOT NULL,
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
         );"""
)

print("Table created successfully")

conn.close()
