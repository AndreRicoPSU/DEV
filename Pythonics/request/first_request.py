# import pyFTPclient
from genericpath import exists
import requests
from pathlib import Path
import os
import sqlite3


"""
1. atualizar a tabela com a versao, tamanho e a data do download
2. Criar uma tabela com os logs
3. Sera que da para implementar o multiprocesso?
4. Disparar um email apos o termino com sucesso ou erros (LOG)
5. A versao sera a data do aquivo (ver se consigo no header)

"""


# adicionar um codigo para nao criar uma nova base se nao achar. pode ser via IF
v_path_base = Path(__file__).parent / "database"
v_db_GE = v_path_base / "GE.db"
conn = sqlite3.connect(v_db_GE)
print("Opened database successfully")

v_path_file = Path(__file__).parent / "Downloads"

"""
         0. ID
         1. ABBREVIATION       VARCHAR(20)     NOT NULL,         
         2. DATABASE           VARCHAR(60)    NOT NULL,
         3. CATEGORY           VARCHAR(60),
         4. DATASET            VARCHAR(60),
         5. SOURCE_PATH        VARCHAR(120),
         6. DESTINATION_PATH   VARCHAR(120),
         7. FILE_FORMAT        CHAR(5),
         8. FILE_NAME          VARCHAR(30),
         9. VERSION            CHAR(5),
         10. UPD_FILE           BLOB,
         11. SIZE               INT,
         12. LAST_UPDATE       INT  
"""

cursor = conn.execute("""SELECT * FROM CONNECTIONS WHERE UPD_FILE = 1""")


for row in cursor:
    print("Download Process:")
    print("ID = ", row[0])
    print("DB = ", row[1])
    print("DATASET = ", row[4])
    print("FILE = ", row[8], "\n")

    v_dir = v_path_file / row[1]

    if not os.path.isdir(v_dir):
        os.makedirs(v_dir)
        v_file = v_dir / row[8]
        v_size_old = 0
        print("created folder : ", v_dir, "\n")
    elif os.path.exists(v_dir / row[8]):
        v_file = v_dir / row[8]
        v_size_old = str(os.stat(v_file).st_size)
    else:
        v_file = v_dir / row[8]
        v_size_old = 0

    file_url = row[5]
    v_size_new = str(requests.get(file_url, stream=True).headers["Content-length"])

    if v_size_new == v_size_old:
        print("File is the same version \n")
    else:
        r = requests.get(file_url, stream=True)

        with open(v_file, "wb") as download:
            for chunk in r.iter_content(chunk_size=1024):

                # writing one chunk at a time to pdf file
                if chunk:
                    download.write(chunk)


print("Operation done successfully")
conn.close()
