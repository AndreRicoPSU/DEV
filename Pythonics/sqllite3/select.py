import sqlite3
con = sqlite3.connect('exempl.db')
cur = con.cursor()
for row in cur.execute('SELECT * FROM Stocks'):
    print(row)