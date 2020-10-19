#!/usr/bin/env python
import sys
from os import listdir
from os.path import isfile, join
import sqlite3
import pandas as pd

strPath = sys.argv[1]
sID = sys.argv[2]

csv = [] 
for f in listdir(strPath):
  if (not isfile(join(strPath, f))) or ("csv" not in f):
    continue
  strCSV = join(strPath,f)
  tab = pd.read_csv(strCSV)
  tab.columns = ["gene","log2fc","pval","qval"]
  tags = f[:-4].split("_")
  tab["contrast"] = tags[0]
  tab["tags"] = ";".join(tags[1:])
  csv += [tab]

data = pd.concat(csv,ignore_index=True)
data = data.dropna()
D = data[["log2fc","pval","qval"]]
D.index = pd.MultiIndex.from_frame(data[["gene","contrast","tags"]])

conn = sqlite3.connect(join(strPath, sID+'.db'))
D.to_sql("DEG",conn,if_exists="replace")
conn.close()

print("sqlite3 db (%s) was created successfully with the first 5 rows:"%join(strPath, sID+'.db'))
conn = sqlite3.connect(join(strPath, sID+'.db'))
df = pd.read_sql_query("select * from DEG limit 5;", conn)
conn.close()
print(df)
