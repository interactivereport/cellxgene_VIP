#!/usr/bin/env python
#Currently VIP only support one method/one model of pre-computed DE, and sqlite3 (python) format is used to store the pre-computed results.

# 1. Organize/Rename the pre-computed DEGs csv file names as: ConditionA.vs.ConditionB_celltype1.csv, each of which contains four columns in order: ["gene","log2fc","pval","qval"]
# 2. Put all pre-coputed DEGs csv files into a folder, such as /path/to/preDEG/
# 3. Activate VIP env, such as: conda activate cellxgeneVIP
# 4. Run the attached python script with two parameters as: (python) ./DEG2sqlite3.py /path/to/preDEG/ AD (you can change the script to avoid do 1&2 steps above )
# 5. If successful, a AD.db will be created in the path of /path/to/preDEG/
# 6. Copy/rename the db file the same as h5ad name in the same folder
# 7. Refresh the browser (clear cache), pre-DEG button should be shown.

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
