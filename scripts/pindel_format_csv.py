import pandas as pd
import numpy as np
import os
import sys
import subprocess
import re

args = sys.argv
bed = args[1]
multi_anno_file = args[2]
outfile = args[3]

# Fix NTLEN format in multi_anno_file
with open(multi_anno_file, "r") as file:
    content = file.read()

content = re.sub(r'NTLEN=(\d+),(\d+)', r'NTLEN=\1;\2', content)

with open(multi_anno_file, "w") as file:
    file.write(content)

# Load BED file
df = pd.read_csv(bed, header=None, sep="\t")
print(df)

if os.stat(multi_anno_file).st_size != 0:
    df.columns = ["A", "B", "C", "D", "E"]
    df1 = pd.read_csv(multi_anno_file)
    df1.insert(loc=5, column='AF', value=0.000)
    df1.insert(loc=6, column='VAF', value=0.000)
    df1.insert(loc=7, column='AR', value=0.000)
    df1.insert(loc=8, column='Ref_Count', value=0.0)
    df1.insert(loc=9, column='Alt_Count', value=0.0)
    
    for i in range(len(df1)):
        df1.at[i, "Alt_Count"] = str(df1["Otherinfo1"][i]).split(",")[1]

    print(df1)
    print(df, "llllllll")
    print(df1["AF"])
    print(df1["VAF"])
    print(df1["Ref_Count"])
    print(df1["Alt_Count"])
    print((df1.dtypes), (df.dtypes))
    
    for i in range(len(df)):
        for j in range(len(df1)):
            if (df1["Start"][j] in range(df["B"][i], df["C"][i])) and df1["Chr"][j] == df["A"][i]:
                df1.at[j, "Ref_Count"] = df["E"][i] + 0.0
                df1.at[j, 'AF'] = df1['Alt_Count'][j] / df1['Ref_Count'][j]
                df1.at[j, "AR"] = df1['AF'][j] / (1 - df1['AF'][j])
                df1.at[j, 'VAF'] = df1['Alt_Count'][j] / df1['Ref_Count'][j] * 100
    
    print(df1)
    print(df1["AF"])
    print(df1["VAF"])
    print(df1["AR"])
    print(df1["Ref_Count"])
    print(df1["Alt_Count"])
    df1.to_csv(outfile, index=False)
else:
    subprocess.run(['touch', outfile], stdout=subprocess.DEVNULL)
