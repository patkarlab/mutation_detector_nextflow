import pandas as pd
import os
import sys

args = sys.argv
combined_out = args[1]
counts_file = args[2]
out_file = args[3]

new_sheetname = os.path.split(combined_out)[1]
new_sheetname = os.path.splitext(new_sheetname)[0]

combined_df = pd.read_csv(combined_out)
count_df = pd.read_csv(counts_file, sep='\t', names=["Chr", "Start", "End","Region","Reads"])

#writer = pd.ExcelWriter(out_file)
with pd.ExcelWriter(out_file) as writer:  
    combined_df.to_excel(writer, sheet_name=new_sheetname, index=False)
    count_df.to_excel(writer, sheet_name='coverage', index=False)

