import pandas as pd
import os, sys

args = sys.argv
sample = args[1]
filepath = args[2]
outfile = args[3]
cava_path = args[4]
coverview_path = args[5]
pindel_path = args[6]

csvfilenames=[filepath+sample+'.concat.ss.combined.csv',cava_path+sample+'.cava.csv',pindel_path,coverview_path]

writer = pd.ExcelWriter(outfile)
for csvfilename in csvfilenames:
    sheetname=os.path.split(csvfilename)[1]
    df = pd.read_csv(csvfilename)
    print('process file:', csvfilename, 'shape:', df.shape)
    df.to_excel(writer,sheet_name=os.path.splitext(sheetname)[0], index=False)
writer.save()

