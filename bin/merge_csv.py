#!/usr/bin/env python3

import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser(description="Take the combined output & coverage output as input files and write them to 2 sheets in 1 excel file")
    
parser.add_argument("-com", "--combined", required=True, help="Path to combined output")
parser.add_argument("-cov", "--coverage", required=True, help="Path to coverage output")
parser.add_argument("-o", "--output", required=True, help="Path to output excel file")

args = parser.parse_args()

combined_df = pd.read_csv(args.combined)
count_df = pd.read_csv(args.coverage, sep='\t', names=["Chr", "Start", "End","Region","Reads"])

new_sheetname = os.path.split(args.combined)[1]
new_sheetname = os.path.splitext(new_sheetname)[0]

with pd.ExcelWriter(args.output) as writer:
    combined_df.to_excel(writer, sheet_name=new_sheetname, index=False)
    count_df.to_excel(writer, sheet_name='coverage', index=False)

