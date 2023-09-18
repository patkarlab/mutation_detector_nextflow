#! /usr/bin/env python3
import sys
import pandas as pd
import os
import xlsxwriter

coverageA1 = sys.argv[1]
coverageB1 = sys.argv[2]
input_file_A = sys.argv[3]
input_file_B = sys.argv[4]
outfile = sys.argv[5]


df1 = pd.read_csv(input_file_A)
df2 = pd.read_csv(input_file_B)
common_columns = ['Chr', 'Start', 'End', 'Ref', 'Alt']
#merged_df = pd.merge(df1, df2, on=common_columns, how='inner', suffixes=('_A1', '_B1'))
#extracted_no_duplicates = merged_df[merged_df['VAF_A1'] != merged_df['VAF_B1']]
#header = merged_df.columns.tolist()

if all(column in df1.columns and column in df2.columns for column in common_columns):
		merged_df = pd.merge(df1, df2, on=common_columns, how='inner', suffixes=('_A1', '_B1'))
		extracted_no_duplicates = merged_df[merged_df['VAF_A1'] != merged_df['VAF_B1']]
		selected_columns =  ['Chr', 'Start', 'End', 'Ref', 'Alt', 'VAF_A1', 'VAF_B1',
		'REF_COUNT_A1', 'REF_COUNT_B1','ALT_COUNT_A1', 'ALT_COUNT_B1',
		'Variant Site_A1', 'Variant Site_B1', 'Gene.refGene_A1', 'Gene.refGene_B1',
		'GeneDetail.refGene_A1', 'GeneDetail.refGene_B1','Variant Function_A1', 'Variant Function_B1',
		'AAChange.refGene_A1','AAChange.refGene_B1',  'cosmic84_A1', 'cosmic84_B1',
		'ExAC_ALL_A1', 'ExAC_ALL_B1', 'ExAC_AFR_A1','ExAC_AFR_B1', 'ExAC_AMR_A1', 'ExAC_AMR_B1',
		'ExAC_EAS_A1', 'ExAC_EAS_B1', 'ExAC_FIN_A1', 'ExAC_FIN_B1', 'ExAC_NFE_A1', 'ExAC_NFE_B1',
		'ExAC_OTH_A1', 'ExAC_OTH_B1', 'ExAC_SAS_A1', 'ExAC_SAS_B1','Mutation_A1', 'Mutation_B1',
		'Genomic_A1', 'Genomic_B1', 'Protein_A1',  'Protein_B1', 'Nucleotide_A1', 'Nucleotide_B1', 'PMID_A1', 'PMID_B1']
		extracted_data = merged_df[selected_columns]
		extracted_data.to_csv('A1B1common.csv', index=False)
		df3 = pd.read_csv("A1B1common.csv", sep=",")

else:
		merged_df = pd.merge(df1, df2, how='inner', suffixes=('_A1', '_B1'))
		selected_columns = df1.columns.tolist() + df2.columns.tolist()
		extracted_data1 = merged_df[selected_columns]
		df3 = pd.read_csv("A1B1common.csv", sep=",")


csv_table1=pd.read_table(coverageA1,sep='\t',names=["Start", "End","Sequence","Reads"])
csv_table2=pd.read_table(coverageB1,sep='\t',names=["Start", "End","Sequence","Reads"])


excel_writer = pd.ExcelWriter(outfile, engine='xlsxwriter')

csv_table1_sheet_name = os.path.splitext(os.path.basename(coverageA1))[0]
csv_table2_sheet_name = os.path.splitext(os.path.basename(coverageB1))[0]
df1_sheet_name = os.path.splitext(os.path.basename(input_file_A))[0]
df2_sheet_name = os.path.splitext(os.path.basename(input_file_B))[0]
df3_sheet_name = 'A1B1Common'

csv_table1.to_excel(excel_writer, sheet_name=csv_table1_sheet_name, index=False)
df1.to_excel(excel_writer, sheet_name=df1_sheet_name, index=False)
csv_table2.to_excel(excel_writer, sheet_name=csv_table2_sheet_name, index=False)
df2.to_excel(excel_writer, sheet_name=df2_sheet_name, index=False)
df3.to_excel(excel_writer, sheet_name='A1B1', index=False)

#csv_table1.to_excel(excel_writer, sheet_name='coverage-A1', index=False)
#df1.to_excel(excel_writer, sheet_name='A1', index=False)
#csv_table2.to_excel(excel_writer, sheet_name='coverage-B1', index=False)
#df2.to_excel(excel_writer, sheet_name='B1', index=False)
#df3.to_excel(excel_writer, sheet_name='A1B1Common', index=False)

