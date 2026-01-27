#!/usr/bin/env python3

import sys
import os
import pandas as pd
from collections import defaultdict

def read_excel_data(filepath):
    xls = pd.ExcelFile(filepath)
    variant_sheet = [s for s in xls.sheet_names if s.lower() != 'coverage'][0]
    df_variants = xls.parse(variant_sheet)
    df_coverage = xls.parse('coverage')
    return df_variants, df_coverage

def extract_prefix(filename):
    base = os.path.basename(filename)
    return base.split('-')[0]

def extract_id(filename):
    base = os.path.basename(filename)
    return base.split('_')[0]

def ensure_columns_exist(df, expected_cols):
    for col in expected_cols:
        if col not in df.columns:
            df[col] = pd.NA
    return df

def main(args):
    grouped_files = defaultdict(list)

    for file in args:
        prefix = extract_prefix(file)
        grouped_files[prefix].append(file)

    for sample_id, files in grouped_files.items():
        if len(files) < 2:
            print(f"Skipping {sample_id}: only {len(files)} file(s) found.")
            continue

        # Ensure fileA is always chosen before fileB
        fileA = [f for f in files if "_A" in f or "-A" in f]
        fileB = [f for f in files if "_B" in f or "-B" in f]

        if fileA and fileB:
            file1, file2 = fileA[0], fileB[0]
        else:
            # fallback: just take first two sorted files
            file1, file2 = sorted(files)[:2]

        df1_variants, df1_coverage = read_excel_data(file1)
        df2_variants, df2_coverage = read_excel_data(file2)

        # Full list of comparison columns to ensure they exist
        comparison_cols = [
            'VAF', 'REF_COUNT', 'ALT_COUNT',
            'Variant Site', 'Gene.refGene', 'GeneDetail.refGene',
            'Variant Function', 'AAChange.refGene',
            'cosmic84', 'ExAC_ALL', 'ExAC_AFR', 'ExAC_AMR',
            'ExAC_EAS', 'ExAC_FIN', 'ExAC_NFE', 'ExAC_OTH', 'ExAC_SAS',
            'Mutation', 'Genomic', 'Protein', 'Nucleotide', 'PMID'
        ]

        df1_variants = ensure_columns_exist(df1_variants, comparison_cols)
        df2_variants = ensure_columns_exist(df2_variants, comparison_cols)

        output_file = f"{sample_id}.xlsx"
        with pd.ExcelWriter(output_file, engine='xlsxwriter') as writer:
            # Order: Variants_A → Coverage_A → Variants_B → Coverage_B
            df1_variants.to_excel(writer, sheet_name=f'Variants_{extract_id(file1)}', index=False)
            df1_coverage.to_excel(writer, sheet_name=f'Coverage_{extract_id(file1)}', index=False)
            df2_variants.to_excel(writer, sheet_name=f'Variants_{extract_id(file2)}', index=False)
            df2_coverage.to_excel(writer, sheet_name=f'Coverage_{extract_id(file2)}', index=False)

            # default: empty A1B1
            df_empty = pd.DataFrame()
            write_A1B1 = False

            if not df1_variants.empty and not df2_variants.empty:
                merge_keys = ['Chr', 'Start', 'End', 'Ref', 'Alt']
                try:
                    merged_df = pd.merge(
                        df1_variants, df2_variants,
                        on=merge_keys, how='inner',
                        suffixes=('_1', '_2')
                    )

                    if 'VAF_1' in merged_df.columns and 'VAF_2' in merged_df.columns:
                        extracted_diff = merged_df

                        selected_columns = merge_keys.copy()
                        for col in comparison_cols:
                            col_1, col_2 = f"{col}_1", f"{col}_2"
                            if col_1 in merged_df.columns and col_2 in merged_df.columns:
                                selected_columns.extend([col_1, col_2])

                        df_final = extracted_diff[selected_columns]
                        df_final.to_excel(writer, sheet_name='A1B1', index=False)
                        write_A1B1 = True
                    else:
                        print(f"⚠ {sample_id}: VAF columns missing, writing empty A1B1.")
                except KeyError as e:
                    print(f"⚠ {sample_id}: Missing merge keys {e}, writing empty A1B1.")

            if not write_A1B1:
                df_empty.to_excel(writer, sheet_name='A1B1', index=False)

        print(f"✔ Merged: {file1}, {file2} → {output_file}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 mergeA1B1.py <xlsx1> <xlsx2> ...")
        sys.exit(1)
    main(sys.argv[1:])

