#!/usr/bin/env python
 
import pandas as pd
import os, re, sys, shutil

def detect_cols(df):
    df.columns = [col.lower().strip() for col in df.columns]
    target_col, barcode_col = None, None

    target_cols = ["target", "crf"]
    barcode_cols = ["bc", "barcode", "sequence", "index"]
    for col in df.columns:
        if target_col is None and col in target_cols:
            target_col = col
        if barcode_col is None and col in barcode_cols:
            barcode_col = col
        if target_col is not None and barcode_col is not None:
            break
    
    if target_col is None or barcode_col is None:
        print(f"Warning: Could not find required columns in barcode input file. Attempting to use first two columns.")
        pat = re.compile(r"^[ACGT]+$", re.IGNORECASE)
        for i in range(min(2, len(df.columns))):
            colname = df.columns[i]
            col_series = df.iloc[:, i].apply(lambda x: str(x).strip() if pd.notnull(x) else "")
            is_seq = col_series.apply(lambda x: bool(pat.fullmatch(x))).all()
            lengths = df.iloc[:, i].apply(lambda x: len(str(x).strip()) if pd.notnull(x) else pd.NA).dropna()
            mode_len = int(lengths.mode().iloc[0]) if not lengths.mode().empty else 0
            equal_len = (lengths == mode_len).all() if mode_len > 0 else False
            if is_seq and equal_len and barcode_col is None:
                barcode_col = colname
            else:
                if target_col is None:
                    target_col = colname
    
    if barcode_col is None:
        raise ValueError("Could not identify barcode column in input file.")
    return target_col, barcode_col

def create_barcode_fasta(input, output_dict):
    try:
        if input.endswith(('.xlsx', '.xls')):
            barcode_df = pd.read_excel(input)
        elif input.endswith(('.tsv')):
            barcode_df = pd.read_csv(input, sep='\t')
        elif input.endswith(('.csv')):
            barcode_df = pd.read_csv(input)
        else:
            raise ValueError("Unsupported file format. Please provide an Excel (.xlsx/.xls) or TSV (.tsv) file.")

        target_col, barcode_col = detect_cols(barcode_df)
        print(f"Detected target column: {target_col}, barcode column: {barcode_col}")
        
        # Clean target names
        barcode_df[target_col] = barcode_df[target_col].str.replace(r'[ /.]', '_', regex=True).str.replace('-', '', regex=False)
        
        # Creat fasta context
        fasta_lines = []
        for _, row in barcode_df.iterrows():
            fasta_lines.append(f">{row[target_col]}")
            fasta_lines.append(row[barcode_col])
        fasta_content = '\n'.join(fasta_lines)
        print(f"Generated FASTA with {len(fasta_lines)//2} barcodes")

        # Same content for both forward and reverse
        output_files = [
            "barcode_fwd.fasta",
            "barcode_rev.fasta"
        ]
        output_paths = [os.path.join(output_dict, filename) for filename in output_files]
        for filepath in output_paths:
            with open(filepath, 'w') as f:
                f.write(fasta_content)

        print(f"Forward barcode FASTA files: {output_paths[0]}")
        print(f"Reverse barcode FASTA files: {output_paths[1]}")
        print("Barcode FASTA creation completed successfully!")

    except Exception as e:
        print(f"Error: Barcode FASTA creation failed with {e}")
        sys.exit(1)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        prog="create_barcode_fasta",
        description="Create FASTA files from barcode input"
    )
    parser.add_argument("-i", "--input", required=True, help="Input barcode file path")
    parser.add_argument("-o", "--output_dict", required=True, help="Output directory path")

    args = parser.parse_args()

    try:
        create_barcode_fasta(args.input, args.output_dict)
    except Exception as e:
        sys.stderr.write(f"[barcode] ERROR: {e}\n")
        sys.exit(1)