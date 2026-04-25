import os
import re
import sys
import subprocess
import gzip
import shutil
import tempfile
import pandas as pd
from typing import Optional, Tuple, List
from .utils import get_files_path, get_cpu_core, drop_empty_files, get_fastq_prefix
from .demux_refined import run_demux

def detect_cols(df:pd.DataFrame) -> Tuple[str, str]:
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

def create_barcode_fasta(input:str, out_dir:str):
  os.makedirs(out_dir, exist_ok=True)

  if input.endswith(('.xlsx', '.xls')):
    barcode_df = pd.read_excel(input)
  elif input.endswith(('.tsv')):
    barcode_df = pd.read_csv(input, sep='\t')
  elif input.endswith(('.csv')):
    barcode_df = pd.read_csv(input)
  else:
    raise ValueError("Unsupported file format. Please provide an Excel (.xlsx/.xls) or TSV (.tsv) file.")
  
  try:
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
    out_path = os.path.join(out_dir, "barcode.fasta")
    with open(out_path, "w") as f:
      f.write(fasta_content)
    print("Barcode FASTQ creation completed successfully!")
    return out_path
  except Exception as e:
    print(f"Error: Barcode FASTQ creation failed with {e}")
    return None


def merge_symmetric_fastq(demux_files: List[str]) -> List[str]:
  unique_combs_dict = get_fastq_prefix(demux_files)
  unique_combs = set(unique_combs_dict.keys())

  print(f"\nGenerated {len(unique_combs)} barcode pairs")
  print("\nStarting file merging process...")
  buf=8 * 1024 * 1024
  merged_fastq = []
  try:
    for comb in list(unique_combs):
      target1, target2 = comb.split("-")
      rev_comb = f"{target2}-{target1}" # reverse orientation counterpart of comb
      
      if target1 > target2:
        has_rev_comb = rev_comb in unique_combs
        r1 = unique_combs_dict[comb]['fwd']
        r2 = unique_combs_dict[comb]['rev']
          
        if has_rev_comb:
          rev_r1 = unique_combs_dict[rev_comb]["fwd"]
          rev_r2 = unique_combs_dict[rev_comb]["rev"]

          tmp1 = tempfile.NamedTemporaryFile(
            delete=False, dir=os.path.dirname(rev_r1), prefix=f"{rev_comb}_R1.tmp."
          )
          tmp2 = tempfile.NamedTemporaryFile(
            delete=False, dir=os.path.dirname(rev_r2), prefix=f"{rev_comb}_R2.tmp."
          )
          tmp1_path, tmp2_path = tmp1.name, tmp2.name
          tmp1.close(); tmp2.close()
          
          with gzip.open(tmp1_path, "wb") as w:
            for src in (rev_r1, r1):
              with gzip.open(src, "rb") as f:
                shutil.copyfileobj(f, w, length=buf)
          with gzip.open(tmp2_path, "wb") as w:
            for src in (rev_r2, r2):
              with gzip.open(src, "rb") as f:
                shutil.copyfileobj(f, w, length=buf)
          os.replace(tmp1_path, rev_r1)
          os.replace(tmp2_path, rev_r2)
          os.remove(r1)
          os.remove(r2)
            
        else:
          new_r1 = os.path.join(os.path.dirname(r1), f"{rev_comb}_R1.fastq.gz")
          new_r2 = os.path.join(os.path.dirname(r2), f"{rev_comb}_R2.fastq.gz")
          os.replace(r1, new_r1)
          os.replace(r2, new_r2)
          merged_fastq.append(new_r1)
          merged_fastq.append(new_r2)
      else:
        merged_fastq.append(unique_combs_dict[comb]['fwd'])
        merged_fastq.append(unique_combs_dict[comb]['rev'])

    print("Symmetric FASTQ merge completed successfully!")
    return merged_fastq
  except Exception as e:
    print(f"Error: symmetric FASTQ merge failed with: {e}")
    return None
  
def demultiplex(
  fastq_r1: str, 
  fastq_r2: str, 
  out_dir: str, 
  barcode: str, 
  error_rate: float = 0,
  threads: Optional[int] = None
):
  os.makedirs(out_dir, exist_ok=True)

  if threads is None:
    threads = get_cpu_core()
  
  if os.path.getsize(fastq_r1) > 1 * 1024 ** 3:
    try:
      run_demux(
        r1_path=fastq_r1,
        r2_path=fastq_r2,
        r1_barcode_file=barcode,
        output_dir=out_dir,
        error_rate=error_rate,
        streaming=True,
        gzip_output=True,
        workers=threads,
      )
      print("Fast demultiplexer execution successful")

      demux_fastq = drop_empty_files(get_files_path(out_dir, ext=".fastq.gz"))
      demux_fastq = [f for f in demux_fastq
                     if "unknown-unknown" not in os.path.basename(f).lower()]
      merged_fastq = merge_symmetric_fastq(demux_fastq)
      if merged_fastq is None:
        print("ERROR: Merging symmetric FASTQ files failed")
      else:
        print("Demultiplexing completed successfully!")
    except Exception as e:
      print(f"ERROR: Demultiplexing failed with {e}")

  else:
    barcode_fasta = None
    if os.path.basename(barcode).endswith('.fasta'):
      barcode_fasta = os.path.abspath(barcode)
    elif os.path.basename(barcode).endswith((".csv", ".tsv", ".xlsx", ".xls")):
      barcode_fasta = create_barcode_fasta(input=barcode, out_dir=os.path.dirname(barcode))
    if barcode_fasta is None or not os.path.exists(barcode_fasta):
      print("ERROR: Failed to create barcode FASTA file. Please check the input file and try again.")
      return
    
    try:
      cmd = [
        "cutadapt",
        "-e", str(error_rate),
        "-j", str(threads),
        "--no-indels",
        "--action", "trim",
        "-g", f"^file:{barcode_fasta}",
        "-G", f"^file:{barcode_fasta}",
        "-o", os.path.join(out_dir, "{name1}-{name2}_R1.fastq.gz"),
        "-p", os.path.join(out_dir, "{name1}-{name2}_R2.fastq.gz"),
        fastq_r1,
        fastq_r2
      ]
      subprocess.run(cmd, check=True)
      print("Cutadapt execution successful")

      demux_fastq = drop_empty_files(get_files_path(out_dir, ext = ".fastq.gz"))
      merged_fastq = merge_symmetric_fastq(demux_fastq)
      if merged_fastq is None:
        print("ERROR: Merging symmetric FASTQ files failed")
      else:
        print("Demultiplexing completed successfully!")
    except Exception as e:
      print(f"ERROR: Demultiplexing failed with {e}")

HELP = """
Demultiplex paired-end FASTQ files by barcode combinations using cutadapt

Required:
  -1, --r1       Path to the input FASTQ R1 file (.fastq.gz)
  -2, --r2       Path to the input FASTQ R2 file (.fastq.gz)
  -o, --out      Directory to save the demultiplexed FASTQ files
  -b, --barcode   Path to the barcode input file (Excel .xlsx/.xls or TSV .tsv or CSV .csv) or pre-generated FASTA file (.fasta)

Optional:
  -e, --error-rate  Max barcode mismatch rate for cutadapt (default: 0)
  -j, --threads     Number of threads to use (default: automatically detect all available CPU cores)

Output:
  OUT_DIR/{name1}-{name2}_R1.fastq.gz
  OUT_DIR/{name1}-{name2}_R2.fastq.gz

Description:
  Demultiplex paired-end FASTQ files by barcode combinations using cutadapt.
  - Uses linked adapter files:
      -g ^file:<FWD_BARCODES.fa>
      -G ^file:<REV_BARCODES.fa>
  - No indels allowed: --no-indels
  - No trimming performed: --action none
  - By default removes empty output FASTQ files (0 lines after gzip decompression)

Example:
  multiEpiPrep demux \\
    -1 ./raw_R1.fastq.gz -2 ./raw_R2.fastq.gz \\
    -b ./barcodes.csv \\
    -o ./demux_out \\
    -e 2 -j 16
"""
  
def register_parser(parser):
    parser.add_argument("-1", "--r1", required=True)
    parser.add_argument("-2", "--r2", required=True)
    parser.add_argument("-o", "--out", required=True)
    parser.add_argument("-b", "--barcode", required=True)
    parser.add_argument("-e", "--error-rate", type=float, default=0.0)
    parser.add_argument("-j", "--threads", type=int, default=None)
    parser.set_defaults(func=run)
  
def run(args):
  demultiplex(
    fastq_r1=args.r1,
    fastq_r2=args.r2,
    out_dir=args.out,
    barcode=args.barcode,
    error_rate=args.error_rate,
    threads=args.threads,
  )