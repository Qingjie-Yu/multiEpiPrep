import os
import subprocess
import gzip
import shutil
import tempfile
from typing import Optional, Dict, List
from .utils import get_files_path, get_cpu_core, drop_empty_files, get_fastq_prefix

def merge_symmetric_fastq(demux_files: List[str]) -> List[str]:
  unique_combs_dict = get_fastq_prefix(demux_files, ext=".fastq.gz")
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
  barcode_fwd_fasta: str, 
  barcode_rev_fasta: str, 
  error_rate: float = 0,
  threads: Optional[int] = None
):
  os.makedirs(out_dir, exist_ok=True)

  if threads is None:
    threads = get_cpu_core()
  
  try:
    cmd = [
      "cutadapt",
      "-e", str(error_rate),
      "-j", str(threads),
      "--no-indels",
      "--action", "trim",
      "-g", f"^file:{barcode_fwd_fasta}",
      "-G", f"^file:{barcode_rev_fasta}",
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
  
  return
  

HELP = """
Demultiplex paired-end FASTQ files by barcode combinations using cutadapt

Required:
  -1, --r1       Path to the input FASTQ R1 file (.fastq.gz)
  -2, --r2       Path to the input FASTQ R2 file (.fastq.gz)
  -o, --out      Directory to save the demultiplexed FASTQ files
  -f, --fwd      Path to FASTA file containing forward barcode sequences
  -r, --rev      Path to FASTA file containing reverse barcode sequences

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
    -f ./barcodes_fwd.fa -r ./barcodes_rev.fa \\
    -o ./demux_out \\
    -e 2 -j 16
"""
  
def register_parser(parser):
    parser.add_argument("-1", "--r1", required=True)
    parser.add_argument("-2", "--r2", required=True)
    parser.add_argument("-o", "--out", required=True)
    parser.add_argument("-f", "--fwd", required=True)
    parser.add_argument("-r", "--rev", required=True)
    parser.add_argument("-e", "--error-rate", type=float, default=0.0)
    parser.add_argument("-j", "--threads", type=int, default=None)
    parser.set_defaults(func=run)
  
def run(args):
  demultiplex(
    fastq_r1=args.r1,
    fastq_r2=args.r2,
    out_dir=args.out,
    barcode_fwd_fasta=args.fwd,
    barcode_rev_fasta=args.rev,
    error_rate=args.error_rate,
    threads=args.threads,
  )