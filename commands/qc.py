import os
import subprocess
from typing import List, Dict, Tuple, Optional
from .utils import get_files_path, get_files_prefix

def count_mapped_reads(bam: str) -> int:
  proc = subprocess.check_output(["samtools", "idxstats", bam], text=True)
  total = 0
  for line in proc.splitlines():
    parts = line.split("\t")
    if len(parts) >= 3:
      try:
        total += int(parts[2])
      except ValueError:
        pass
  return total

def write_readcount_tsv(path: str, rows: List[Tuple[str, str, int]], target_col = "pair", read_count_col = "read_count") -> None:
  with open(path, "w", encoding="utf-8") as f:
    f.write(f"{target_col}\t{read_count_col}\n")
    for _, prefix, count in rows:
      f.write(f"{prefix}\t{count}\n")

def qc_by_percentile(
  bam_dir: str,
  out_dir: str,
  percentile: float = 0.25,
  exclude: List[str] = ["unknown", "IgG_control"]
):
  if not (0.0 <= percentile <= 1.0):
    raise ValueError("percentile must be between 0 and 1")
  
  os.makedirs(out_dir, exist_ok=True)
  all_tsv = os.path.join(out_dir, "all_read_count.tsv")
  filtered_tsv = os.path.join(out_dir, "filtered_read_count.tsv")

  bam_list = get_files_path(bam_dir, ext = ".bam")
  prefix_list = get_files_prefix(bam_list,'.bam')
  has_dash = ["-" in prefix for prefix in prefix_list]

  if all(has_dash):
    mode = "all_dash"
    should_exclude = [any(p in exclude for p in prefix.split("-")) for prefix in prefix_list]
  elif not any(has_dash):
    mode = "no_dash"
    should_exclude = [prefix in exclude for prefix in prefix_list]
  else:
    mode = "mixed"
    should_exclude = [any(item in prefix for item in exclude) for prefix in prefix_list]

  prefix_to_bam = {
      prefix: bam
      for prefix, bam, _exclude in zip(prefix_list, bam_list, should_exclude) if not _exclude
  }

  try: 
    # Get all read counts
    rows = [
        (bam, prefix, count_mapped_reads(bam))
        for prefix, bam in prefix_to_bam.items()
      ]
    rows = sorted(rows, key=lambda x: x[2], reverse=True) 
    write_readcount_tsv(all_tsv, rows)
    
    # Filter read counts
    counts_sorted = sorted([c for _, _, c in rows])
    threshold = counts_sorted[int((len(counts_sorted)-1)*percentile)]
    filtered_rows = [(bam, prefix, c) for bam, prefix, c in rows if c >= threshold]
    write_readcount_tsv(filtered_tsv, filtered_rows)

    # Get passed .bam files
    filtered_bam_list = [bam for bam, _, _ in filtered_rows]
  
  except Exception as e:
    print(e)
    return None

  return filtered_bam_list

HELP = """
Perform read-count-based quality control on BAM files.

Required:
  -i, --input       Directory containing input BAM files
  -o, --output      Output directory for QC results

Optional:
  -p, --percentile  Percentile threshold for filtering (default: 0.25)
  -e, --exclude     Comma-separated list of CRF names to exclude (default: unknown,IgG_control)
                    Exclusion is performed by matching substrings in BAM filenames. If -e/--exclude is provided without a value, the exclude list is empty.

Description:
  Computes the total number of mapped reads for each BAM file belonging to one specific sample by running `samtools idxstats`. Based on the resulting read counts, BAM files are filtered using a rank-based percentile threshold, discarding low-coverage files while retaining those above the cutoff.

Example:
  multiEpiPrep qc -i BAM_DIR -o OUT_DIR
"""

def register_parser(parser):
  parser.add_argument("-i", "--input", required=True)
  parser.add_argument("-o", "--output", required=True)
  parser.add_argument("-p", "--percentile", type=float, default=0.25)
  parser.add_argument(
    "-e", "--exclude",
    type=lambda s: s.split(",") if s else [],
    default="unknown,IgG_control",
    metavar="NAMES"
  )
  parser.set_defaults(func=run)

def run(args):
  qc_by_percentile(
    bam_dir=args.input,
    out_dir=args.output,
    percentile=args.percentile,
    exclude=args.exclude
  )
