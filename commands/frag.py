import os
import shlex
import subprocess
from typing import Optional
from multiprocessing import Pool
from .utils import get_files_path, get_cpu_core

def is_paired_bam(bam):
  result = subprocess.run(
    ["samtools", "view", "-c", "-f", "1", bam],
    capture_output=True,
    text=True,
    check=True
  )
  return int(result.stdout.strip()) > 0

def bam_to_bed(bam, prefix, bed_dir):
  sort_bam = os.path.join(bed_dir, f"{prefix}.sort.bam")
  bedpe = os.path.join(bed_dir, f"{prefix}.bedpe")
  tmp_bed = os.path.join(bed_dir, f"{prefix}.tmp.bed")
  bed = os.path.join(bed_dir, f"{prefix}.bed")

  subprocess.run(
    ["samtools", "sort", "-n", bam, "-o", sort_bam], 
    check=True
  )
  paired = is_paired_bam(bam)

  if paired: 
    # paired-end: name sort required
    subprocess.run(
      ["bash", "-c",
      f"bedtools bamtobed -i {shlex.quote(sort_bam)} -bedpe > {shlex.quote(bedpe)}"],
      check=True
    )
    subprocess.run(
      [
        "bash", "-c",
        f"""awk 'BEGIN{{OFS="\\t"}} NF>=6 && $1==$4 {{print $1,$2,$6}}' {shlex.quote(bedpe)} > {shlex.quote(tmp_bed)}"""
      ],
      check=True
    )
  else:
    # single-end: coordinate sort not required
    subprocess.run(
      [
        "bash", "-c",
        f"""bedtools bamtobed -i {shlex.quote(sort_bam)} | awk 'BEGIN{{OFS="\\t"}} {{print $1,$2,$3}}' > {shlex.quote(tmp_bed)}"""
      ],
      check=True
    )
  
  subprocess.run(
    ["sort", "-k1,1", "-k2,2n", tmp_bed, "-o", bed],
    check=True
  )

  # cleanup
  for p in [sort_bam, bedpe, tmp_bed]:
    if os.path.exists(p):
      os.remove(p)

  return bed


def convert_bam_to_bed(
  bam_dir: str,
  out_dir: str,
  threads: Optional[int] = None
):
  os.makedirs(out_dir, exist_ok=True)

  if threads is None:
    threads = min(8, get_cpu_core())

  bam_list = get_files_path(bam_dir, ext=".bam")
  with Pool(processes=threads) as pool:
    async_results = []
    for bam in bam_list:
      prefix = os.path.basename(bam).removesuffix('.bam')
      ar = pool.apply_async(bam_to_bed, args=(bam, prefix, out_dir))
      async_results.append((prefix, ar))

    bed_list = []
    for prefix, ar in async_results:
      try:
        bed_list.append(ar.get())
      except Exception as e:
        print(f"Failed: {prefix} -> {e}")
        pool.terminate()
        return None

  return bed_list


HELP = """
Convert BAM files into sorted BED

Required:
  -i, --input     Directory containing BAM files (.bam)
  -o, --output    Output directory where final BED files will be written

Optional:
  -j, --threads   Number of parallel BAM conversion jobs
                  (default: automatically detect all available CPU cores)

Output:
  OUT_DIR/{prefix}.bed

Description:
  Convert BAM files into sorted BED files by:
    1) name-sorting BAM
    2) converting to BEDPE
    3) extracting fragment coordinates
    4) sorting BED by genomic position

Example:
  multiEpiPrep frag -i ./bed -o ./output
"""

def register_parser(parser):
  parser.add_argument("-i", "--input", required=True)
  parser.add_argument("-o", "--output", required=True)
  parser.add_argument("-j", "--threads", type=int, default=None)
  parser.set_defaults(func=run)

def run(args):
  convert_bam_to_bed(
    bam_dir=args.input,
    out_dir=args.output,
    threads=args.threads
  )

