import os
import shlex
import subprocess
from typing import Optional
from multiprocessing import Pool
from .ref_prep import ref_prep
from .utils import get_files_path, get_cpu_core

def bed_to_bigwig(bed, prefix, bw_dir, chrom_sizes, normalized = True):
  bg = os.path.join(bw_dir, f"{prefix}.bedGraph")
  sort_bg = os.path.join(bw_dir, f"{prefix}.sorted.bedGraph")
  bw = os.path.join(bw_dir, f"{prefix}.bw")

  with open(bed, "r") as f:
    total_reads = sum(1 for line in f if line.strip())
  if total_reads == 0:
    print(f"Empty BED file: {bed}")
    return None
  
  scale_opt = ""
  if normalized:
    scale_factor = 1e6 / total_reads
    scale_opt = f"-scale {scale_factor}"
  
  subprocess.run(
    [
      "bash", "-c",
      f"bedtools genomecov -bg {scale_opt} -i {shlex.quote(bed)} -g {shlex.quote(chrom_sizes)} > {shlex.quote(bg)}"
    ],
     check=True
  )
  subprocess.run(
    [
      "bash", "-c",
      f"LC_ALL=C sort -k1,1 -k2,2n {shlex.quote(bg)} > {shlex.quote(sort_bg)}"
    ],
     check=True
  )
  subprocess.run(
    ["bedGraphToBigWig", sort_bg, chrom_sizes, bw],
    check=True
  )

  for p in [bg, sort_bg]:
    if os.path.exists(p):
      os.remove(p)
  return bw

def convert_bed_to_bigwig(
  bed_dir: str,
  ref_genome:str,
  out_dir:str,
  normalized: bool = True,
  threads: Optional[int] = None
):
  os.makedirs(out_dir, exist_ok=True)

  if threads is None:
    threads = min(8, get_cpu_core())

  ref_preparation = ref_prep()
  chrom_sizes = ref_preparation.get_chromsize(ref_genome)

  bed_list = get_files_path(bed_dir, ext = ".bed")
  with Pool(processes=threads) as pool:
    bigwig_list = []
    async_results = []
    for bed in bed_list:
      prefix = os.path.basename(bed).removesuffix('.bed')
      ar = pool.apply_async(bed_to_bigwig, args=(bed, prefix, out_dir, chrom_sizes, normalized))
      async_results.append((prefix, ar))
    
    for prefix, ar in async_results:
      try:
        bigwig_list.append(ar.get())
      except Exception as e:
        print(f"Failed: {prefix} -> {e}")
        pool.terminate()
        return None
      
  return bigwig_list

HELP = """
Convert BED files to bigWig tracks

Required:
  -i, --input     Directory containing BED files (.bed)
  -o, --output    Output directory where final BIGWIG files will be written
  -g, --genome    Reference genome name (hg38 or mm10)

Optional:
  -n, --normalized     Apply RPM normalization (default: true)
      --no-normalized  Disable RPM normalization
  -j, --threads        Number of parallel BED conversion jobs
                       (default: automatically detect all available CPU cores)

Output:
  OUT_DIR/{prefix}.bw

Description:
  Convert BED files to bigWig files by:
    1) name-sorting BAM
    2) converting to BEDPE
    3) extracting fragment coordinates
    4) sorting BED by genomic position

    
Example:
  multiEpiPrep track -i ./bed -o ./output -g hg38
  multiEpiPrep track -i ./bed -o ./output -g mm10 --no-normalized
"""

def register_parser(parser):
  parser.add_argument("-i", "--input", required=True)
  parser.add_argument("-o", "--output", required=True)
  parser.add_argument("-g", "--genome", required=True)
  parser.add_argument(
    "--no-normalized",
    dest="normalized",
    action="store_false",
    default=True,
    help="Disable RPM normalization"
  )
  parser.add_argument("-j", "--threads", type=int, default=None)
  parser.set_defaults(func=run)

def run(args):
  convert_bed_to_bigwig(
    bed_dir=args.input,
    out_dir=args.output,
    ref_genome=args.genome,
    normalized=args.normalized,
    threads=args.threads
  )
