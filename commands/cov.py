import os
import shlex
import subprocess
from typing import Optional
from multiprocessing import Pool
from .ref_prep import ref_prep
from .utils import get_files_path, get_cpu_core

def bed_to_bedgraph(bed, prefix, bg_dir, chrom_sizes):
  bg = os.path.join(bg_dir, f"{prefix}.bedGraph")
  subprocess.run(
    ["bash", "-c", f"bedtools genomecov -bg -i {shlex.quote(bed)} -g {shlex.quote(chrom_sizes)} > {shlex.quote(bg)}"],
    check=True
  )
  return bg
  
def convert_bed_to_bedgraph(
  bed_dir: str,
  ref_genome:str,
  out_dir:str,
  threads: Optional[int] = None
):
  os.makedirs(out_dir, exist_ok=True)

  if threads is None:
    threads = min(8, get_cpu_core())

  ref_preparation = ref_prep()
  chrom_sizes = ref_preparation.get_chromsize(ref_genome)

  bed_list = get_files_path(bed_dir, ext = ".bed")
  with Pool(processes=threads) as pool:
    async_results = []
    for bed in bed_list:
      prefix = os.path.basename(bed).removesuffix('.bed')
      ar = pool.apply_async(bed_to_bedgraph, args=(bed, prefix, out_dir, chrom_sizes))
      async_results.append((prefix, ar))
    
    bedgraph_list = []
    for prefix, ar in async_results:
      try:
        bedgraph_list.append(ar.get())
      except Exception as e:
        print(f"Failed: {prefix} -> {e}")
        pool.terminate()
        return None
      
  return bedgraph_list


HELP = """
Convert BED files to bedGraph coverage tracks

Required:
  -i, --input     Directory containing BED files (.bed)
  -o, --output    Output directory where final BEDGRAPH files will be written
  -g, --genome    Reference genome name (hg38 or mm10)

Optional:
  -j, --threads   Number of parallel BED conversion jobs
                  (default: automatically detect all available CPU cores)

Output:
  OUT_DIR/{prefix}.bedGraph

Example:
  multiEpiPrep cov -i ./bed -o ./output -g hg38
"""

def register_parser(parser):
  parser.add_argument("-i", "--input", required=True)
  parser.add_argument("-o", "--output", required=True)
  parser.add_argument("-g", "--genome", required=True)
  parser.add_argument("-j", "--threads", type=int, default=None)
  parser.set_defaults(func=run)

def run(args):
  convert_bed_to_bedgraph(
    bed_dir=args.input,
    out_dir=args.output,
    ref_genome=args.genome,
    threads=args.threads
  )