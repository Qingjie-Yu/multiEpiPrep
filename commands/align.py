import os
import pysam
import shlex
import subprocess
from typing import Optional, Tuple, List
from .ref_prep import ref_prep
from .utils import get_files_path, get_fastq_prefix, get_picard_jar_path, get_cpu_core, get_memory

def generate_clean_bam(
  fastq_r1: str,
  fastq_r2: str,
  group: str,
  bam_out: str,
  index_prefix: str,
  min_length: int,
  max_length: int,
  threads: int = 1,
):
  cmd_str = (
    "set -euo pipefail; "
    f"bowtie2 "
    f"-x {shlex.quote(index_prefix)} "
    f"-I {min_length} "
    f"-X {max_length} "
    f"-p {threads} "
    f"-q "
    f"--local "
    f"--very-sensitive-local "
    f"--no-unal "
    f"--no-mixed "
    f"--no-discordant "
    f"--rg-id {group} "
    f"--rg SM:{group} "
    f"-1 {shlex.quote(fastq_r1)} "
    f"-2 {shlex.quote(fastq_r2)} "
    f"| samtools view -@ {threads} -h -q 10 - "
    r"""| awk 'BEGIN{OFS="\t"} /^@/ {print; next} ($3!="chrM" && $3 !~ /random/ && $3 !~ /^chrUn/) {print}' """
    f"| samtools view -@ {threads} -b - "
    f"> {shlex.quote(bam_out)}"
  )
  subprocess.run(["bash", "-c", cmd_str], check=True)
  return bam_out


def _extract_adapter_from_read(n):
  parts = n.split("::")
  if len(parts) < 2:
    raise ValueError(
        f"No adapter information found in read name '{n}'. "
        f"Expected format '...::[TAG1]-[TAG2]::[NAME1]=[TAG1]::[NAME2]=[TAG2]...'."
    )
  tags = {}
  for part in parts[2:]:
    if "=" in part:
      key, value = part.split("=", 1)
      tags[key] = value
  return tags


def extract_adapter_as_tags(
  bam: str,
  bam_out: str,
  cb_tag: Optional[str] = "CB",
  umi_tag: Optional[str] = "UB",
):
  assert (cb_tag is not None) or (umi_tag is not None)
  n_total, n_written, n_skipped = 0, 0, 0
  with pysam.AlignmentFile(bam, mode="rb") as fi, \
      pysam.AlignmentFile(bam_out, mode="wb", header=fi.header) as fo:
    for read in fi.fetch(until_eof=True):
      n_total += 1
      try:
        tags = _extract_adapter_from_read(read.query_name)
        if cb_tag is not None and tags.get("CB"):
          read.set_tag(cb_tag, tags["CB"], value_type="Z", replace=True)
        if umi_tag is not None and tags.get("UMI"):
          read.set_tag(umi_tag, tags["UMI"], value_type="Z", replace=True)
        n_written += 1
        fo.write(read)
      except ValueError:
        n_skipped += 1
      if n_total % 100000 == 0:
        print("Reads processed:", n_total)
  print("Read counts:", dict(input=n_total, written=n_written, skipped=n_skipped))
  os.remove(bam)
  return bam_out


def _dedup_bam_by_cb(bam_in: str, bam_out: str, cb_tag: str = "CB"):
  def _extract_primary_pair(
    bucket: List[pysam.AlignedSegment],
  ) -> Tuple[Optional[pysam.AlignedSegment], Optional[pysam.AlignedSegment]]:
    r1 = None
    r2 = None
    for read in bucket:
      if read.is_secondary or read.is_supplementary:
        continue
      if read.is_read1 and r1 is None:
        r1 = read
      elif read.is_read2 and r2 is None:
        r2 = read
    return r1, r2

  n_total, n_written, n_skipped, n_duplicated = 0, 0, 0, 0
  with pysam.AlignmentFile(bam_in, "rb") as fi, \
      pysam.AlignmentFile(bam_out, "wb", header=fi.header) as fo:
    seen = set()
    current_qname = None
    bucket: List[pysam.AlignedSegment] = []

    def _process_bucket(bucket: List[pysam.AlignedSegment]):
      nonlocal n_total, n_written, n_skipped, n_duplicated
      if not bucket:
        return
      n_total += 1
      r1, r2 = _extract_primary_pair(bucket)
      if r1 is None or r2 is None:
        n_skipped += 1
        return
      try:
        cb = r1.get_tag(cb_tag)
      except KeyError:
        n_skipped += 1
        return
      if cb is None:
        n_skipped += 1
        return
      _start = min(r1.reference_start, r2.reference_start)
      _end = max(r1.reference_end, r2.reference_end)
      key = (cb, r1.reference_id, _start, _end)
      if key in seen:
        n_duplicated += 1
        return
      seen.add(key)
      n_written += 1
      fo.write(r1)
      fo.write(r2)

    for read in fi.fetch(until_eof=True):
      qname = read.query_name
      if current_qname is None:
        current_qname = qname
        bucket = [read]
      elif qname == current_qname:
        bucket.append(read)
      else:
        _process_bucket(bucket)
        current_qname = qname
        bucket = [read]
    _process_bucket(bucket)  # final bucket

  print("Read counts:", dict(input=n_total, written=n_written, skipped=n_skipped, duplicated=n_duplicated))


def deduplicate_bam(
  bam: str,
  bam_out: str,
  cb_tag: Optional[str] = "CB",
  umi_tag: Optional[str] = "UB",
  threads: int = 1,
  java_mem: Optional[str] = None,
  picard_jar: Optional[str] = None,
):
  if umi_tag is not None:
    umi_dedup_log = f"{bam_out}.umi_dedup.log"
    if cb_tag is not None:
      # (cb-aware, umi-aware) dedup based on chrom coordination, cell barcode and UMI
      subprocess.run([
        "umi_tools", "dedup",
        "--paired",
        "-I", bam,
        "-S", bam_out,
        f"--cell-tag={cb_tag}",
        f"--umi-tag={umi_tag}",
        "--log", umi_dedup_log,
      ], check=True)
    else:
        # (umi-aware) dedup based on chrom coordination and UMI
      subprocess.run([
        "umi_tools", "dedup",
        "--paired",
        "-I", bam,
        "-S", bam_out,
        f"--umi-tag={umi_tag}",
        "--log", umi_dedup_log,
      ], check=True)

  else:
      if cb_tag is not None:
        # (cb-aware) dedup based on chrom coordination and cell barcode
        name_sort_bam = f"{bam_out}.namesort.tmp.bam"
        subprocess.run(["samtools", "sort", "-n", "-@", str(threads), "-o", name_sort_bam, bam], check=True)
        _dedup_bam_by_cb(name_sort_bam, bam_out, cb_tag)
        os.remove(name_sort_bam)
      else:
        # dedup based on only chrom coordination (Picard)
        assert java_mem is not None
        assert picard_jar is not None
        sort_bam = f"{bam_out}.sort.tmp.bam"
        subprocess.run(["samtools", "sort", "-@", str(threads), "-o", sort_bam, bam], check=True)
        metrics = f"{bam_out}.picard_metrics.log"
        subprocess.run([
          "java", f"-Xmx{java_mem}",
          "-jar", picard_jar,
          "MarkDuplicates",
          "--INPUT", sort_bam,
          "--OUTPUT", bam_out,
          "--METRICS_FILE", metrics,
          "--REMOVE_DUPLICATES", "true",
        ], check=True)
        os.remove(sort_bam)

  os.remove(bam)
  return bam_out


def sort_and_index_bam(
    bam: str,
    bam_out: str,
    threads: int = 1,
):
  result = subprocess.run(
    ["samtools", "view", "-H", bam],
    capture_output=True, text=True, check=True,
  )
  is_sorted = any("SO:coordinate" in line for line in result.stdout.splitlines())

  if is_sorted:
    os.replace(bam, bam_out)
  else:
    subprocess.run(["samtools", "sort", "-@", str(threads), "-o", bam_out, bam], check=True)
    os.remove(bam)

  subprocess.run(["samtools", "index", bam_out], check=True)
  return bam_out


def run_alignment(
  fastq_r1: str,
  fastq_r2: str,
  group: str,
  out_prefix: str,
  index_prefix: str,
  min_length: int,
  max_length: int,
  cb_tag: Optional[str] = "CB",
  umi_tag: Optional[str] = "UB",
  threads: int = 1,
  java_mem: Optional[str] = None,
  picard_jar: Optional[str] = None,
):
  # Step 1: align + filter
  clean_bam = generate_clean_bam(
    fastq_r1=fastq_r1,
    fastq_r2=fastq_r2,
    group=group,
    bam_out=f"{out_prefix}.clean.bam",
    index_prefix=index_prefix,
    min_length=min_length,
    max_length=max_length,
    threads=threads,
  )

  # Step 2: extract CB/UMI from read name into BAM tags
  if (cb_tag is not None) or (umi_tag is not None):
    tagged_bam = extract_adapter_as_tags(
        bam=clean_bam,
        bam_out=f"{out_prefix}.tagged.bam",
        cb_tag=cb_tag,
        umi_tag=umi_tag,
    )
  else:
    tagged_bam = clean_bam

  # Step 3: deduplication
  dedup_bam = deduplicate_bam(
    bam=tagged_bam,
    bam_out=f"{out_prefix}.dedup.bam",
    cb_tag=cb_tag,
    umi_tag=umi_tag,
    threads=threads,
    java_mem=java_mem,
    picard_jar=picard_jar,
  )

  # Step 4: coordinate sort + index
  final_bam = sort_and_index_bam(
    bam=dedup_bam,
    bam_out=f"{out_prefix}.bam",
    threads=threads,
  )

  return final_bam

def align(
  fastq_dir: str,
  out_dir: str,
  ref_genome: str,
  has_cb: bool = False,
  has_umi: bool = False,
  min_length: int = 10,
  max_length: int = 800,
  picard_jar: Optional[str] = None,
  java_mem: Optional[str] = None,
  threads: Optional[str] = None
):
  # Create output dir
  os.makedirs(out_dir, exist_ok=True)

  # Get reference index
  ref_preparation = ref_prep() 
  index_prefix = ref_preparation.get_bowtie2_index(ref_genome)

  # Check adapter structure
  if has_cb:
    cb_tag = "CB"
  else:
    cb_tag = None
  if has_umi:
    umi_tag = "UB"
  else:
    umi_tag = None

  # Get resource setting
  if picard_jar is None:
    picard_jar = get_picard_jar_path()
  if threads is None:
    threads = get_cpu_core()
  if java_mem is None:
    java_mem = get_memory()
  
  # Run
  fastq_list = get_files_path(fastq_dir, ext = ".fastq.gz")
  prefix_dict = get_fastq_prefix(fastq_list)
  bam_list = []

  for prefix in prefix_dict.keys():
    print(f"Processing target: {prefix}")
    r1 = prefix_dict[prefix]['fwd']
    r2 = prefix_dict[prefix]['rev']
    out_prefix = os.path.join(out_dir, prefix)

    try:
      final_bam = run_alignment(fastq_r1=r1, fastq_r2=r2, group=prefix, out_prefix=out_prefix, index_prefix=index_prefix, min_length=min_length, max_length=max_length, cb_tag=cb_tag, umi_tag=umi_tag, threads=threads, java_mem=java_mem, picard_jar=picard_jar)
      bam_list.append(final_bam)
      print(f"Successfully processed {prefix}")
    except Exception as e:
      print(f"Error processing {prefix}: {str(e)}")
      return 
  return


HELP = """
Align preprocessed paired-end FASTQ files to a reference genome and generate deduplicated BAM files.

Required:
  -i, --input     Directory containing paired-end FASTQ files (.fastq.gz)
  -o, --output    Output directory where final BAM files will be written
  -g, --genome    Reference genome name (hg38 or mm10)

Optional:
  --cb            Reads contain cell barcodes (CB) encoded in the read name.
                  If set, deduplication becomes CB-aware.
  --umi           Reads contain UMIs encoded in the read name.
                  If set, deduplication becomes UMI-aware.

  --min-len       Minimum fragment length accepted by Bowtie2. Pairs with inferred insert size below this threshold are discarded.
                  (default: 10)
  --max-len       Maximum fragment length accepted by Bowtie2 (-X flag). Pairs with inferred insert size above this threshold are discarded.
                  (default: 800)

  --java-mem      Java memory for Picard
                  (default: automatically detect available memory)
  --picard-jar    Path to picard.jar 
                  (default: automatically detect the picard.jar in conda env)
  -j, --threads   Number of CPU threads for Bowtie2 and samtools operations
                  (default: automatically detect all available CPU cores)

Output:
  {out_dir}/{prefix}.bam

Description:
  Processes all FASTQ file pairs found in input directory through a four-step pipeline: alignment with Bowtie2, optional CB/UMI tag extraction, deduplication, and coordinate sorting with indexing. Intermediate files are removed at each step to minimize disk usage.

Usage Example:
  1) Regular Hiplex CUT&Tag
     multiEpiPrep align \\
       -i ./adapter \\
       -o ./bam \\
       -g hg38

  2) UMI-containing Hiplex CUT&Tag
     multiEpiPrep align \\
       -i ./adapter \\
       -o ./bam \\
       -g hg38 \\
       --umi
"""

def register_parser(parser):
  parser.add_argument("-i", "--input", required=True)
  parser.add_argument("-o", "--output", required=True)
  parser.add_argument("-g", "--genome", required=True)
  parser.add_argument('--cb', action='store_true')
  parser.add_argument('--umi', action='store_true')
  parser.add_argument("--min-len", type=int, default=10)
  parser.add_argument("--max-len", type=int, default=800)
  parser.add_argument("--java-mem", type=str, default=None)
  parser.add_argument("--picard-jar", type=str, default=None)
  parser.add_argument("-j", "--threads", type=int, default=None)
  parser.set_defaults(func=run)

def run(args):
  align(
    fastq_dir=args.input,
    out_dir=args.output,
    ref_genome=args.genome,
    has_cb=args.cb,
    has_umi=args.umi,
    min_length=args.min_len,
    max_length=args.max_len,
    java_mem=args.java_mem,
    picard_jar=args.picard_jar,
    threads=args.threads
  )
  