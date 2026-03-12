import os, re, gzip, math
import subprocess
from typing import List, Dict, Optional
from collections import deque
from multiprocessing import Pool
from .utils import get_files_path, get_fastq_prefix, get_cpu_core

class AdapterConfig:
  def __init__(self, tag_list, label_list, error_rate, laxity, not_found_symbol="[Not Found]", bases="ATCGN"):
    self.tag_list = tag_list
    self.label_list = label_list
    self.error_rate = error_rate
    self.laxity = laxity
    self.not_found_symbol = not_found_symbol
    self.bases = bases
    self.config_list = self.process_tag()

  def _hamming_mismatches(self, s1: str, s2: str) -> int:
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))
  
  def check_ambiguity(self, config: Dict[str, int]):
    seqs = list(config.keys())
    for i, seq1 in enumerate(seqs[:-1]):
      mm1 = config[seq1]
      for seq2 in seqs[i+1:]:
        if len(seq1) != len(seq2):
          continue
        mm2 = config[seq2]
        dist = self._hamming_mismatches(seq1, seq2)
        if dist <= mm1 + mm2:
          print(
            f"Warning: potential ambiguity detected between candidate "
            f"sequences:\n"
            f"  - seq1: '{seq1}' (length={len(seq1)}, max_mismatch={mm1})\n"
            f"  - seq2: '{seq2}' (length={len(seq2)}, max_mismatch={mm2})\n"
            f"  - Hamming distance: {dist}\n"
            f"  - Overlap threshold (mm1 + mm2): {mm1 + mm2}\n"
            f"  Interpretation: because distance <= allowed mismatch budget sum, "
            f"these two candidates may match overlapping observed sequences.\n"
            f"  Consequence: read assignment may depend on candidate order in the "
            f"matching loop.\n"
            f"  Suggestion: reduce error_rate, remove one of the candidates, or "
            f"use a more specific sequence set."
          )

  def process_tag(self):
    config_list = []
    for tag in self.tag_list:
      if isinstance(tag, int):
        config_list.append(tag)
        continue
      elif isinstance(tag, str):
        tag = [tag]

      config = {}
      seen = set()
      for t in tag:
        if t in seen:
          print(
            f"Warning: duplicate candidate sequence detected: '{t}'. "
            f"This sequence appears more than once in the candidate list."
          )
        else:
          seen.add(t)
          config[t] = math.floor(len(t) * self.error_rate)

      self.check_ambiguity(config)
      config_list.append(config)
    return config_list
  
  def check_read_for_tag(self, read: str):
    queue = deque()
    queue.append((0, 0, 0, ""))
    read_len = len(read)

    while queue:
      offset, config_idx, advance_used, out = queue.popleft()

      if config_idx >= len(self.config_list):
        return out

      if offset >= read_len:
        continue

      if advance_used > self.laxity:
        continue

      config = self.config_list[config_idx]
      matched = False

      if isinstance(config, int):
        if offset + config <= read_len:
          captured = read[offset:offset + config]
          new_out = f"{out}-[{captured}]"
          queue.append((offset + config, config_idx + 1, 0, new_out))
          matched = True
      else:
        # exact match first
        for seq, max_mm in config.items():
          if read.startswith(seq, offset):
            new_out = f"{out}-[{seq}]"
            queue.append((offset + len(seq), config_idx + 1, 0, new_out))
            matched = True
            break

        # fuzzy match if exact failed
        if not matched:
          for seq, max_mm in config.items():
            if max_mm <= 0:
              continue
            end = offset + len(seq)
            if end > read_len:
              continue
            sub = read[offset:end]
            mm = self._hamming_mismatches(sub, seq)
            if mm <= max_mm:
              # keep the actually observed sequence segment, for downstream trim / annotation
              new_out = f"{out}-[{sub}]"
              queue.append((end, config_idx + 1, 0, new_out))
              matched = True
              break

      if not matched:
        new_out = f"{out}-{read[offset]}"
        queue.append((offset + 1, config_idx, advance_used + 1, new_out))

    return self.not_found_symbol
  
def run_adapter_identification(fastq_r1_in: str, fastq_r2_in: str, r1_out: str, r2_out: str, config: AdapterConfig, record_idx: Optional[List[int]] = None, trim: bool = True, concat_paired_tags: bool = False, buffer_size = 8 * 1024 * 1024):

  def open_text_maybe_gz(path):
    if str(path).endswith(".gz"):
      return gzip.open(path, "rt", encoding="ascii", newline="")
    return open(path, "rt", encoding="ascii", newline="")

  def open_write_gz(path):
    outfile = open(path, 'wb')
    p = subprocess.Popen(
        ["pigz", "-p", "1", "-c", "-1"],
        stdin=subprocess.PIPE,
        stdout=outfile
    )
    outfile.close()
    return p

  def read_fastq_records(path):
    with open_text_maybe_gz(path) as fi:
      while True:
        n = fi.readline().rstrip("\n")
        if not n:
          break
        s = fi.readline().rstrip("\n")
        p = fi.readline().rstrip("\n")
        q = fi.readline().rstrip("\n")
        if not s or not p or not q:
          raise ValueError(f"Incomplete FASTQ record in {path}")
        yield n, s, p, q

  def add_adapter_to_name(n, adapter_str, tag_str):
    nc = n.split(" ", 1)[0]
    return f"{nc}::{adapter_str}::{tag_str}"

  def add_both_adapter_to_name(n, adapter_str_1, adapter_str_2, tag_str):
    nc = n.split(" ", 1)[0]
    return f"{nc}::{adapter_str_1}::{adapter_str_2}::{tag_str}"

  def trim_adapter_from_read(s, q, adapter_str):
    trim_len = sum(c.isalpha() for c in adapter_str)
    return s[trim_len:], q[trim_len:]

  read_count, fail_count = 0, 0
  prefix = os.path.basename(fastq_r1_in)
  proc1 = open_write_gz(r1_out)
  proc2 = open_write_gz(r2_out)

  try:
    buf1, buf2 = [], []
    size1 = 0
    if concat_paired_tags:
      # --- Paired-tag mode: identify adapters on BOTH R1 and R2 ---
      for (n1, s1, p1, q1), (n2, s2, p2, q2) in zip(read_fastq_records(fastq_r1_in), read_fastq_records(fastq_r2_in)):
        read_count += 1

        adapter_str_1 = config.check_read_for_tag(s1).lstrip("-")
        adapter_str_2 = config.check_read_for_tag(s2).lstrip("-")
        if adapter_str_1 == config.not_found_symbol or adapter_str_2 == config.not_found_symbol:
          fail_count += 1
          continue
        
        if record_idx is not None:
          tags_1 = re.findall(r'\[([^\]]+)\]', adapter_str_1)
          tags_2 = re.findall(r'\[([^\]]+)\]', adapter_str_2)
          tag_str = "::".join(
              f"{name}={tag1}-{tag2}"
              for i, (tag1, tag2, name) in enumerate(zip(tags_1, tags_2, config.label_list))
              if i in record_idx
          )
        else:
          tag_str = ""
        new_n1 = add_both_adapter_to_name(
            n1, adapter_str_1, adapter_str_2, tag_str)
        new_n2 = add_both_adapter_to_name(
            n2, adapter_str_1, adapter_str_2, tag_str)
        if trim:
          new_s1, new_q1 = trim_adapter_from_read(s1, q1, adapter_str_1)
          new_s2, new_q2 = trim_adapter_from_read(s2, q2, adapter_str_2)
          rec1 = f"{new_n1}\n{new_s1}\n{p1}\n{new_q1}\n"
          rec2 = f"{new_n2}\n{new_s2}\n{p2}\n{new_q2}\n"
        else:
          rec1 = f"{new_n1}\n{s1}\n{p1}\n{q1}\n"
          rec2 = f"{new_n2}\n{s2}\n{p2}\n{q2}\n"
        buf1.append(rec1)
        buf2.append(rec2)

        size1 += len(rec1)
        if size1 > buffer_size:
          proc1.stdin.write("".join(buf1).encode())
          proc2.stdin.write("".join(buf2).encode())
          buf1.clear()
          buf2.clear()
          size1 = 0

      if buf1:
        proc1.stdin.write("".join(buf1).encode())
        proc2.stdin.write("".join(buf2).encode())

    else:
      # --- Single-tag mode: identify adapter on R1 only, pass R2 through unchanged ---
      for (n1, s1, p1, q1), (n2, s2, p2, q2) in zip(read_fastq_records(fastq_r1_in), read_fastq_records(fastq_r2_in)):
        read_count += 1

        adapter_str = config.check_read_for_tag(s1)
        adapter_str = adapter_str.lstrip("-")

        if adapter_str == config.not_found_symbol:
          fail_count += 1
          continue
        
        if record_idx is not None:
          tags = re.findall(r'\[([^\]]+)\]', adapter_str)
          tag_str = "::".join(
              f"{name}={tag}"
              for i, (tag, name) in enumerate(zip(tags, config.label_list))
              if i in record_idx
          )
        else:
          tag_str = ""

        new_n1 = add_adapter_to_name(n1, adapter_str, tag_str)
        new_n2 = add_adapter_to_name(n2, adapter_str, tag_str)
        if trim:
          new_s1, new_q1 = trim_adapter_from_read(s1, q1, adapter_str)
          rec1 = f"{new_n1}\n{new_s1}\n{p1}\n{new_q1}\n"
          rec2 = f"{new_n2}\n{s2}\n{p2}\n{q2}\n"
        else:
          rec1 = f"{new_n1}\n{s1}\n{p1}\n{q1}\n"
          rec2 = f"{new_n2}\n{s2}\n{p2}\n{q2}\n"
        buf1.append(rec1)
        buf2.append(rec2)

        size1 += len(rec1)
        if size1 > buffer_size:
          proc1.stdin.write("".join(buf1).encode())
          proc2.stdin.write("".join(buf2).encode())
          buf1.clear()
          buf2.clear()
          size1 = 0

      if buf1:
        proc1.stdin.write("".join(buf1).encode())
        proc2.stdin.write("".join(buf2).encode())

    proc1.stdin.close()
    proc2.stdin.close()

  finally:
    if proc1.stdin and not proc1.stdin.closed:
      proc1.stdin.close()
    if proc2.stdin and not proc2.stdin.closed:
      proc2.stdin.close()
    rc1 = proc1.wait()
    rc2 = proc2.wait()
    if rc1 != 0:
        raise RuntimeError(f"pigz failed for {r1_out} with exit code {rc1}")
    if rc2 != 0:
        raise RuntimeError(f"pigz failed for {r2_out} with exit code {rc2}")

  msg = f"  [{prefix}] Done: {read_count:,} reads"
  if read_count == 0:
    msg += ", no reads found."
  elif fail_count > 0:
    msg += f", failed: {fail_count:,} ({fail_count/read_count:.2%})"
  else:
    msg += ", all reads identified."
  print(msg)

  return r1_out, r2_out

def identify_adapter(
  fastq_dir: str,
  out_dir: str,
  cb: Optional[int | str | List[str]] = None,
  sp: Optional[int | str | List[str]] = None,
  umi: Optional[int | str | List[str]] = None,
  linker: Optional[int | str | List[str]] = None,
  error_rate: float = 0.1,
  laxity: int = 0,
  threads: Optional[int] = None,
  exclude: List[str] = ["unknown", "IgG_control"]
):
  # Basic info
  label_list = []
  tag_list = []
  for label, value in [
      ("CB", cb),
      ("SPACER", sp),
      ("UMI", umi),
      ("LINKER", linker)
  ]:
    if value is None:
      continue
    label_list.append(label)

    if isinstance(value, int):
      tag_list.append(value)
      print(f"{label}: fixed length = {value}")
    elif isinstance(value, str):
      tag_list.append(value)
      print(f"{label}: fixed sequence = {value}")
    elif isinstance(value, list) and all(isinstance(x, str) for x in value):
      val_dedup = list(dict.fromkeys(value))  # Remove duplicates
      tag_list.append(val_dedup)
      print(f"{label}: candidate sequences = {','.join(val_dedup)}")
  
  if not label_list:
      raise ValueError("No adapter structure information provided.")
  label_str = '-'.join(f"[{name}]" for name in label_list)
  print(f"Adapter structure: {label_str}")

  # Determine which tag indices correspond to CB and UMI for downstream record keeping.
  record_idx = []
  for i, label in enumerate(label_list):
    if label == "CB":
      record_idx.append(i)
    elif label == "UMI":
      record_idx.append(i)
  record_idx = record_idx if record_idx else None

  # Create output dir
  os.makedirs(out_dir, exist_ok=True)

  # Extract .fastq.gz
  fastq_list = get_files_path(fastq_dir, ext = ".fastq.gz")
  prefix_dict = get_fastq_prefix(fastq_list)
  if not prefix_dict:
    raise ValueError("No valid FASTQ pairs were detected from input fastq_files.")
  
  if exclude is not None:
    keys = list(prefix_dict.keys())
    has_dash = ["-" in prefix for prefix in keys]
    if all(has_dash):
      keys_to_remove = [prefix for prefix in keys if any(p in exclude for p in prefix.split("-"))]
    elif not any(has_dash):
      keys_to_remove = [prefix for prefix in keys if prefix in exclude]
    else:
      keys_to_remove = [prefix for prefix in keys if any(item in prefix for item in exclude)]
    for prefix in keys_to_remove:
      prefix_dict.pop(prefix, None)
  if not prefix_dict:
    raise ValueError("No FASTQ pairs remained after applying exclude filter.")
  
  # Get avaiable cores
  if threads is None:
    threads = get_cpu_core()
  core_per_task = 3
  workers = max(1, min(len(prefix_dict), int(threads // core_per_task)))

  # Run
  trim = True
  concat_paired_tags = True
  config = AdapterConfig(tag_list, label_list, error_rate, laxity)

  with Pool(processes=workers) as pool:
    trimmed = []
    async_results = []
    for prefix in prefix_dict.keys():
      r1 = prefix_dict[prefix]['fwd']
      r2 = prefix_dict[prefix]['rev']
      r1_out = os.path.join(out_dir, f"{prefix}_R1.trimmed.fastq.gz")
      r2_out = os.path.join(out_dir, f"{prefix}_R2.trimmed.fastq.gz")
      ar = pool.apply_async(run_adapter_identification, args=(r1, r2, r1_out, r2_out, config, record_idx, trim, concat_paired_tags))
      async_results.append((prefix, ar))
  
    for prefix, ar in async_results:
      try:
        trimmed.extend(list(ar.get()))
      except Exception as e:
        print(f"Failed: {prefix} -> {e}")
        pool.terminate()
        break
  return


HELP = """
Identify adapter structures in paired-end FASTQ files and annotate read names with parsed tag information.

Required:
  -i, --input     Directory containing merged demultiplexed FASTQ files (.fastq.gz)
  -o, --output    Output directory for trimmed FASTQ files

Adapter structure options:
  --cb             Cell barcode definition
  --sp             Spacer definition
  --umi            UMI definition 
  --linker         Linker definition
                   Supported input formats for all four options:
                     1) Integer: fixed length
                        Example: --umi 8
                     2) String: fixed sequence
                        Example: --linker GCGATCGAGGACGGCAGATGTGTATAAGAGACAG
                     3) Comma-separated sequences: candidate sequence list
                        Example: --cb AAAA,CCCC,GGGG,TTTT
                     4) @file: one candidate sequence per line
                        Example: --cb @barcodes.txt

Matching options:
  -r, --error-rate  Float mismatch rate allowed for sequence-based tag matching
                    Maximum mismatches are calculated as floor(len(tag) * error_rate)
                    (default: 0.1)
  -l, --laxity      Maximum number of bases allowed to skip when searching for the next tag
                     (default: 0)

Other options:
  -e, --exclude     List of FASTQ prefixes to skip before adapter identification
                    (default: ["unknown", "IgG_control"])
  -j, --threads     Number of FASTQ prefixes to process in parallel
                    (default: automatically detect all available CPU cores)
Output:
  OUT_DIR/{prefix}_R1.trimmed.fastq.gz
  OUT_DIR/{prefix}_R2.trimmed.fastq.gz

Description:
  This command scans paired-end FASTQ files for a user-defined adapter structure composed of CB, SPACER, UMI, and LINKER components.

  [cell barcode] - [spacer] - [umi] - [linker]  

  Identified adapter information from both R1 and R2 is appended to the read name. CB and UMI values are additionally recorded as explicit key-value fields in the read name for downstream extraction.

Read name example:
  Input:
    @A00123:45:H3F7MDSX2:1:1101:10000:1000 1:N:0:ATCGTAGC
  Output:
    @A00123:45:H3F7MDSX2:1:1101:10000:1000::[GGGG]-[CCCCCC]::[AAAA]-[TTTTTT]::CB=GGGG-AAAA::UMI=CCCCCC-TTTTTT

Usage Example:
  1) Regular Hiplex CUT&Tag
     multiEpiPrep adapt \\
       -i ./demux \\
       -o ./adapter \\
       --linker GCGATCGAGGACGGCAGATGTGTATAAGAGACAG

  2) UMI-containing Hiplex CUT&Tag
     multiEpiPrep adapt \\
       -i ./demux \\
       -o ./adapter \\
       --sp 8 \\
       --umi 8 \\
       --linker GCGATCGAGGACGGCAGATGTGTATAAGAGACAG

  3) Candidate CB sequence list
     multiEpiPrep adapt \\
       -i ./demux \\
       -o ./adapter \\
       --cb AAAA,CCCC,GGGG,TTTT \\
       --umi 8 \\
       --linker GCGATCGAGGACGGCAGATGTGTATAAGAGACAG
  
  4) Strict Mode ===
    multiEpiPrep adapt \\
      -i ./demux \\
      -o ./adapt \\
      -r 0 \\
      -l 0
"""

def register_parser(parser):
  parser.add_argument("-i", "--input", required=True)
  parser.add_argument("-o", "--output", required=True)
  parser.add_argument("--cb", default = None)
  parser.add_argument("--sp", default=None)
  parser.add_argument("--umi", default=None)
  parser.add_argument("--linker", default=None)
  parser.add_argument("-r", "--error-rate", type=float, default=0.1)
  parser.add_argument("-l", "--laxity", type=int, default=0)
  parser.add_argument("-j", "--threads", type=int, default=None)
  parser.add_argument("-e", "--exclude", default="unknown,IgG_control")
  parser.set_defaults(func=run)

def run(args):
  def parse_tag_arg(value: str) -> "int | str | List[str] | None":
    """
    Unified parser for adapter component arguments.
    None       -> None       (component not provided)
    "8"        -> int 8      (fixed length)
    "@file"    -> List[str]  (sequences read from file, one per line)
    "A,B,C"    -> List[str]  (inline comma-separated sequences)
    "ATCG"     -> str        (single fixed sequence)
    """
    if value is None:
      return None
    value = value.strip()
    if not value:
      return None
    
    if value.isdigit():
      return int(value)
    if value.startswith("@"):
      path = value[1:]
      if not os.path.isfile(path):
        raise FileNotFoundError(f"Sequence file not found: {path}")
      with open(path, "r") as f:
          seqs = [line.strip() for line in f if line.strip()]
      if not seqs:
        raise ValueError(f"Sequence file is empty: {path}")
      return seqs
    
    if "," in value:
        seqs = [s.strip() for s in value.split(",") if s.strip()]
        if not seqs:
            raise ValueError(f"Could not parse any sequences from: {value}")
        return seqs
    return value
    
  identify_adapter(
    fastq_dir=args.input,
    out_dir=args.output,
    cb=parse_tag_arg(args.cb),
    sp=parse_tag_arg(args.sp),
    umi=parse_tag_arg(args.umi),
    linker=parse_tag_arg(args.linker),
    error_rate=args.error_rate,
    laxity=args.laxity,
    threads=args.threads,
    exclude=[x.strip() for x in args.exclude.split(",") if x.strip()]
  )