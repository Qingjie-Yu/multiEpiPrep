import os, re, gzip, math
import pandas as pd
from typing import List, Dict, Optional
from collections import deque
import threading
from multiprocessing import Pool, Manager
from .utils import get_files_path, get_fastq_prefix, get_cpu_core

class AdapterConfig:
  def __init__(self, tag_list, label_list, error_rate, laxity, not_found_symbol = "[Not Found]", bases = "ATCGN"):
    self.tag_list = tag_list
    self.label_list = label_list
    self.error_rate = error_rate
    self.laxity = laxity
    self.not_found_symbol = not_found_symbol
    self.bases = bases
    self.config_list = self.process_tag()

  def generate_seqs_with_mismatches(self, seq: str):
      mismatches = math.floor(len(seq) * self.error_rate)
      if mismatches <= 0:
        return []  # No mismatches allowed, return empty set

      bases = set(self.bases)
      seq_list = list(seq)
      n = len(seq_list)
      results = set()

      def dfs(i: int, used: int):
        if i == n:
          if used > 0:
            results.add("".join(seq_list))
          return
        original = seq_list[i]
        # Branch 1: keep the original base at position i
        dfs(i + 1, used)
        # Branch 2: substitute base at position i (only if budget remains)
        if used < mismatches:
          for b in bases:
            if b == original:
              continue
            seq_list[i] = b
            dfs(i + 1, used + 1)
          seq_list[i] = original  # Restore after backtracking

      dfs(0, 0)
      return list(results)
  
  def check_ambiguity(self, config: Dict[str, List[str]]):
    seqs = list(config.keys())
    for i, seq in enumerate(seqs[:-1]):
      for other_seq in seqs[i+1:]:
        if (set(config[seq]) | {seq}) & (set(config[other_seq]) | {other_seq}):
          print(f"Warning: Ambiguity detected between {seq} and {other_seq} with error rate {self.error_rate}")
          print(f"  {seq} variants: {config[seq]}")
          print(f"  {other_seq} variants: {config[other_seq]}")
  
  def process_tag(self):
    config_list = []
    for tag in self.tag_list:
      if isinstance(tag, int):
        config_list.append(tag)
        continue
      elif isinstance(tag, str):
        tag = [tag]
      
      config = {}
      for t in tag:
        if self.error_rate == 0 :
          config[t] = []  # No variants when error_rate is zero
        else:
          config[t] = self.generate_seqs_with_mismatches(t)
      self.check_ambiguity(config)
      config_list.append(config)
    return config_list
  
  def check_read_for_tag(self, read: str):
    queue = deque()
    queue.append((0, 0, 0, ""))  #Initial state: start of read, first tag, no advances used, empty output
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
        # Fixed-length tag: capture exactly `config` bases unconditionally
        if offset + config <= read_len:
          new_out = f"{out}-[{read[offset:offset + config]}]"
          queue.append((offset + config, config_idx + 1, 0, new_out))
          matched = True
      else:
        # Sequence-based tag: try exact matches first
        for seq in config.keys():
          if read.startswith(seq, offset):
            new_out = f"{out}-[{seq}]"
            queue.append((offset + len(seq), config_idx + 1, 0, new_out))
            matched = True
            break
        # If no exact match, try mismatch variants
        if not matched:
          for variants in config.values():
            for seq in variants:
              if read.startswith(seq, offset):
                new_out = f"{out}-[{seq}]"
                queue.append((offset + len(seq), config_idx + 1, 0, new_out))
                matched = True
                break
            if matched:
              break

      # No match at this offset: advance one base
      if not matched:
        new_out = f"{out}-{read[offset]}"
        queue.append((offset + 1, config_idx, advance_used + 1, new_out))
    
    return f"{self.not_found_symbol}"

def run_adapter_identification(fastq_r1_in: str, fastq_r2_in: str, r1_out: str, r2_out: str, config: AdapterConfig, record_idx: Optional[List[int]] = None, trim: bool = True, concat_paired_tags: bool = False, progress_queue=None):

  def open_text_maybe_gz(path):
    if str(path).endswith(".gz"):
      return gzip.open(path, "rt", encoding="ascii", newline="")
    return open(path, "rt", encoding="ascii", newline="")
  
  def read_fastq_records(path):
    with open_text_maybe_gz(path) as fi:
      while True:
        n = fi.readline().strip()
        if not n:
          break
        s = fi.readline().strip()
        p = fi.readline().strip()
        q = fi.readline().strip()
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

  PRINT_INTERVAL = 10000
  read_count, fail_count = 0,0
  prefix = os.path.basename(fastq_r1_in)
  print(f"Processing {prefix}")

  with gzip.open(r1_out, "wt", encoding="ascii", newline="") as w1, \
       gzip.open(r2_out, "wt", encoding="ascii", newline="") as w2:
    
    if concat_paired_tags:
      # --- Paired-tag mode: identify adapters on BOTH R1 and R2 ---
      for (n1, s1, p1, q1), (n2, s2, p2, q2) in zip(read_fastq_records(fastq_r1_in), read_fastq_records(fastq_r2_in)):
        read_count += 1
        if progress_queue and read_count % PRINT_INTERVAL == 0:
          progress_queue.put(f"  [{prefix}] {read_count:,} reads, failed: {fail_count:,}")

        adapter_str_1 = config.check_read_for_tag(s1).lstrip("-")
        adapter_str_2 = config.check_read_for_tag(s2).lstrip("-")
        if adapter_str_1 == config.not_found_symbol or adapter_str_2 == config.not_found_symbol:
          fail_count += 1
          continue
        
        tags_1 = re.findall(r'\[([^\]]+)\]', adapter_str_1)
        tags_2 = re.findall(r'\[([^\]]+)\]', adapter_str_2)
        tag_str = "::".join(
            f"{name}={tag1}-{tag2}"
            for i, (tag1, tag2, name) in enumerate(zip(tags_1, tags_2, config.label_list))
            if record_idx is not None and i in record_idx
        )
        new_n1 = add_both_adapter_to_name(n1, adapter_str_1, adapter_str_2, tag_str)
        new_n2 = add_both_adapter_to_name(n2, adapter_str_1, adapter_str_2, tag_str)
        if trim:
          new_s1, new_q1 = trim_adapter_from_read(s1, q1, adapter_str_1)
          new_s2, new_q2 = trim_adapter_from_read(s2, q2, adapter_str_2)
          w1.write(f"{new_n1}\n{new_s1}\n{p1}\n{new_q1}\n")
          w2.write(f"{new_n2}\n{new_s2}\n{p2}\n{new_q2}\n")
        else:
          w1.write(f"{new_n1}\n{s1}\n{p1}\n{q1}\n")
          w2.write(f"{new_n2}\n{s2}\n{p2}\n{q2}\n")

    else:
      # --- Single-tag mode: identify adapter on R1 only, pass R2 through unchanged ---
      for (n1, s1, p1, q1), (n2, s2, p2, q2) in zip(read_fastq_records(fastq_r1_in), read_fastq_records(fastq_r2_in)):
        read_count += 1
        if progress_queue and read_count % PRINT_INTERVAL == 0:
          progress_queue.put(f"  [{prefix}] {read_count:,} reads, failed: {fail_count:,}")

        adapter_str = config.check_read_for_tag(s1)
        adapter_str = adapter_str.lstrip("-")

        if adapter_str == config.not_found_symbol:
          fail_count += 1
          continue
        
        tags = re.findall(r'\[([^\]]+)\]', adapter_str)
        tag_str = "::".join(
           f"{name}={tag}"
           for i, (tag, name) in enumerate(zip(tags, config.label_list))
           if record_idx is not None and i in record_idx
        )
        new_n1 = add_adapter_to_name(n1, adapter_str, tag_str)
        new_n2 = add_adapter_to_name(n2, adapter_str, tag_str)
        if trim:
          new_s1, new_q1 = trim_adapter_from_read(s1, q1, adapter_str)
          w1.write(f"{new_n1}\n{new_s1}\n{p1}\n{new_q1}\n")
          w2.write(f"{new_n2}\n{s2}\n{p2}\n{q2}\n")
        else:
          w1.write(f"{new_n1}\n{s1}\n{p1}\n{q1}\n")
          w2.write(f"{new_n2}\n{s2}\n{p2}\n{q2}\n")

  msg = f"  [{prefix}] Done: {read_count:,} reads"
  if fail_count > 0:
      msg += f", failed: {fail_count:,} ({fail_count/read_count:.2%})"
  else:
      msg += ", all reads identified."
  if progress_queue:
      progress_queue.put(msg)
  
  return r1_out, r2_out

def adapter_identification(
  fastq_r1_in: str,
  fastq_r2_in: str,
  r1_out: str,
  r2_out: str,
  cb: Optional[int | str | List[str]] = None,
  sp: Optional[int | str | List[str]] = None,
  umi: Optional[int | str | List[str]] = None,
  linker: Optional[int | str | List[str]] = None,
  error_rate: float = 0.1,
  laxity: int = 1,
  trim: bool = True,
  concat_paired_tags: bool = False
):
  """
  Identify adapter structures in FASTQ files and annotate read names with identified tags.

  Example Input read name:
  @A00123:45:H3F7MDSX2:1:1101:10000:1000 1:N:0:ATCGTAGC
  Example Output read name:
  @A00123:45:H3F7MDSX2:1:1101:10000:1000::[GGGG]-[CCCCCC]::CB=GGGG::UMI=CCCCCC
  """
  # Basic information
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
      val_dedup = list(set(value))  # Remove duplicates
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

  # Identify adapter sequence
  config = AdapterConfig(tag_list, label_list, error_rate, laxity)
  run_adapter_identification(fastq_r1_in, fastq_r2_in, r1_out, r2_out, config, record_idx, trim, concat_paired_tags)
  return (r1_out, r2_out)

def identify_adapter(
  fastq_dir: str,
  out_dir: str,
  cb: Optional[int | str | List[str]] = None,
  sp: Optional[int | str | List[str]] = None,
  umi: Optional[int | str | List[str]] = None,
  linker: Optional[int | str | List[str]] = None,
  error_rate: float = 0.1,
  laxity: int = 0,
  threads: Optional[int] = None
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
      val_dedup = list(set(value))  # Remove duplicates
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

  # Get avaiable cores
  if threads is None:
    threads = min(get_cpu_core(), 16)

  # Run
  trim = True
  concat_paired_tags = True
  config = AdapterConfig(tag_list, label_list, error_rate, laxity)

  def listener(q):
    while True:
      msg = q.get()
      if msg == "DONE":
        break
      print(msg, flush=True)
  
  manager = Manager()
  progress_queue = manager.Queue()
  listener_thread = threading.Thread(target=listener, args=(progress_queue,))
  listener_thread.start()

  try:
    with Pool(processes=threads) as pool:
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
  finally:
    progress_queue.put("DONE")
    listener_thread.join()
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
  -e, --error-rate  Float mismatch rate allowed for sequence-based tag matching
                    Maximum mismatches are calculated as floor(len(tag) * error_rate)
                    (default: 0.1)
  -l, --laxity      Maximum number of bases allowed to skip when searching for the next tag
                     (default: 0)

Other options:
  -j, --threads     Number of FASTQ pairs to process in parallel
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
      -e 0 \\
      -l 0
"""

def register_parser(parser):
  parser.add_argument("-i", "--input", required=True)
  parser.add_argument("-o", "--output", required=True)
  parser.add_argument("--cb", default = None)
  parser.add_argument("--sp", default=None)
  parser.add_argument("--umi", default=None)
  parser.add_argument("--linker", default=None)
  parser.add_argument("-e", "--error-rate", type=float, default=0.1)
  parser.add_argument("-l", "--laxity", type=int, default=0)
  parser.add_argument("-j", "--threads", type=int, default=None)
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
    threads=args.threads
  )