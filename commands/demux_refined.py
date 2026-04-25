#!/usr/bin/env python3
"""
================================================================================
FAST DNA BARCODE DEMULTIPLEXER
================================================================================

A high-performance tool for demultiplexing paired-end FASTQ files based on
DNA barcodes. Optimized for speed using hash-based O(1) barcode lookup and
parallel processing.

================================================================================
"""

import gzip
import json
import os
import sys
from collections import defaultdict, OrderedDict
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from multiprocessing import Pool, cpu_count
from itertools import combinations, product
import re
import time
import argparse


# =============================================================================
# VERSION
# =============================================================================

__version__ = "2.1.0"


# =============================================================================
# DATA STRUCTURES
# =============================================================================

@dataclass
class BarcodeConfig:
    """
    Configuration for a single barcode.

    Attributes:
        barcode: The DNA barcode sequence (e.g., "ATCACG")
        name: Human-readable name for this barcode (e.g., "H3K4me3")
    """
    barcode: str
    name: str


# =============================================================================
# FILE UTILITIES
# =============================================================================

def open_file(filepath: str, mode: str = 'r'):
    """
    Open a file, automatically handling gzip compression.
    """
    filepath = str(filepath)
    is_gzipped = filepath.endswith('.gz')

    if mode == 'r':
        if is_gzipped:
            return gzip.open(filepath, 'rt', encoding='utf-8')
        else:
            return open(filepath, 'r', encoding='utf-8')
    elif mode == 'w':
        if is_gzipped:
            return gzip.open(filepath, 'wt', encoding='utf-8')
        else:
            return open(filepath, 'w', encoding='utf-8')
    elif mode == 'a':
        if is_gzipped:
            return gzip.open(filepath, 'at', encoding='utf-8')
        else:
            return open(filepath, 'a', encoding='utf-8')
    else:
        raise ValueError(f"Unsupported mode: {mode}")


def get_file_size_str(filepath: str) -> str:
    """Get human-readable file size."""
    if not os.path.exists(filepath):
        return "N/A"
    size = os.path.getsize(filepath)
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if size < 1024:
            return f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} PB"


def is_gzipped(filepath: str) -> bool:
    """Check if file is gzipped."""
    return str(filepath).endswith('.gz')


def detect_file_format(filepath: str) -> dict:
    """Detect file format and properties."""
    filepath = str(filepath)

    info = {
        'path': filepath,
        'exists': os.path.exists(filepath),
        'is_gzipped': is_gzipped(filepath),
        'size': None,
        'size_str': None,
        'extension': None
    }

    if info['exists']:
        info['size'] = os.path.getsize(filepath)
        info['size_str'] = get_file_size_str(filepath)

    if filepath.endswith('.fastq.gz'):
        info['extension'] = '.fastq.gz'
    elif filepath.endswith('.fq.gz'):
        info['extension'] = '.fq.gz'
    elif filepath.endswith('.fastq'):
        info['extension'] = '.fastq'
    elif filepath.endswith('.fq'):
        info['extension'] = '.fq'
    else:
        info['extension'] = Path(filepath).suffix

    return info


def estimate_memory_usage(num_reads: int, num_barcodes: int = 100, max_mm: int = 1) -> dict:
    """Estimate memory usage for given parameters."""
    fastq_gb = (num_reads * 2 * 300) / (1024**3)
    output_gb = fastq_gb

    variants = {0: 1, 1: 19, 2: 154, 3: 694}.get(max_mm, 1000)
    hash_mb = (num_barcodes * 2 * variants * 100) / (1024**2)

    total_gb = fastq_gb + output_gb + (hash_mb / 1024)

    return {
        'fastq_gb': fastq_gb,
        'output_gb': output_gb,
        'hash_mb': hash_mb,
        'total_gb': total_gb,
        'recommended_mode': 'streaming' if total_gb > 16 else 'standard'
    }


def calculate_buffer_size_for_memory(
    target_memory_gb: float,
    num_r1_barcodes: int,
    num_r2_barcodes: int,
    bytes_per_record: int = 300,
    active_pair_ratio: float = 0.3
) -> int:
    """
    Calculate buffer size to achieve target memory usage.
    """
    target_bytes = target_memory_gb * (1024 ** 3)

    # Reserve 500MB for hash maps, Python overhead, etc.
    reserved_bytes = 500 * 1024 * 1024
    available_bytes = max(target_bytes - reserved_bytes, 100 * 1024 * 1024)

    # Estimate number of active pairs
    total_pairs = num_r1_barcodes * num_r2_barcodes + 1
    active_pairs = max(int(total_pairs * active_pair_ratio), 100)

    # Calculate buffer size
    buffer_size = int(available_bytes / (active_pairs * bytes_per_record))

    # Set reasonable bounds (min 1K, max 10M)
    buffer_size = max(1_000, min(buffer_size, 10_000_000))

    return buffer_size


# =============================================================================
# HASH MAP BUILDER
# =============================================================================

class BarcodeHashMap:
    """
    Pre-computed hash map for O(1) barcode lookup.
    """

    __slots__ = ['lookup_table', 'barcode_len', 'max_mismatches']

    BASES = 'ACGTN'

    def __init__(self, barcodes: List[BarcodeConfig], max_mismatches: int = 1):
        self.lookup_table: Dict[str, Tuple[str, int]] = {}
        self.barcode_len = len(barcodes[0].barcode) if barcodes else 0
        self.max_mismatches = max_mismatches

        self._build(barcodes, max_mismatches)

    def _generate_variants(self, barcode: str, max_mismatches: int) -> List[Tuple[str, int]]:
        """Generate all variants of a barcode up to max_mismatches."""
        variants = [(barcode, 0)]

        if max_mismatches == 0:
            return variants

        barcode_len = len(barcode)

        for num_mismatches in range(1, max_mismatches + 1):
            for positions in combinations(range(barcode_len), num_mismatches):
                original_bases = [barcode[pos] for pos in positions]

                for replacement_bases in product(self.BASES, repeat=num_mismatches):
                    if all(orig == repl for orig, repl in zip(original_bases, replacement_bases)):
                        continue

                    if any(orig != repl for orig, repl in zip(original_bases, replacement_bases)):
                        variant_list = list(barcode)
                        actual_mismatches = 0

                        for pos, new_base in zip(positions, replacement_bases):
                            if variant_list[pos] != new_base:
                                actual_mismatches += 1
                            variant_list[pos] = new_base

                        if actual_mismatches > 0:
                            variant = ''.join(variant_list)
                            variants.append((variant, actual_mismatches))

        return variants

    def _build(self, barcodes: List[BarcodeConfig], max_mismatches: int):
        """Build the hash map with all variants."""
        print(f"  Building hash map with max {max_mismatches} mismatch(es)...")

        total_variants = 0
        collisions = 0

        for config in barcodes:
            variants = self._generate_variants(config.barcode, max_mismatches)

            for variant, distance in variants:
                if variant in self.lookup_table:
                    existing_name, existing_distance = self.lookup_table[variant]

                    if distance < existing_distance:
                        self.lookup_table[variant] = (config.name, distance)
                    elif distance == existing_distance and existing_name != config.name:
                        collisions += 1
                        del self.lookup_table[variant]
                else:
                    self.lookup_table[variant] = (config.name, distance)
                    total_variants += 1

        print(f"  Total variants in hash map: {total_variants:,}")
        if collisions > 0:
            print(f"  Ambiguous variants removed: {collisions:,}")

    def lookup(self, sequence: str) -> Optional[Tuple[str, int]]:
        """O(1) barcode lookup."""
        candidate = sequence[:self.barcode_len]
        return self.lookup_table.get(candidate)


# =============================================================================
# FASTQ PARSER
# =============================================================================

def read_fastq_records(filepath: str, progress_callback=None) -> List[Tuple[str, str, str]]:
    """Read all FASTQ records from a file into memory."""
    records = []
    record_count = 0

    file_info = detect_file_format(filepath)
    print(f"  Reading: {filepath}")
    print(f"  Format:  {'gzipped' if file_info['is_gzipped'] else 'plain text'}")
    print(f"  Size:    {file_info['size_str']}")

    start_time = time.time()

    with open_file(filepath, 'r') as f:
        while True:
            header = f.readline()
            if not header:
                break
            seq = f.readline().strip()
            f.readline()  # Skip +
            qual = f.readline().strip()

            records.append((header.strip(), seq, qual))
            record_count += 1

            if progress_callback and record_count % 500_000 == 0:
                elapsed = time.time() - start_time
                rate = record_count / elapsed
                progress_callback(record_count, rate)

    elapsed = time.time() - start_time
    print(f"  Records: {record_count:,}")
    print(f"  Time:    {elapsed:.1f}s ({record_count/elapsed:,.0f} records/sec)")

    return records


# =============================================================================
# WORKER FUNCTION FOR PARALLEL PROCESSING
# =============================================================================

def process_chunk(args: Tuple) -> Tuple[Dict[str, int], Dict[str, List], Dict[str, int]]:
    """Process a chunk of read pairs."""
    (chunk, r1_lookup, r2_lookup, barcode_len, trim) = args

    stats = defaultdict(int)
    outputs = defaultdict(list)
    mismatch_stats = defaultdict(int)

    for r1_h, r1_s, r1_q, r2_h, r2_s, r2_q in chunk:

        r1_bc = r1_s[:barcode_len]
        r2_bc = r2_s[:barcode_len]

        r1_match = r1_lookup.get(r1_bc)
        r2_match = r2_lookup.get(r2_bc)

        if r1_match and r2_match:
            r1_name, r1_dist = r1_match
            r2_name, r2_dist = r2_match

            pair_name = f"{r1_name}-{r2_name}"
            stats[pair_name] += 1

            total_mismatches = r1_dist + r2_dist
            mismatch_stats[f'r1_mm_{r1_dist}'] += 1
            mismatch_stats[f'r2_mm_{r2_dist}'] += 1
            mismatch_stats[f'total_mm_{total_mismatches}'] += 1

            if trim:
                r1_s = r1_s[barcode_len:]
                r1_q = r1_q[barcode_len:]
                r2_s = r2_s[barcode_len:]
                r2_q = r2_q[barcode_len:]
        else:
            r1_label = r1_match[0] if r1_match else 'unknown'
            r2_label = r2_match[0] if r2_match else 'unknown'
            pair_name = f"{r1_label}-{r2_label}"
            stats['unmatched'] += 1

            if not r1_match:
                mismatch_stats['r1_no_match'] += 1
            if not r2_match:
                mismatch_stats['r2_no_match'] += 1

        outputs[pair_name].append((r1_h, r1_s, r1_q, r2_h, r2_s, r2_q))

    stats['total'] = len(chunk)

    return dict(stats), dict(outputs), dict(mismatch_stats)


# =============================================================================
# STANDARD IN-MEMORY DEMULTIPLEXER
# =============================================================================

class FastDemultiplexer:
    """
    High-performance parallel demultiplexer (in-memory mode).
    """

    def __init__(
        self,
        r1_barcodes: List[BarcodeConfig],
        r2_barcodes: List[BarcodeConfig],
        output_dir: str = "demultiplexed",
        num_workers: int = None,
        trim_barcode: bool = True,
        chunk_size: int = 100_000,
        r1_mismatches: int = 1,
        r2_mismatches: int = 1,
        gzip_output: bool = False,
        compression_level: int = 6
    ):
        self.output_dir = Path(output_dir)
        self.num_workers = num_workers or cpu_count()
        self.trim_barcode = trim_barcode
        self.chunk_size = chunk_size
        self.r1_mismatches = r1_mismatches
        self.r2_mismatches = r2_mismatches
        self.gzip_output = gzip_output
        self.compression_level = compression_level

        self.r1_barcodes = r1_barcodes
        self.r2_barcodes = r2_barcodes

        # Build hash maps
        print(f"\n{'='*70}")
        print("BUILDING BARCODE HASH MAPS")
        print(f"{'='*70}")
        print(f"R1 mismatches allowed: {self.r1_mismatches}")
        print(f"R2 mismatches allowed: {self.r2_mismatches}")
        print()

        print("Building R1 hash map...")
        self.r1_hashmap = BarcodeHashMap(r1_barcodes, self.r1_mismatches)

        print("\nBuilding R2 hash map...")
        self.r2_hashmap = BarcodeHashMap(r2_barcodes, self.r2_mismatches)

        self.barcode_len = self.r1_hashmap.barcode_len

        self.stats = defaultdict(int)
        self.mismatch_stats = defaultdict(int)

        self.output_dir.mkdir(parents=True, exist_ok=True)

        print(f"\n{'='*70}")
        print("DEMULTIPLEXER INITIALIZED (Standard Mode)")
        print(f"{'='*70}")
        print(f"R1 barcodes:          {len(r1_barcodes)}")
        print(f"R2 barcodes:          {len(r2_barcodes)}")
        print(f"R1 hash map entries:  {len(self.r1_hashmap.lookup_table):,}")
        print(f"R2 hash map entries:  {len(self.r2_hashmap.lookup_table):,}")
        print(f"Parallel workers:     {self.num_workers}")
        print(f"Chunk size:           {self.chunk_size:,}")
        print(f"Trim barcodes:        {self.trim_barcode}")
        print(f"Gzip output:          {self.gzip_output}")
        if self.gzip_output:
            print(f"Compression level:    {self.compression_level}")

    def _get_output_extension(self) -> str:
        return '.fastq.gz' if self.gzip_output else '.fastq'

    def demultiplex(self, r1_path: str, r2_path: str) -> Dict:
        start_time = time.time()

        r1_info = detect_file_format(r1_path)
        r2_info = detect_file_format(r2_path)

        print(f"\n{'='*70}")
        print("STARTING DEMULTIPLEXING (Standard Mode)")
        print(f"{'='*70}")
        print(f"\nINPUT FILES:")
        print(f"  R1: {r1_path}")
        print(f"      Format: {'gzipped' if r1_info['is_gzipped'] else 'plain text'}")
        print(f"      Size:   {r1_info['size_str']}")
        print(f"  R2: {r2_path}")
        print(f"      Format: {'gzipped' if r2_info['is_gzipped'] else 'plain text'}")
        print(f"      Size:   {r2_info['size_str']}")
        print(f"\nOUTPUT:")
        print(f"  Directory: {self.output_dir}")
        print(f"  Format:    {'gzipped (.fastq.gz)' if self.gzip_output else 'plain text (.fastq)'}")
        print(f"{'='*70}\n")

        # Load data
        print("Step 1/5: Loading R1 FASTQ file...")
        load_start = time.time()

        def progress_cb(count, rate):
            print(f"           {count:,} records ({rate:,.0f}/sec)")

        r1_records = read_fastq_records(r1_path, progress_cb)

        print("\n         Loading R2 FASTQ file...")
        r2_records = read_fastq_records(r2_path, progress_cb)

        total_reads = len(r1_records)
        load_time = time.time() - load_start

        if len(r1_records) != len(r2_records):
            print(f"\nWARNING: R1 ({len(r1_records):,}) and R2 ({len(r2_records):,}) have different record counts!")
            total_reads = min(len(r1_records), len(r2_records))

        print(f"\n         Total: {total_reads:,} read pairs loaded in {load_time:.1f}s")

        # Create chunks
        print("\nStep 2/5: Creating processing chunks...")

        chunks = []
        for i in range(0, total_reads, self.chunk_size):
            end = min(i + self.chunk_size, total_reads)
            chunk = [
                (r1_records[j][0], r1_records[j][1], r1_records[j][2],
                 r2_records[j][0], r2_records[j][1], r2_records[j][2])
                for j in range(i, end)
            ]
            chunks.append(chunk)

        print(f"         Created {len(chunks)} chunks of up to {self.chunk_size:,} reads each")

        del r1_records
        del r2_records

        worker_args = [
            (
                chunk,
                self.r1_hashmap.lookup_table,
                self.r2_hashmap.lookup_table,
                self.barcode_len,
                self.trim_barcode
            )
            for chunk in chunks
        ]

        print(f"\nStep 3/5: Processing with {self.num_workers} parallel workers...")
        process_start = time.time()

        all_outputs = defaultdict(list)

        with Pool(self.num_workers) as pool:
            results = pool.map(process_chunk, worker_args)

        process_time = time.time() - process_start
        print(f"         Completed in {process_time:.1f}s ({total_reads/process_time:,.0f} reads/sec)")

        print("\nStep 4/5: Merging results...")
        merge_start = time.time()

        for stats, outputs, mm_stats in results:
            for key, value in stats.items():
                self.stats[key] += value
            for pair_name, records in outputs.items():
                all_outputs[pair_name].extend(records)
            for key, value in mm_stats.items():
                self.mismatch_stats[key] += value

        merge_time = time.time() - merge_start
        print(f"         Completed in {merge_time:.1f}s")

        print(f"\nStep 5/5: Writing output files ({'gzipped' if self.gzip_output else 'plain text'})...")
        write_start = time.time()

        self._write_outputs(all_outputs)

        write_time = time.time() - write_start
        print(f"         Completed in {write_time:.1f}s")

        total_time = time.time() - start_time
        self._print_stats(total_time, load_time, process_time, write_time)
        self._save_report()

        return dict(self.stats)

    def _write_outputs(self, outputs: Dict[str, List]):
        ext = self._get_output_extension()
        num_pairs = len(outputs)

        for idx, (pair_name, records) in enumerate(outputs.items(), 1):
            r1_path = self.output_dir / f"{pair_name}_R1{ext}"
            r2_path = self.output_dir / f"{pair_name}_R2{ext}"

            if self.gzip_output:
                with gzip.open(r1_path, 'wt', compresslevel=self.compression_level) as r1_f, \
                     gzip.open(r2_path, 'wt', compresslevel=self.compression_level) as r2_f:

                    for r1_h, r1_s, r1_q, r2_h, r2_s, r2_q in records:
                        r1_f.write(f"{r1_h}\n{r1_s}\n+\n{r1_q}\n")
                        r2_f.write(f"{r2_h}\n{r2_s}\n+\n{r2_q}\n")
            else:
                buffer_size = 10 * 1024 * 1024
                with open(r1_path, 'w', buffering=buffer_size) as r1_f, \
                     open(r2_path, 'w', buffering=buffer_size) as r2_f:

                    for r1_h, r1_s, r1_q, r2_h, r2_s, r2_q in records:
                        r1_f.write(f"{r1_h}\n{r1_s}\n+\n{r1_q}\n")
                        r2_f.write(f"{r2_h}\n{r2_s}\n+\n{r2_q}\n")

            if idx % 100 == 0 or idx == num_pairs:
                print(f"         Written {idx}/{num_pairs} barcode pairs...")

    def _print_stats(self, total_time, load_time, process_time, write_time):
        total_reads = self.stats.get('total', 0)

        print(f"\n{'='*70}")
        print("DEMULTIPLEXING COMPLETE")
        print(f"{'='*70}")

        print(f"\n{'TIMING':^70}")
        print("-" * 70)
        print(f"  Loading:         {load_time:>10.1f}s")
        print(f"  Processing:      {process_time:>10.1f}s")
        print(f"  Writing:         {write_time:>10.1f}s")
        print(f"  Total:           {total_time:>10.1f}s")

        print(f"\n{'PERFORMANCE':^70}")
        print("-" * 70)
        reads_per_sec = total_reads / total_time if total_time > 0 else 0
        print(f"  Total reads:     {total_reads:>15,}")
        print(f"  Reads/second:    {reads_per_sec:>15,.0f}")

        print(f"\n{'RESULTS SUMMARY':^70}")
        print("-" * 70)
        undetermined = self.stats.get('unmatched', 0)
        matched = total_reads - undetermined

        pct_matched = 100 * matched / total_reads if total_reads > 0 else 0
        pct_undetermined = 100 * undetermined / total_reads if total_reads > 0 else 0

        print(f"  Matched:         {matched:>15,} ({pct_matched:>6.2f}%)")
        print(f"  Undetermined:    {undetermined:>15,} ({pct_undetermined:>6.2f}%)")

        print(f"\n{'MISMATCH STATISTICS':^70}")
        print("-" * 70)

        print("\n  R1 Mismatches:")
        for i in range(self.r1_mismatches + 1):
            key = f'r1_mm_{i}'
            count = self.mismatch_stats.get(key, 0)
            pct = 100 * count / total_reads if total_reads > 0 else 0
            bar = '█' * int(pct / 2)
            print(f"    {i} mismatch(es): {count:>12,} ({pct:>6.2f}%) {bar}")

        print("\n  R2 Mismatches:")
        for i in range(self.r2_mismatches + 1):
            key = f'r2_mm_{i}'
            count = self.mismatch_stats.get(key, 0)
            pct = 100 * count / total_reads if total_reads > 0 else 0
            bar = '█' * int(pct / 2)
            print(f"    {i} mismatch(es): {count:>12,} ({pct:>6.2f}%) {bar}")

        print(f"\n{'TOP 10 BARCODE PAIRS':^70}")
        print("-" * 70)

        pair_stats = {k: v for k, v in self.stats.items() if k not in ['total', 'unmatched']}
        sorted_pairs = sorted(pair_stats.items(), key=lambda x: -x[1])[:10]

        for pair_name, count in sorted_pairs:
            pct = 100 * count / total_reads if total_reads > 0 else 0
            bar = '█' * int(pct)
            print(f"  {pair_name:<35} {count:>10,} ({pct:>5.2f}%) {bar}")

        if len(pair_stats) > 10:
            print(f"\n  ... and {len(pair_stats) - 10} more pairs")

        ext = self._get_output_extension()
        print(f"\n{'OUTPUT FILES':^70}")
        print("-" * 70)
        print(f"  Directory:       {self.output_dir}")
        print(f"  Format:          {ext}")
        print(f"  Barcode pairs:   {len(pair_stats) + (1 if undetermined > 0 else 0)}")

        print(f"\n{'='*70}")

    def _save_report(self):
        report_path = self.output_dir / "demux_report.json"
        total_reads = self.stats.get('total', 0)

        report = {
            'version': __version__,
            'mode': 'standard',
            'settings': {
                'r1_mismatches': self.r1_mismatches,
                'r2_mismatches': self.r2_mismatches,
                'num_r1_barcodes': len(self.r1_barcodes),
                'num_r2_barcodes': len(self.r2_barcodes),
                'barcode_length': self.barcode_len,
                'trim_barcode': self.trim_barcode,
                'num_workers': self.num_workers,
                'chunk_size': self.chunk_size,
                'gzip_output': self.gzip_output,
                'compression_level': self.compression_level if self.gzip_output else None
            },
            'summary': {
                'total_reads': total_reads,
                'matched': total_reads - self.stats.get('unmatched', 0),
                'unmatched': self.stats.get('unmatched', 0)
            },
            'mismatch_stats': dict(self.mismatch_stats),
            'pair_counts': {k: v for k, v in self.stats.items() if k not in ['total']},
            'barcodes': {
                'r1': [asdict(bc) for bc in self.r1_barcodes],
                'r2': [asdict(bc) for bc in self.r2_barcodes]
            }
        }

        with open(report_path, 'w') as f:
            json.dump(report, f, indent=2)

        print(f"Report saved to: {report_path}")


# =============================================================================
# STREAMING DEMULTIPLEXER WITH LRU FILE HANDLE CACHE
# =============================================================================

class StreamingDemultiplexer:
    """
    Memory-efficient streaming demultiplexer with LRU file handle caching.
    """

    def __init__(
        self,
        r1_barcodes: List[BarcodeConfig],
        r2_barcodes: List[BarcodeConfig],
        output_dir: str = "demultiplexed",
        trim_barcode: bool = True,
        r1_mismatches: int = 1,
        r2_mismatches: int = 1,
        gzip_output: bool = False,
        compression_level: int = 6,
        buffer_size: int = 10_000,
        max_open_files: int = 200
    ):
        self.output_dir = Path(output_dir)
        self.trim_barcode = trim_barcode
        self.r1_mismatches = r1_mismatches
        self.r2_mismatches = r2_mismatches
        self.gzip_output = gzip_output
        self.compression_level = compression_level
        self.buffer_size = buffer_size
        self.max_open_files = max_open_files

        self.r1_barcodes = r1_barcodes
        self.r2_barcodes = r2_barcodes

        # Build hash maps
        print(f"\n{'='*70}")
        print("BUILDING BARCODE HASH MAPS")
        print(f"{'='*70}")
        print(f"R1 mismatches allowed: {self.r1_mismatches}")
        print(f"R2 mismatches allowed: {self.r2_mismatches}")
        print()

        print("Building R1 hash map...")
        self.r1_hashmap = BarcodeHashMap(r1_barcodes, self.r1_mismatches)

        print("\nBuilding R2 hash map...")
        self.r2_hashmap = BarcodeHashMap(r2_barcodes, self.r2_mismatches)

        self.barcode_len = self.r1_hashmap.barcode_len

        self.output_buffers: Dict[str, List] = defaultdict(list)
        self.output_handles: OrderedDict[str, Tuple] = OrderedDict()
        self.created_files: set = set()

        self.stats = defaultdict(int)
        self.mismatch_stats = defaultdict(int)
        self.file_handle_stats = {
            'opens': 0,
            'closes': 0,
            'cache_hits': 0,
            'cache_misses': 0
        }

        self.output_dir.mkdir(parents=True, exist_ok=True)

        num_possible_pairs = len(r1_barcodes) * len(r2_barcodes) + 1
        max_buffer_memory = num_possible_pairs * buffer_size * 300

        print(f"\n{'='*70}")
        print("DEMULTIPLEXER INITIALIZED (Streaming Mode)")
        print(f"{'='*70}")
        print(f"R1 barcodes:          {len(r1_barcodes)}")
        print(f"R2 barcodes:          {len(r2_barcodes)}")
        print(f"R1 hash map entries:  {len(self.r1_hashmap.lookup_table):,}")
        print(f"R2 hash map entries:  {len(self.r2_hashmap.lookup_table):,}")
        print(f"Buffer size:          {buffer_size:,} records per pair")
        print(f"Max open files:       {max_open_files} pairs ({max_open_files * 2} handles)")
        print(f"Max buffer memory:    {max_buffer_memory / 1e6:.1f} MB")
        print(f"Trim barcodes:        {self.trim_barcode}")
        print(f"Gzip output:          {self.gzip_output}")
        if self.gzip_output:
            print(f"Compression level:    {self.compression_level}")

    def _get_output_extension(self) -> str:
        return '.fastq.gz' if self.gzip_output else '.fastq'

    def _get_file_paths(self, pair_name: str) -> Tuple[Path, Path]:
        ext = self._get_output_extension()
        r1_path = self.output_dir / f"{pair_name}_R1{ext}"
        r2_path = self.output_dir / f"{pair_name}_R2{ext}"
        return r1_path, r2_path

    def _open_file_pair(self, pair_name: str) -> Tuple:
        r1_path, r2_path = self._get_file_paths(pair_name)

        if pair_name in self.created_files:
            mode = 'at' if self.gzip_output else 'a'
        else:
            mode = 'wt' if self.gzip_output else 'w'
            self.created_files.add(pair_name)

        if self.gzip_output:
            r1_f = gzip.open(r1_path, mode, compresslevel=self.compression_level)
            r2_f = gzip.open(r2_path, mode, compresslevel=self.compression_level)
        else:
            r1_f = open(r1_path, mode, buffering=1024*1024)
            r2_f = open(r2_path, mode, buffering=1024*1024)

        self.file_handle_stats['opens'] += 2

        return (r1_f, r2_f)

    def _close_file_pair(self, pair_name: str):
        if pair_name in self.output_handles:
            r1_f, r2_f = self.output_handles[pair_name]
            r1_f.close()
            r2_f.close()
            del self.output_handles[pair_name]
            self.file_handle_stats['closes'] += 2

    def _evict_lru_handle(self):
        if self.output_handles:
            oldest_pair = next(iter(self.output_handles))

            if oldest_pair in self.output_buffers and self.output_buffers[oldest_pair]:
                self._write_buffer_to_handles(oldest_pair)

            self._close_file_pair(oldest_pair)

    def _get_output_handles(self, pair_name: str) -> Tuple:
        if pair_name in self.output_handles:
            self.output_handles.move_to_end(pair_name)
            self.file_handle_stats['cache_hits'] += 1
            return self.output_handles[pair_name]

        self.file_handle_stats['cache_misses'] += 1

        while len(self.output_handles) >= self.max_open_files:
            self._evict_lru_handle()

        handles = self._open_file_pair(pair_name)
        self.output_handles[pair_name] = handles

        return handles

    def _write_buffer_to_handles(self, pair_name: str):
        if pair_name not in self.output_buffers:
            return

        records = self.output_buffers[pair_name]
        if not records:
            return

        r1_f, r2_f = self.output_handles[pair_name]

        for r1_h, r1_s, r1_q, r2_h, r2_s, r2_q in records:
            r1_f.write(f"{r1_h}\n{r1_s}\n+\n{r1_q}\n")
            r2_f.write(f"{r2_h}\n{r2_s}\n+\n{r2_q}\n")

        self.output_buffers[pair_name] = []

    def _flush_buffer(self, pair_name: str):
        if pair_name not in self.output_buffers:
            return

        records = self.output_buffers[pair_name]
        if not records:
            return

        self._get_output_handles(pair_name)
        self._write_buffer_to_handles(pair_name)

    def _flush_all_buffers(self):
        for pair_name in list(self.output_buffers.keys()):
            if self.output_buffers[pair_name]:
                self._flush_buffer(pair_name)

    def _close_all_handles(self):
        for pair_name in list(self.output_buffers.keys()):
            if self.output_buffers[pair_name]:
                if pair_name in self.output_handles:
                    self._write_buffer_to_handles(pair_name)
                else:
                    self._flush_buffer(pair_name)

        for pair_name in list(self.output_handles.keys()):
            self._close_file_pair(pair_name)

    def demultiplex(self, r1_path: str, r2_path: str) -> Dict:
        start_time = time.time()

        r1_info = detect_file_format(r1_path)
        r2_info = detect_file_format(r2_path)

        print(f"\n{'='*70}")
        print("STARTING DEMULTIPLEXING (Streaming Mode)")
        print(f"{'='*70}")
        print(f"\nINPUT FILES:")
        print(f"  R1: {r1_path}")
        print(f"      Format: {'gzipped' if r1_info['is_gzipped'] else 'plain text'}")
        print(f"      Size:   {r1_info['size_str']}")
        print(f"  R2: {r2_path}")
        print(f"      Format: {'gzipped' if r2_info['is_gzipped'] else 'plain text'}")
        print(f"      Size:   {r2_info['size_str']}")
        print(f"\nOUTPUT:")
        print(f"  Directory: {self.output_dir}")
        print(f"  Format:    {'gzipped (.fastq.gz)' if self.gzip_output else 'plain text (.fastq)'}")
        print(f"\nPROCESSING:")
        print(f"  Mode:           Streaming (low memory)")
        print(f"  Buffer:         {self.buffer_size:,} records per barcode pair")
        print(f"  Max open files: {self.max_open_files} pairs")
        print(f"{'='*70}\n")

        print("Processing reads...")

        try:
            with open_file(r1_path, 'r') as r1_f, open_file(r2_path, 'r') as r2_f:

                while True:
                    r1_h = r1_f.readline()
                    if not r1_h:
                        break
                    r1_s = r1_f.readline().strip()
                    r1_f.readline()
                    r1_q = r1_f.readline().strip()

                    r2_h = r2_f.readline()
                    if not r2_h:
                        break
                    r2_s = r2_f.readline().strip()
                    r2_f.readline()
                    r2_q = r2_f.readline().strip()

                    self.stats['total'] += 1

                    r1_match = self.r1_hashmap.lookup_table.get(r1_s[:self.barcode_len])
                    r2_match = self.r2_hashmap.lookup_table.get(r2_s[:self.barcode_len])

                    if r1_match and r2_match:
                        r1_name, r1_dist = r1_match
                        r2_name, r2_dist = r2_match
                        pair_name = f"{r1_name}-{r2_name}"

                        self.stats[pair_name] += 1
                        self.mismatch_stats[f'r1_mm_{r1_dist}'] += 1
                        self.mismatch_stats[f'r2_mm_{r2_dist}'] += 1
                        self.mismatch_stats[f'total_mm_{r1_dist + r2_dist}'] += 1

                        if self.trim_barcode:
                            r1_s = r1_s[self.barcode_len:]
                            r1_q = r1_q[self.barcode_len:]
                            r2_s = r2_s[self.barcode_len:]
                            r2_q = r2_q[self.barcode_len:]
                    else:
                        r1_label = r1_match[0] if r1_match else 'unknown'
                        r2_label = r2_match[0] if r2_match else 'unknown'
                        pair_name = f"{r1_label}-{r2_label}"
                        self.stats['unmatched'] += 1

                        if not r1_match:
                            self.mismatch_stats['r1_no_match'] += 1
                        if not r2_match:
                            self.mismatch_stats['r2_no_match'] += 1

                    self.output_buffers[pair_name].append(
                        (r1_h.strip(), r1_s, r1_q, r2_h.strip(), r2_s, r2_q)
                    )

                    if len(self.output_buffers[pair_name]) >= self.buffer_size:
                        self._flush_buffer(pair_name)

                    if self.stats['total'] % 1_000_000 == 0:
                        elapsed = time.time() - start_time
                        rate = self.stats['total'] / elapsed
                        open_files = len(self.output_handles) * 2
                        print(f"  Processed {self.stats['total']:,} reads "
                              f"({rate:,.0f}/sec, {open_files} open files)")

            print("\nFlushing remaining buffers...")
            self._flush_all_buffers()

        finally:
            print("Closing file handles...")
            self._close_all_handles()

        total_time = time.time() - start_time
        self._print_stats(total_time)
        self._save_report()

        return dict(self.stats)

    def _print_stats(self, total_time: float):
        total_reads = self.stats.get('total', 0)

        print(f"\n{'='*70}")
        print("DEMULTIPLEXING COMPLETE (Streaming Mode)")
        print(f"{'='*70}")

        print(f"\n{'TIMING':^70}")
        print("-" * 70)
        print(f"  Total time:      {total_time:>10.1f}s")

        print(f"\n{'PERFORMANCE':^70}")
        print("-" * 70)
        reads_per_sec = total_reads / total_time if total_time > 0 else 0
        print(f"  Total reads:     {total_reads:>15,}")
        print(f"  Reads/second:    {reads_per_sec:>15,.0f}")

        print(f"\n{'FILE HANDLE STATISTICS':^70}")
        print("-" * 70)
        print(f"  File opens:      {self.file_handle_stats['opens']:>15,}")
        print(f"  File closes:     {self.file_handle_stats['closes']:>15,}")
        print(f"  Cache hits:      {self.file_handle_stats['cache_hits']:>15,}")
        print(f"  Cache misses:    {self.file_handle_stats['cache_misses']:>15,}")

        cache_total = self.file_handle_stats['cache_hits'] + self.file_handle_stats['cache_misses']
        if cache_total > 0:
            hit_rate = 100 * self.file_handle_stats['cache_hits'] / cache_total
            print(f"  Cache hit rate:  {hit_rate:>14.1f}%")

        print(f"\n{'RESULTS SUMMARY':^70}")
        print("-" * 70)
        undetermined = self.stats.get('unmatched', 0)
        matched = total_reads - undetermined

        pct_matched = 100 * matched / total_reads if total_reads > 0 else 0
        pct_undetermined = 100 * undetermined / total_reads if total_reads > 0 else 0

        print(f"  Matched:         {matched:>15,} ({pct_matched:>6.2f}%)")
        print(f"  Undetermined:    {undetermined:>15,} ({pct_undetermined:>6.2f}%)")

        print(f"\n{'MISMATCH STATISTICS':^70}")
        print("-" * 70)

        print("\n  R1 Mismatches:")
        for i in range(self.r1_mismatches + 1):
            key = f'r1_mm_{i}'
            count = self.mismatch_stats.get(key, 0)
            pct = 100 * count / total_reads if total_reads > 0 else 0
            bar = '█' * int(pct / 2)
            print(f"    {i} mismatch(es): {count:>12,} ({pct:>6.2f}%) {bar}")

        print("\n  R2 Mismatches:")
        for i in range(self.r2_mismatches + 1):
            key = f'r2_mm_{i}'
            count = self.mismatch_stats.get(key, 0)
            pct = 100 * count / total_reads if total_reads > 0 else 0
            bar = '█' * int(pct / 2)
            print(f"    {i} mismatch(es): {count:>12,} ({pct:>6.2f}%) {bar}")

        print(f"\n{'TOP 10 BARCODE PAIRS':^70}")
        print("-" * 70)

        pair_stats = {k: v for k, v in self.stats.items() if k not in ['total', 'unmatched']}
        sorted_pairs = sorted(pair_stats.items(), key=lambda x: -x[1])[:10]

        for pair_name, count in sorted_pairs:
            pct = 100 * count / total_reads if total_reads > 0 else 0
            bar = '█' * int(pct)
            print(f"  {pair_name:<35} {count:>10,} ({pct:>5.2f}%) {bar}")

        if len(pair_stats) > 10:
            print(f"\n  ... and {len(pair_stats) - 10} more pairs")

        ext = self._get_output_extension()
        print(f"\n{'OUTPUT FILES':^70}")
        print("-" * 70)
        print(f"  Directory:       {self.output_dir}")
        print(f"  Format:          {ext}")
        print(f"  Barcode pairs:   {len(self.created_files)}")

        print(f"\n{'='*70}")

    def _save_report(self):
        report_path = self.output_dir / "demux_report.json"
        total_reads = self.stats.get('total', 0)

        report = {
            'version': __version__,
            'mode': 'streaming',
            'settings': {
                'r1_mismatches': self.r1_mismatches,
                'r2_mismatches': self.r2_mismatches,
                'num_r1_barcodes': len(self.r1_barcodes),
                'num_r2_barcodes': len(self.r2_barcodes),
                'barcode_length': self.barcode_len,
                'trim_barcode': self.trim_barcode,
                'buffer_size': self.buffer_size,
                'max_open_files': self.max_open_files,
                'gzip_output': self.gzip_output,
                'compression_level': self.compression_level if self.gzip_output else None
            },
            'summary': {
                'total_reads': total_reads,
                'matched': total_reads - self.stats.get('unmatched', 0),
                'unmatched': self.stats.get('unmatched', 0)
            },
            'file_handle_stats': self.file_handle_stats,
            'mismatch_stats': dict(self.mismatch_stats),
            'pair_counts': {k: v for k, v in self.stats.items() if k not in ['total']},
            'barcodes': {
                'r1': [asdict(bc) for bc in self.r1_barcodes],
                'r2': [asdict(bc) for bc in self.r2_barcodes]
            }
        }

        with open(report_path, 'w') as f:
            json.dump(report, f, indent=2)

        print(f"Report saved to: {report_path}")


# =============================================================================
# TEST DATA GENERATOR
# =============================================================================

def create_test_data(
    num_reads: int = 1_000_000,
    num_r1_barcodes: int = 100,
    num_r2_barcodes: int = 100,
    r1_mismatches: int = 1,
    r2_mismatches: int = 1,
    mismatch_rate: float = 0.1,
    output_dir: str = ".",
    gzip_output: bool = False
):
    """Create test FASTQ data for benchmarking."""
    import random

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    ext = '.fastq.gz' if gzip_output else '.fastq'

    print(f"\n{'='*70}")
    print("GENERATING TEST DATA")
    print(f"{'='*70}")
    print(f"  Reads:           {num_reads:,}")
    print(f"  R1 barcodes:     {num_r1_barcodes}")
    print(f"  R2 barcodes:     {num_r2_barcodes}")
    print(f"  R1 mismatches:   {r1_mismatches}")
    print(f"  R2 mismatches:   {r2_mismatches}")
    print(f"  Mismatch rate:   {mismatch_rate:.1%}")
    print(f"  Output format:   {'gzipped' if gzip_output else 'plain text'}")
    print(f"{'='*70}\n")

    def random_barcode(length=6):
        return ''.join(random.choices('ACGT', k=length))

    def introduce_mismatches(barcode: str, num_mm: int) -> str:
        if num_mm == 0:
            return barcode
        barcode_list = list(barcode)
        positions = random.sample(range(len(barcode)), min(num_mm, len(barcode)))
        for pos in positions:
            original = barcode_list[pos]
            alternatives = [b for b in 'ACGT' if b != original]
            barcode_list[pos] = random.choice(alternatives)
        return ''.join(barcode_list)

    r1_barcode_seqs = set()
    while len(r1_barcode_seqs) < num_r1_barcodes:
        r1_barcode_seqs.add(random_barcode())

    r2_barcode_seqs = set()
    while len(r2_barcode_seqs) < num_r2_barcodes:
        r2_barcode_seqs.add(random_barcode())

    r1_barcodes = [BarcodeConfig(barcode=bc, name=f"Ab{i:03d}") for i, bc in enumerate(sorted(r1_barcode_seqs))]
    r2_barcodes = [BarcodeConfig(barcode=bc, name=f"Sample{i:03d}") for i, bc in enumerate(sorted(r2_barcode_seqs))]

    # Save barcode files
    r1_barcode_file = output_path / "r1_barcodes.tsv"
    r2_barcode_file = output_path / "r2_barcodes.tsv"
    barcode_json_file = output_path / "barcodes.json"

    print("Saving barcode files...")

    with open(r1_barcode_file, 'w') as f:
        f.write("barcode\tname\n")
        for bc in r1_barcodes:
            f.write(f"{bc.barcode}\t{bc.name}\n")

    with open(r2_barcode_file, 'w') as f:
        f.write("barcode\tname\n")
        for bc in r2_barcodes:
            f.write(f"{bc.barcode}\t{bc.name}\n")

    barcode_data = {
        "r1_barcodes": [asdict(bc) for bc in r1_barcodes],
        "r2_barcodes": [asdict(bc) for bc in r2_barcodes],
        "metadata": {
            "num_reads": num_reads,
            "barcode_length": 6,
            "r1_mismatches": r1_mismatches,
            "r2_mismatches": r2_mismatches,
            "mismatch_rate": mismatch_rate,
            "generated_at": time.strftime("%Y-%m-%d %H:%M:%S")
        }
    }

    with open(barcode_json_file, 'w') as f:
        json.dump(barcode_data, f, indent=2)

    # Print barcodes
    print(f"\n{'R1 BARCODES':^70}")
    print("-" * 70)
    for bc in r1_barcodes[:10]:
        print(f"  {bc.barcode:<12} {bc.name}")
    if len(r1_barcodes) > 10:
        print(f"  ... and {len(r1_barcodes) - 10} more")

    print(f"\n{'R2 BARCODES':^70}")
    print("-" * 70)
    for bc in r2_barcodes[:10]:
        print(f"  {bc.barcode:<12} {bc.name}")
    if len(r2_barcodes) > 10:
        print(f"  ... and {len(r2_barcodes) - 10} more")

    # Generate FASTQ files
    print(f"\n{'GENERATING FASTQ FILES':^70}")
    print("-" * 70)

    r1_fastq = output_path / f"test_R1{ext}"
    r2_fastq = output_path / f"test_R2{ext}"

    start_time = time.time()

    if gzip_output:
        r1_f = gzip.open(r1_fastq, 'wt', compresslevel=6)
        r2_f = gzip.open(r2_fastq, 'wt', compresslevel=6)
    else:
        r1_f = open(r1_fastq, 'w')
        r2_f = open(r2_fastq, 'w')

    try:
        for i in range(num_reads):
            r1_bc_config = random.choice(r1_barcodes)
            r2_bc_config = random.choice(r2_barcodes)

            r1_bc = r1_bc_config.barcode
            r2_bc = r2_bc_config.barcode

            if random.random() < mismatch_rate and r1_mismatches > 0:
                r1_bc = introduce_mismatches(r1_bc, random.randint(1, r1_mismatches))

            if random.random() < mismatch_rate and r2_mismatches > 0:
                r2_bc = introduce_mismatches(r2_bc, random.randint(1, r2_mismatches))

            seq1 = r1_bc + random_barcode(94)
            seq2 = r2_bc + random_barcode(94)
            qual = 'I' * 100

            r1_f.write(f"@READ_{i}/1\n{seq1}\n+\n{qual}\n")
            r2_f.write(f"@READ_{i}/2\n{seq2}\n+\n{qual}\n")

            if (i + 1) % 100_000 == 0:
                elapsed = time.time() - start_time
                rate = (i + 1) / elapsed
                print(f"  Written {i+1:,} reads ({rate:,.0f}/sec)")
    finally:
        r1_f.close()
        r2_f.close()

    elapsed = time.time() - start_time
    print(f"\n  Completed in {elapsed:.1f}s")
    print(f"  Files: {r1_fastq}, {r2_fastq}")

    return r1_barcodes, r2_barcodes, str(r1_fastq), str(r2_fastq)


# =============================================================================
# BARCODE FILE LOADERS
# =============================================================================


def load_barcodes_from_tsv(filepath: str) -> List[BarcodeConfig]:
    """Load barcodes from a TSV file with intelligent column detection.

    Column resolution order:
      1. Match header names: barcode col → bc/barcode/sequence/index,
                             name col   → name/target/crf/sample/id
      2. Fallback: scan each column for uniform-length pure-ACGT content
         and use that as the barcode column; first non-barcode column as name.
    """
    BARCODE_COL_NAMES = {"bc", "barcode", "sequence", "index"}
    NAME_COL_NAMES    = {"name", "target", "crf", "sample", "id"}
    ACGT_PAT          = re.compile(r"^[ACGT]+$", re.IGNORECASE)

    with open(filepath, 'r', encoding='utf-8') as f:
        raw_lines = [ln.rstrip('\n') for ln in f if ln.strip()]

    if not raw_lines:
        raise ValueError(f"Empty barcode file: {filepath}")

    header_parts = raw_lines[0].split('\t')
    col_names    = [c.lower().strip() for c in header_parts]
    data_rows    = [ln.split('\t') for ln in raw_lines[1:] if ln.strip()]

    if not data_rows:
        raise ValueError(f"No data rows in barcode file: {filepath}")

    # --- Pass 1: match by column name ---
    barcode_idx, name_idx = None, None
    for i, col in enumerate(col_names):
        if barcode_idx is None and col in BARCODE_COL_NAMES:
            barcode_idx = i
        if name_idx is None and col in NAME_COL_NAMES:
            name_idx = i
        if barcode_idx is not None and name_idx is not None:
            break

    # --- Pass 2: content-based fallback ---
    if barcode_idx is None or name_idx is None:
        print(f"Warning: Could not find required columns in {filepath}. "
              f"Attempting content-based detection.")
        for i in range(len(col_names)):
            values = [row[i].strip() for row in data_rows if len(row) > i and row[i].strip()]
            if not values:
                continue
            is_seq   = all(ACGT_PAT.fullmatch(v) for v in values)
            lengths  = [len(v) for v in values]
            mode_len = max(set(lengths), key=lengths.count) if lengths else 0
            equal_len = mode_len > 0 and all(l == mode_len for l in lengths)

            if is_seq and equal_len and barcode_idx is None:
                barcode_idx = i
            elif name_idx is None:
                name_idx = i

            if barcode_idx is not None and name_idx is not None:
                break

    if barcode_idx is None:
        raise ValueError(f"Could not identify barcode column in {filepath}.")

    barcodes = []
    for row in data_rows:
        if len(row) <= barcode_idx:
            continue
        barcode = row[barcode_idx].strip()
        if not barcode:
            continue
        name = (row[name_idx].strip()
                if name_idx is not None and len(row) > name_idx
                else barcode)
        barcodes.append(BarcodeConfig(barcode=barcode, name=name))

    return barcodes


def load_barcodes_from_fasta(filepath: str) -> List[BarcodeConfig]:
    """Load barcodes from a FASTA file.

    Record name (text after '>') → BarcodeConfig.name
    Sequence lines               → BarcodeConfig.barcode (uppercased, multi-line joined)
    """
    barcodes = []
    current_name: Optional[str] = None
    seq_parts: List[str] = []

    with open_file(filepath, 'r') as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_name is not None and seq_parts:
                    barcodes.append(BarcodeConfig(
                        barcode=''.join(seq_parts).upper(),
                        name=current_name
                    ))
                current_name = line[1:].strip()
                seq_parts = []
            else:
                seq_parts.append(line)

    # flush last record
    if current_name is not None and seq_parts:
        barcodes.append(BarcodeConfig(
            barcode=''.join(seq_parts).upper(),
            name=current_name
        ))

    if not barcodes:
        raise ValueError(f"No valid FASTA records found in {filepath}.")

    return barcodes


def load_barcodes_auto(filepath: str, json_key: str = 'r1_barcodes') -> List[BarcodeConfig]:
    """Detect barcode file format (FASTA / TSV / JSON) and load accordingly.

    Detection order:
      1. Extension: .json              → JSON  (use json_key to pick r1/r2 list)
                    .fa / .fasta / .fa.gz / .fasta.gz → FASTA
                    .tsv / .txt / .csv → TSV
      2. Content peek (ambiguous ext): '{' → JSON, '>' → FASTA, else TSV

    Args:
        filepath:  Path to the barcode file.
        json_key:  Key to extract from a JSON file — 'r1_barcodes' or 'r2_barcodes'.
                   Ignored for FASTA/TSV.
    """
    lower = str(filepath).lower()
    fasta_exts = ('.fa', '.fasta', '.fa.gz', '.fasta.gz')
    tsv_exts   = ('.tsv', '.txt', '.csv')

    if lower.endswith('.json'):
        return _load_barcodes_from_json_key(filepath, json_key)
    if any(lower.endswith(ext) for ext in fasta_exts):
        return load_barcodes_from_fasta(filepath)
    if any(lower.endswith(ext) for ext in tsv_exts):
        return load_barcodes_from_tsv(filepath)

    # Ambiguous extension — peek at first non-whitespace character
    with open_file(filepath, 'r') as f:
        for line in f:
            stripped = line.strip()
            if stripped:
                if stripped.startswith('{'):
                    return _load_barcodes_from_json_key(filepath, json_key)
                if stripped.startswith('>'):
                    return load_barcodes_from_fasta(filepath)
                return load_barcodes_from_tsv(filepath)

    raise ValueError(f"Cannot determine barcode file format for {filepath}.")


def _load_barcodes_from_json_key(filepath: str, key: str) -> List[BarcodeConfig]:
    """Extract one list (r1_barcodes or r2_barcodes) from a JSON barcode file."""
    with open(filepath, 'r') as f:
        data = json.load(f)
    if key not in data:
        raise KeyError(f"Key '{key}' not found in {filepath}. Available keys: {list(data.keys())}")
    return [BarcodeConfig(barcode=bc['barcode'], name=bc['name']) for bc in data[key]]


# =============================================================================
# PUBLIC API
# =============================================================================

def run_demux(
    r1_path: str,
    r2_path: str,
    r1_barcode_file: str,
    output_dir: str,
    r2_barcode_file: str = None,
    error_rate: float = None,
    r1_mismatches: int = 1,
    r2_mismatches: int = 1,
    streaming: bool = False,
    trim_barcode: bool = True,
    gzip_output: bool = True,
    compression_level: int = 6,
    workers: int = None,
    chunk_size: int = 100_000,
    buffer_size: int = 10_000,
    max_open_files: int = 200,
    memory_target: float = None,
) -> Dict:
    """Demultiplex paired-end FASTQ files from raw barcode files.

    Importable entry point for external callers (e.g. pipeline wrappers).
    Handles barcode loading internally.

    Args:
        r1_path:           Path to R1 FASTQ (plain or gzipped).
        r2_path:           Path to R2 FASTQ (plain or gzipped).
        r1_barcode_file:   Barcode file for R1 — FASTA, TSV, JSON; auto-detected.
        output_dir:        Directory for output files.
        r2_barcode_file:   Barcode file for R2. Defaults to r1_barcode_file when
                           the same barcodes are used on both ends.
        error_rate:        Fractional mismatch rate (e.g. 0.1). When set, overrides
                           r1_mismatches/r2_mismatches via round(error_rate × barcode_len).
        r1_mismatches:     Max substitutions in R1 barcode (default 1). Ignored
                           when error_rate is provided.
        r2_mismatches:     Max substitutions in R2 barcode (default 1). Ignored
                           when error_rate is provided.
        streaming:         Use StreamingDemultiplexer — low memory, recommended
                           for files > 1 GB (default False).
        trim_barcode:      Remove barcode prefix from reads (default True).
        gzip_output:       Write .fastq.gz output (default True).
        compression_level: gzip compression level 1–9 (default 6).
        workers:           Parallel workers for standard mode (None = all CPUs).
                           Unused in streaming mode.
        chunk_size:        Reads per chunk in standard mode (default 100 000).
                           Unused in streaming mode.
        buffer_size:       Records buffered per barcode pair in streaming mode
                           (default 10 000). Overridden when memory_target is set.
        max_open_files:    Max simultaneous output file handles in streaming mode
                           (default 200). Raise if barcode combinations > 100.
        memory_target:     Target RAM in GB; auto-calculates buffer_size when set
                           (streaming mode only). Recommended: available RAM × 0.7.

    Returns:
        dict mapping pair names to read counts (plus 'total' and 'unmatched').
    """
    for label, path in (("R1 FASTQ", r1_path), ("R2 FASTQ", r2_path),
                        ("R1 barcode file", r1_barcode_file)):
        if not os.path.exists(path):
            raise FileNotFoundError(f"{label} not found: {path}")

    r1_barcodes = load_barcodes_auto(r1_barcode_file, json_key='r1_barcodes')
    r2_barcodes = load_barcodes_auto(r2_barcode_file if r2_barcode_file else r1_barcode_file,
                                     json_key='r2_barcodes')

    if error_rate is not None:
        barcode_len = len(r1_barcodes[0].barcode)
        mm = round(error_rate * barcode_len)
        print(f"  error_rate={error_rate} × barcode_len={barcode_len} → mismatches={mm}")
        r1_mismatches = r2_mismatches = mm

    if streaming and memory_target is not None:
        buffer_size = calculate_buffer_size_for_memory(
            memory_target, len(r1_barcodes), len(r2_barcodes)
        )
        print(f"  Memory target: {memory_target} GB → buffer_size = {buffer_size:,}")

    if streaming:
        demuxer = StreamingDemultiplexer(
            r1_barcodes=r1_barcodes,
            r2_barcodes=r2_barcodes,
            output_dir=output_dir,
            trim_barcode=trim_barcode,
            r1_mismatches=r1_mismatches,
            r2_mismatches=r2_mismatches,
            gzip_output=gzip_output,
            compression_level=compression_level,
            buffer_size=buffer_size,
            max_open_files=max_open_files,
        )
    else:
        demuxer = FastDemultiplexer(
            r1_barcodes=r1_barcodes,
            r2_barcodes=r2_barcodes,
            output_dir=output_dir,
            num_workers=workers,
            trim_barcode=trim_barcode,
            chunk_size=chunk_size,
            r1_mismatches=r1_mismatches,
            r2_mismatches=r2_mismatches,
            gzip_output=gzip_output,
            compression_level=compression_level,
        )

    return demuxer.demultiplex(r1_path, r2_path)


# =============================================================================
# MAIN FUNCTION
# =============================================================================

HELP_DESCRIPTION = """
================================================================================
                        FAST DNA BARCODE DEMULTIPLEXER
================================================================================

A high-performance tool for demultiplexing paired-end FASTQ files.

MODES:
  Standard: Fast, in-memory processing (~500K-1M reads/sec)
  Streaming: Low memory, handles large files (~200K-400K reads/sec)

================================================================================
"""

def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        prog='demux.py',
        description=HELP_DESCRIPTION,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('--version', '-V', action='version', version=f'%(prog)s {__version__}')

    # Input/Output
    io_group = parser.add_argument_group('Input/Output')
    io_group.add_argument('--r1', metavar='FILE', help='R1 FASTQ file')
    io_group.add_argument('--r2', metavar='FILE', help='R2 FASTQ file')
    io_group.add_argument('--output', '-o', metavar='DIR', default='demultiplexed', help='Output directory')
    io_group.add_argument('--gzip-output', '-z', action='store_true', help='Compress output')
    io_group.add_argument('--compression-level', metavar='N', type=int, default=6, choices=range(1, 10))

    # Barcode files
    barcode_group = parser.add_argument_group('Barcode Files')
    barcode_group.add_argument('--r1-barcodes', metavar='FILE', help='R1 barcodes (FASTA / TSV / JSON)')
    barcode_group.add_argument('--r2-barcodes', metavar='FILE', help='R2 barcodes (FASTA / TSV / JSON); defaults to --r1-barcodes')

    # Mismatch settings
    mismatch_group = parser.add_argument_group('Mismatch Settings')
    mismatch_group.add_argument('--r1-mismatches', metavar='N', type=int, default=1, help='R1 max mismatches (default: 1)')
    mismatch_group.add_argument('--r2-mismatches', metavar='N', type=int, default=1, help='R2 max mismatches (default: 1)')

    # Processing mode
    mode_group = parser.add_argument_group('Processing Mode')
    mode_group.add_argument('--streaming', '-s', action='store_true', help='Use streaming mode')
    mode_group.add_argument('--buffer-size', metavar='N', type=int, default=10_000, help='Buffer size')
    mode_group.add_argument('--max-open-files', metavar='N', type=int, default=200, help='Max open files')
    mode_group.add_argument('--memory-target', metavar='GB', type=float, default=None, help='Target memory (GB)')
    mode_group.add_argument('--workers', '-w', metavar='N', type=int, default=None, help='Parallel workers')
    mode_group.add_argument('--chunk-size', metavar='N', type=int, default=100_000, help='Chunk size')

    # Benchmark
    bench_group = parser.add_argument_group('Benchmark Mode')
    bench_group.add_argument('--benchmark', action='store_true', help='Run benchmark')
    bench_group.add_argument('--reads', metavar='N', type=int, default=1_000_000, help='Test reads')
    bench_group.add_argument('--r1-barcode-count', metavar='N', type=int, default=100, help='R1 barcode count')
    bench_group.add_argument('--r2-barcode-count', metavar='N', type=int, default=100, help='R2 barcode count')
    bench_group.add_argument('--mismatch-rate', metavar='RATE', type=float, default=0.1, help='Mismatch rate')
    bench_group.add_argument('--gzip-test-data', action='store_true', help='Gzip test files')

    # Other
    other_group = parser.add_argument_group('Other Options')
    other_group.add_argument('--no-trim', action='store_true', help='Do not trim barcodes')
    other_group.add_argument('--quiet', '-q', action='store_true', help='Suppress messages')

    args = parser.parse_args()

    # Print banner
    if not args.quiet:
        mode_str = "Streaming" if args.streaming else "Standard"
        print(f"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                      FAST DNA BARCODE DEMULTIPLEXER                          ║
║                              Version {__version__}                                     ║
║  Mode: {mode_str:<72}║
╚══════════════════════════════════════════════════════════════════════════════╝
        """)

    _common = dict(
        streaming=args.streaming,
        trim_barcode=not args.no_trim,
        gzip_output=args.gzip_output,
        compression_level=args.compression_level,
        workers=args.workers,
        chunk_size=args.chunk_size,
        buffer_size=args.buffer_size,
        max_open_files=args.max_open_files,
        memory_target=args.memory_target,
    )

    def _run_loaded(r1_path, r2_path, r1_barcodes, r2_barcodes, r1_mm, r2_mm):
        """Create demuxer from pre-loaded barcodes and run — for benchmark/JSON paths."""
        buf = _common['buffer_size']
        if _common['streaming'] and _common['memory_target'] is not None:
            buf = calculate_buffer_size_for_memory(
                _common['memory_target'], len(r1_barcodes), len(r2_barcodes)
            )
            print(f"  Memory target: {_common['memory_target']} GB → buffer_size = {buf:,}")
        if _common['streaming']:
            demuxer = StreamingDemultiplexer(
                r1_barcodes=r1_barcodes, r2_barcodes=r2_barcodes,
                output_dir=args.output, trim_barcode=_common['trim_barcode'],
                r1_mismatches=r1_mm, r2_mismatches=r2_mm,
                gzip_output=_common['gzip_output'],
                compression_level=_common['compression_level'],
                buffer_size=buf, max_open_files=_common['max_open_files'],
            )
        else:
            demuxer = FastDemultiplexer(
                r1_barcodes=r1_barcodes, r2_barcodes=r2_barcodes,
                output_dir=args.output, num_workers=_common['workers'],
                trim_barcode=_common['trim_barcode'], chunk_size=_common['chunk_size'],
                r1_mismatches=r1_mm, r2_mismatches=r2_mm,
                gzip_output=_common['gzip_output'],
                compression_level=_common['compression_level'],
            )
        return demuxer.demultiplex(r1_path, r2_path)

    if args.benchmark:
        r1_barcodes, r2_barcodes, r1_fastq, r2_fastq = create_test_data(
            num_reads=args.reads,
            num_r1_barcodes=args.r1_barcode_count,
            num_r2_barcodes=args.r2_barcode_count,
            r1_mismatches=args.r1_mismatches,
            r2_mismatches=args.r2_mismatches,
            mismatch_rate=args.mismatch_rate,
            output_dir=".",
            gzip_output=args.gzip_test_data
        )
        _run_loaded(r1_fastq, r2_fastq, r1_barcodes, r2_barcodes,
                    args.r1_mismatches, args.r2_mismatches)

    elif args.r1 and args.r2:
        if not os.path.exists(args.r1):
            print(f"ERROR: R1 file not found: {args.r1}")
            sys.exit(1)
        if not os.path.exists(args.r2):
            print(f"ERROR: R2 file not found: {args.r2}")
            sys.exit(1)

        if args.r1_barcodes:
            run_demux(
                r1_path=args.r1, r2_path=args.r2,
                r1_barcode_file=args.r1_barcodes,
                r2_barcode_file=args.r2_barcodes,
                output_dir=args.output,
                r1_mismatches=args.r1_mismatches,
                r2_mismatches=args.r2_mismatches,
                **_common,
            )
        else:
            print("ERROR: Please provide barcode files via --r1-barcodes (FASTA / TSV / JSON).")
            sys.exit(1)

    else:
        parser.print_help()
        print("\nQuick start: python demux.py --benchmark")


if __name__ == "__main__":
    main()
