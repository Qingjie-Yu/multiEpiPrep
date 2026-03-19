import os
import pysam
from multiprocessing import Pool
from typing import List, Optional
import pandas as pd
import pyranges as pr
from .utils import download_bed, get_files_path, get_files_prefix, get_cpu_core


def load_reference_bed(bed_path) -> set:
    """Load reference BED as a set of (chrom, start, end) intervals for fast lookup,
    backed by PyRanges for overlap queries."""
    ref = pd.read_csv(
        bed_path,
        sep="\t",
        header=None,
        usecols=[0, 1, 2],
        names=["Chromosome", "Start", "End"]
    )
    return pr.PyRanges(ref)


def compute_overlap_ratio(bam_path: str, prefix: str, ref_gr: pr.PyRanges) -> dict:
    """Count reads in BAM and compute fraction overlapping reference regions.
    """
    rows = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            rows.append({
                "Chromosome": bam.get_reference_name(read.reference_id),
                "Start": read.reference_start,
                "End": read.reference_end,
                "ReadID": read.query_name
            })

    total_reads = len(rows)
    if total_reads == 0:
        return {"prefix": prefix, "total_reads": 0, "overlap_reads": 0, "overlap_ratio": 0.0}

    read_gr = pr.PyRanges(pd.DataFrame(rows))
    overlap_reads = len(read_gr.join(ref_gr).df["ReadID"].unique())

    return {
        "prefix": prefix,
        "total_reads": total_reads,
        "overlap_reads": overlap_reads,
        "overlap_ratio": overlap_reads / total_reads
    }


def qc(
    bam_dir: str,
    out_dir: str,
    ref_genome: str,
    source: str = "ccre",
    exclude: List[str] = ["unknown", "IgG_control"],
    min_reads: int = 1000,
    min_ratio: float = 0.0,
    threads: Optional[int] = None
):
    """Compute per-sample read overlap QC and filter by read count and overlap ratio.

    For each BAM in `bam_dir` (excluding specified prefixes), counts total mapped
    reads and the fraction overlapping a reference feature set (cCRE or DHS).
    Writes two TSVs: one with all samples, one with samples passing both filters.

    Parameters
    ----------
    bam_dir : str
        Directory containing input BAM files.
    out_dir : str
        Directory for output TSV files.
    ref_genome : str
        Reference genome; must be 'hg38' or 'mm10'.
    source : str
        Feature set to overlap against: 'ccre' or 'dhs'. Default 'ccre'.
    exclude : list of str
        Sample name tokens to exclude from analysis.
    min_reads : int
        Minimum total read count to pass QC. Default 1000.
    min_ratio : float
        Minimum overlap ratio to pass QC. Default 0.0 (no filter).
    threads : int, optional
        Worker processes for parallel BAM processing. Defaults to min(32, nCPU).
    """
    if ref_genome not in ["hg38", "mm10"]:
        raise ValueError("Only hg38 or mm10 supported for ref_genome")

    if source.lower() == "dhs":
        ref_gr = load_reference_bed(download_bed(f"DHS_{ref_genome}.bed"))
    elif source.lower() == "ccre":
        ref_gr = load_reference_bed(download_bed(f"ENCODE_cCRE_v4_{ref_genome}.bed"))
    else:
        raise ValueError(f"Unknown source '{source}': choose 'ccre' or 'dhs'")

    os.makedirs(out_dir, exist_ok=True)
    all_tsv = os.path.join(out_dir, "all_qc.tsv")
    filtered_tsv = os.path.join(out_dir, "filtered_qc.tsv")

    bam_list = get_files_path(bam_dir, ext=".bam")
    prefix_list = get_files_prefix(bam_list, ".bam")
    has_dash = ["-" in p for p in prefix_list]

    if all(has_dash):
        should_exclude = [any(tok in exclude for tok in p.split("-")) for p in prefix_list]
    elif not any(has_dash):
        should_exclude = [p in exclude for p in prefix_list]
    else:
        should_exclude = [any(item in p for item in exclude) for p in prefix_list]

    prefix_to_bam = {
        prefix: bam
        for prefix, bam, skip in zip(prefix_list, bam_list, should_exclude)
        if not skip
    }

    if threads is None:
        threads = min(32, get_cpu_core())

    with Pool(processes=threads) as pool:
        async_results = {
            prefix: pool.apply_async(compute_overlap_ratio, args=(bam, prefix, ref_gr))
            for prefix, bam in prefix_to_bam.items()
        }

        all_rows = []
        for prefix, ar in async_results.items():
            try:
                all_rows.append(ar.get())
            except Exception as e:
                print(f"Failed: {prefix} -> {e}")
                pool.terminate()
                return None

    all_rows.sort(key=lambda x: x["total_reads"], reverse=True)
    all_df = pd.DataFrame(all_rows, columns=["prefix", "total_reads", "overlap_reads", "overlap_ratio"])
    all_df.to_csv(all_tsv, sep="\t", index=False)

    filtered_df = all_df[
        (all_df["total_reads"] >= min_reads) &
        (all_df["overlap_ratio"] >= min_ratio)
    ]
    filtered_df.to_csv(filtered_tsv, sep="\t", index=False)

    print(f"Total samples: {len(all_df)}, passed QC: {len(filtered_df)}")
    return

HELP = """
Perform read-overlap-based quality control on BAM files.

Required:
  -i, --input       Directory containing input BAM files
  -o, --output      Output directory for QC results
  -g, --genome      Reference genome: hg38 or mm10

Optional:
  -s, --source      Reference feature set to overlap against: ccre or dhs (default: ccre)
  -r, --min-reads   Minimum total mapped reads to pass QC (default: 1000)
  -R, --min-ratio   Minimum overlap ratio to pass QC (default: 0.0, no filter)
  -e, --exclude     Comma-separated list of sample name tokens to exclude
                    (default: unknown,IgG_control). Exclusion matches substrings
                    in BAM filename prefixes. Pass an empty string to disable.
  -t, --threads     Number of parallel worker processes (default: min(32, nCPU))

Description:
  For each BAM file in the input directory, counts total mapped reads and computes
  the fraction overlapping a reference feature set (ENCODE cCRE v4 or DHS). Two
  TSV files are written to the output directory:
    all_qc.tsv      - QC stats for all non-excluded samples
    filtered_qc.tsv - Samples passing both --min-reads and --min-ratio thresholds

Example:
  multiEpiPrep qc -i BAM_DIR -o OUT_DIR -g hg38
  multiEpiPrep qc -i BAM_DIR -o OUT_DIR -g mm10 -s dhs -r 5000 -R 0.2 -e IgG_control
"""

def register_parser(parser):
    parser.add_argument("-i", "--input", required=True, metavar="BAM_DIR")
    parser.add_argument("-o", "--output", required=True, metavar="OUT_DIR")
    parser.add_argument("-g", "--genome", required=True, choices=["hg38", "mm10"], metavar="GENOME")
    parser.add_argument("-s", "--source", default="ccre", choices=["ccre", "dhs"], metavar="SOURCE")
    parser.add_argument("-r", "--min-reads", type=int, default=1000, metavar="N", dest="min_reads")
    parser.add_argument("-R", "--min-ratio", type=float, default=0.0, metavar="RATIO", dest="min_ratio")
    parser.add_argument("-e", "--exclude",
                        type=lambda s: s.split(",") if s else [],
                        default="unknown,IgG_control",
                        metavar="NAMES")
    parser.add_argument("-t", "--threads", type=int, default=None, metavar="N")
    parser.set_defaults(func=run)


def run(args):
    exclude = args.exclude if isinstance(args.exclude, list) else args.exclude.split(",")
    qc(
        bam_dir=args.input,
        out_dir=args.output,
        ref_genome=args.genome,
        source=args.source,
        exclude=exclude,
        min_reads=args.min_reads,
        min_ratio=args.min_ratio,
        threads=args.threads
    )