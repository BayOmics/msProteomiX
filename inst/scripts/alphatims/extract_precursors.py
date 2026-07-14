#!/usr/bin/env python3
"""
extract_precursors.py — Extract MS1 precursors from timsTOF .d/.hdf files
and compute Contamination Ratio (CR) per RT slice.

Reproduces the analysis from:
  "Multimodal single cell-resolved spatial proteomics reveal
   pancreatic tumor heterogeneity" (Fig. 2f, 2g)

Strategy:
  1. Load timsTOF data via alphatims
  2. Select n_slices equally-spaced MS1 frames across the LC gradient
  3. For each frame, extract all MS1 precursors and compute CR
  4. Export summary CSV + optional per-slice CSVs (for heatmap)

Usage:
  # Single file
  python extract_precursors.py sample.d --output_dir ./cr_output

  # Single HDF file (Mac-friendly, after converting on Windows)
  python extract_precursors.py sample.hdf --output_dir ./cr_output

  # Batch: all .d folders in a directory
  python extract_precursors.py /path/to/data/ --output_dir ./cr_output

  # Export per-slice CSVs for heatmap at RT=32 min
  python extract_precursors.py sample.d --export_slices --heatmap_rts 32,48

Dependencies:
  pip install alphatims

Author: msProteomiX project (BayOmics)
"""

import argparse
import logging
import os
import platform
import sys
import time

import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger(__name__)


# ==============================================================================
# Core functions
# ==============================================================================

def parse_line_points(line_points_str):
    """Parse dividing line from 'mz1,im1,mz2,im2' string.

    Returns (k, b) where: IM = k * mz + b
    """
    pts = [float(x) for x in line_points_str.split(",")]
    if len(pts) != 4:
        raise ValueError(
            f"--line_points requires exactly 4 values (mz1,im1,mz2,im2), "
            f"got {len(pts)}"
        )
    mz1, im1, mz2, im2 = pts
    k = (im2 - im1) / (mz2 - mz1)
    b = im1 - k * mz1
    return k, b


def compute_cr(mz_values, mobility_values, intensity_values, k, b):
    """Compute Contamination Ratio for a set of precursors.

    Convention (cf. Meier et al. diaPASEF, Nat Methods 2020):
      - Points WHERE mobility >= k*mz + b are singly-charged contaminants
        (higher 1/K0 at same m/z due to CCS/z with z=1).
      - Points WHERE mobility <  k*mz + b are multi-charge peptide signal.

    CR = sum(contaminant_intensity) / sum(signal_intensity) * 100
    """
    line_values = k * mz_values + b
    is_contam = mobility_values >= line_values

    contam_inten = float(intensity_values[is_contam].sum())
    signal_inten = float(intensity_values[~is_contam].sum())

    if signal_inten == 0:
        return np.nan
    return (contam_inten / signal_inten) * 100.0


def select_slice_frames(data, n_slices=10):
    """Select n_slices equally-spaced MS1 frames across the LC gradient.

    Strategy (matching the paper):
      1. Get total RT range from MS1 frames only (MsMsType == 0)
      2. Generate n_slices equally-spaced target RT values via linspace
      3. For each target RT, find the closest MS1 frame

    Parameters
    ----------
    data : alphatims.bruker.TimsTOF
        Loaded TimsTOF dataset.
    n_slices : int
        Number of slices to select (default: 10).

    Returns
    -------
    list of (frame_index, rt_min)
        Selected frame indices and their RT in minutes.
    """
    ms1_mask = data.frames["MsMsType"] == 0
    ms1_frame_indices = data.frames.index[ms1_mask].values
    ms1_rt_sec = data.frames.loc[ms1_mask, "Time"].values
    ms1_rt_min = ms1_rt_sec / 60.0

    target_rts = np.linspace(ms1_rt_min.min(), ms1_rt_min.max(), n_slices)

    selected = []
    for target_rt in target_rts:
        closest_idx = np.argmin(np.abs(ms1_rt_min - target_rt))
        frame_idx = int(ms1_frame_indices[closest_idx])
        actual_rt = float(ms1_rt_min[closest_idx])
        selected.append((frame_idx, actual_rt))

    return selected


def extract_frame_data(data, frame_idx):
    """Extract all MS1 precursors from a single frame.

    Uses alphatims 5D slicing: data[frame_idx, :, 0]
      - frame_idx: specific frame
      - ':': all mobility scans
      - 0: precursor_index=0 (MS1 only)

    Parameters
    ----------
    data : alphatims.bruker.TimsTOF
        Loaded TimsTOF dataset.
    frame_idx : int
        Frame index to extract.

    Returns
    -------
    pd.DataFrame
        Columns: rt_values_min, mz_values, mobility_values,
                 corrected_intensity_values
    """
    df = data[frame_idx, :, 0]
    cols = ["rt_values_min", "mz_values", "mobility_values",
            "corrected_intensity_values"]
    # Ensure all required columns exist
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise KeyError(
            f"Missing columns in alphatims output: {missing}. "
            f"Available: {list(df.columns)}"
        )
    return df[cols].copy()


def process_single_file(input_path, output_dir, n_slices, k, b,
                         export_slices, heatmap_rts):
    """Process a single .d or .hdf file.

    Parameters
    ----------
    input_path : str
        Path to .d folder or .hdf file.
    output_dir : str
        Output directory.
    n_slices : int
        Number of RT slices.
    k, b : float
        CR dividing line parameters (IM = k * mz + b).
    export_slices : bool
        If True, export per-slice CSV files.
    heatmap_rts : list of float or None
        If provided, also export slices closest to these RT values.

    Returns
    -------
    pd.DataFrame
        Summary dataframe with CR per slice.
    """
    import alphatims.bruker

    # macOS warning
    if platform.system() == "Darwin" and input_path.endswith(".d"):
        logger.warning(
            "macOS detected: mz_values and mobility_values from .d files "
            "are ESTIMATED (errors up to 6 Th possible). "
            "For accurate values, convert to .hdf on Windows first: "
            "alphatims export hdf sample.d"
        )

    logger.info(f"Loading: {input_path}")
    t0 = time.time()
    data = alphatims.bruker.TimsTOF(input_path)
    logger.info(
        f"Loaded in {time.time() - t0:.1f}s: "
        f"{len(data)} detector events, "
        f"{data.frame_max_index} frames, "
        f"RT range: {data.rt_values[0]/60:.1f}-{data.rt_values[-1]/60:.1f} min"
    )

    sample_name = data.sample_name

    # Select frames
    selected_frames = select_slice_frames(data, n_slices)
    logger.info(
        f"Selected {len(selected_frames)} slices: "
        f"RT = {[f'{rt:.1f}' for _, rt in selected_frames]} min"
    )

    # Determine which slices to export as CSV
    export_rt_set = set()
    if heatmap_rts:
        for target_rt in heatmap_rts:
            # Find closest selected frame
            closest = min(selected_frames, key=lambda x: abs(x[1] - target_rt))
            export_rt_set.add(closest[1])
            logger.info(
                f"Heatmap RT={target_rt:.1f} -> closest slice RT={closest[1]:.1f}"
            )

    # Create output directories
    sample_dir = os.path.join(output_dir, sample_name)
    if export_slices or export_rt_set:
        os.makedirs(sample_dir, exist_ok=True)

    # Process each slice
    summary_rows = []
    for i, (frame_idx, rt_min) in enumerate(selected_frames):
        logger.info(
            f"  Slice {i+1}/{len(selected_frames)}: "
            f"frame={frame_idx}, RT={rt_min:.2f} min"
        )

        try:
            df_slice = extract_frame_data(data, frame_idx)
        except Exception as e:
            logger.error(f"  Failed to extract frame {frame_idx}: {e}")
            summary_rows.append({
                "file_name": f"{sample_name}_slice_rt{rt_min:.1f}",
                "rt_min": rt_min,
                "tot_inten": np.nan,
                "contamination_ratio": np.nan,
            })
            continue

        n_precursors = len(df_slice)
        if n_precursors == 0:
            logger.warning(f"  No precursors in frame {frame_idx}")
            summary_rows.append({
                "file_name": f"{sample_name}_slice_rt{rt_min:.1f}",
                "rt_min": rt_min,
                "tot_inten": 0,
                "contamination_ratio": np.nan,
            })
            continue

        # Compute CR
        mz = df_slice["mz_values"].values
        mob = df_slice["mobility_values"].values
        inten = df_slice["corrected_intensity_values"].values.astype(np.float64)

        cr = compute_cr(mz, mob, inten, k, b)
        tot_inten = float(np.log2(inten.sum())) if inten.sum() > 0 else 0

        slice_name = f"{sample_name}_slice_rt{rt_min:.1f}"
        summary_rows.append({
            "file_name": slice_name,
            "rt_min": rt_min,
            "tot_inten": tot_inten,
            "contamination_ratio": cr,
        })

        logger.info(
            f"    {n_precursors} precursors, CR={cr:.1f}%, "
            f"TIC(log2)={tot_inten:.1f}"
        )

        # Export per-slice CSV if requested
        should_export = export_slices or (rt_min in export_rt_set)
        if should_export:
            slice_csv = os.path.join(sample_dir, f"{slice_name}.csv")
            df_slice.to_csv(slice_csv, index=False)
            logger.info(f"    Exported: {slice_csv}")

    # Write summary CSV
    summary_df = pd.DataFrame(summary_rows)
    summary_csv = os.path.join(output_dir, f"{sample_name}_cr_summary.csv")
    summary_df.to_csv(summary_csv, index=False)
    logger.info(f"Summary saved: {summary_csv}")

    return summary_df


def find_input_files(input_path):
    """Find all .d folders and .hdf files in the input path.

    Parameters
    ----------
    input_path : str
        A .d folder, .hdf file, or directory containing them.

    Returns
    -------
    list of str
        Paths to process.
    """
    if input_path.endswith(".d") or input_path.endswith(".hdf"):
        if os.path.exists(input_path):
            return [input_path]
        else:
            logger.error(f"File not found: {input_path}")
            return []

    # Directory: scan for .d folders and .hdf files
    results = []
    if os.path.isdir(input_path):
        for item in sorted(os.listdir(input_path)):
            full_path = os.path.join(input_path, item)
            if item.endswith(".d") and os.path.isdir(full_path):
                results.append(full_path)
            elif item.endswith(".hdf") and os.path.isfile(full_path):
                results.append(full_path)

    if not results:
        logger.error(
            f"No .d folders or .hdf files found in: {input_path}"
        )

    return results


# ==============================================================================
# Main
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Extract MS1 precursors from timsTOF .d/.hdf files "
            "and compute Contamination Ratio (CR) per RT slice."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Single file
  python extract_precursors.py sample.d

  # Batch processing
  python extract_precursors.py /path/to/data/

  # Export heatmap slices at RT=32 and RT=48
  python extract_precursors.py sample.d --export_slices --heatmap_rts 32,48

  # Custom dividing line and 20 slices
  python extract_precursors.py sample.d --n_slices 20 --line_points 400,0.85,900,1.25
        """,
    )
    parser.add_argument(
        "input",
        help="Path to .d folder, .hdf file, or directory containing them",
    )
    parser.add_argument(
        "--output_dir",
        default="./cr_output",
        help="Output directory (default: ./cr_output)",
    )
    parser.add_argument(
        "--n_slices",
        type=int,
        default=10,
        help="Number of equally-spaced RT slices per sample (default: 10)",
    )
    parser.add_argument(
        "--line_points",
        default="350,0.8,950,1.3",
        help=(
            "CR dividing line as mz1,im1,mz2,im2 "
            "(default: 350,0.8,950,1.3 per paper)"
        ),
    )
    parser.add_argument(
        "--export_slices",
        action="store_true",
        help="Export per-slice CSV files (needed for heatmap in R)",
    )
    parser.add_argument(
        "--heatmap_rts",
        default=None,
        help=(
            "Comma-separated RT values (min) for heatmap export. "
            "Implies --export_slices for those RTs. Example: 32,48"
        ),
    )

    args = parser.parse_args()

    # Parse parameters
    k, b = parse_line_points(args.line_points)
    logger.info(
        f"CR dividing line: IM = {k:.6f} * m/z + {b:.4f} "
        f"(from points {args.line_points})"
    )

    heatmap_rts = None
    if args.heatmap_rts:
        heatmap_rts = [float(x) for x in args.heatmap_rts.split(",")]
        logger.info(f"Heatmap export at RT: {heatmap_rts} min")

    os.makedirs(args.output_dir, exist_ok=True)

    # Find input files
    input_files = find_input_files(args.input)
    if not input_files:
        sys.exit(1)

    logger.info(f"Found {len(input_files)} file(s) to process")

    # Process each file
    all_summaries = []
    for i, input_path in enumerate(input_files):
        logger.info(f"\n{'='*60}")
        logger.info(f"Processing [{i+1}/{len(input_files)}]: {input_path}")
        logger.info(f"{'='*60}")

        try:
            summary = process_single_file(
                input_path=input_path,
                output_dir=args.output_dir,
                n_slices=args.n_slices,
                k=k,
                b=b,
                export_slices=args.export_slices,
                heatmap_rts=heatmap_rts,
            )
            all_summaries.append(summary)
        except Exception as e:
            logger.error(f"Failed to process {input_path}: {e}")
            import traceback
            traceback.print_exc()
            continue

    # Combined summary
    if len(all_summaries) > 1:
        combined = pd.concat(all_summaries, ignore_index=True)
        combined_csv = os.path.join(args.output_dir, "all_cr_summary.csv")
        combined.to_csv(combined_csv, index=False)
        logger.info(f"\nCombined summary: {combined_csv}")

    logger.info(f"\nDone! Output directory: {args.output_dir}")


if __name__ == "__main__":
    main()
