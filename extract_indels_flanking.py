#!/usr/bin/env python3
from Bio import SeqIO, pairwise2
import argparse

def extract_region(fasta_path, start, end):
    """Extract a subsequence from a FASTA file using 0-based indexing."""
    for record in SeqIO.parse(fasta_path, "fasta"):
        return record.seq[start:end], record.id
    raise ValueError("Could not read sequence from FASTA.")

def main():
    parser = argparse.ArgumentParser(description="Compare a specific region in REF vs QRY")
    parser.add_argument("--ref", required=True, help="Reference FASTA file")
    parser.add_argument("--qry", required=True, help="Query FASTA file")
    parser.add_argument("--ref_start", type=int, required=True, help="Start (1-based) in REF")
    parser.add_argument("--ref_end", type=int, required=True, help="End (1-based) in REF")
    parser.add_argument("--qry_center", type=int, required=True, help="Center position (1-based) in QRY")
    parser.add_argument("--flank", type=int, default=50, help="Number of flanking bases")
    parser.add_argument("--outfile", help="Save output to file")

    args = parser.parse_args()

    # Convert to 0-based indexing for Python slicing
    ref_start = args.ref_start - args.flank - 1  # -1 to convert to 0-based
    ref_end = args.ref_end + args.flank
    qry_start = args.qry_center - args.flank - 1
    qry_end = args.qry_center + args.flank

    # Extract regions
    ref_seq, ref_id = extract_region(args.ref, ref_start, ref_end)
    qry_seq, qry_id = extract_region(args.qry, qry_start, qry_end)

    # Perform global alignment with high gap penalty
    alignment = pairwise2.align.globalms(
        ref_seq, qry_seq,
        2,    # match score
        -1,   # mismatch penalty
        -20,  # gap open penalty (high to discourage gaps)
        -1,   # gap extension penalty
        one_alignment_only=True
    )[0]

    ref_aln, qry_aln, score, begin, end = alignment

    # Generate match line: "|" = match, "*" = mismatch, " " = gap
    match_line = "".join(
        "|" if r == q else "*" if r != "-" and q != "-" else " "
        for r, q in zip(ref_aln, qry_aln)
    )

    # Determine dynamic label width for clean formatting
    ref_label_base = f"REF ({ref_id}): "
    qry_label_base = f"QRY ({qry_id}): "
    label_width = max(len(ref_label_base), len(qry_label_base))
    ref_label_fmt = f"{{:<{label_width}}}"

    # Calculate actual coordinate ranges for the header
    actual_ref_start = args.ref_start - args.flank
    actual_ref_end = args.ref_end + args.flank
    actual_qry_start = args.qry_center - args.flank
    actual_qry_end = args.qry_center + args.flank

    # Assemble output
    lines = []
    chunk = 80
    lines.append(
        f"\n--- Alignment: REF({actual_ref_start}-{actual_ref_end}) vs "
        f"QRY({actual_qry_start}-{actual_qry_end}, center={args.qry_center}) ---"
    )

    for i in range(0, len(ref_aln), chunk):
        ref_chunk = ref_aln[i:i+chunk]
        match_chunk = match_line[i:i+chunk]
        qry_chunk = qry_aln[i:i+chunk]

        lines.append(ref_label_fmt.format(ref_label_base) + ref_chunk)
        lines.append(" " * label_width + match_chunk)
        lines.append(ref_label_fmt.format(qry_label_base) + qry_chunk)
        lines.append("")

    result = "\n".join(lines)
    print(result)

    # Save to file if requested
    if args.outfile:
        with open(args.outfile, "w") as f:
            f.write(result)

if __name__ == "__main__":
    main()

