from Bio import SeqIO, pairwise2
from Bio.pairwise2 import format_alignment
import argparse

def extract_region(fasta_path, start, end):
    for record in SeqIO.parse(fasta_path, "fasta"):
        return record.seq[start-1:end], record.id  # 1-based indexing fix
    raise ValueError("FASTA is empty or not found.")

def main():
    parser = argparse.ArgumentParser(description="Visualize variation in target region.")
    parser.add_argument("--ref", required=True, help="Reference genome FASTA")
    parser.add_argument("--qry", required=True, help="Query genome FASTA")
    parser.add_argument("--start", type=int, required=True, help="Start position in reference (1-based)")
    parser.add_argument("--end", type=int, required=True, help="End position in reference (1-based)")
    parser.add_argument("--flank", type=int, default=50, help="Number of flanking bases")
    parser.add_argument("--outfile", help="Save output to file")

    args = parser.parse_args()

    # Define extraction region with flanks
    ref_start = max(1, args.start - args.flank)
    ref_end = args.end + args.flank

    # Extract sequences
    ref_seq, ref_id = extract_region(args.ref, ref_start, ref_end)
    qry_seq, qry_id = extract_region(args.qry, ref_start, ref_end)  # Assuming same coordinates

    # Align sequences
    alignment = pairwise2.align.globalms(ref_seq, qry_seq, 2, -1, -5, -0.5, one_alignment_only=True)[0]
    ref_aligned, qry_aligned, score, begin, end = alignment

    # Format alignment with mismatches marked
    match_line = "".join("|" if r == q else "*" if r != "-" and q != "-" else " " for r, q in zip(ref_aligned, qry_aligned))

    # Chunked output for readability
    output_lines = []
    chunk_size = 80
    output_lines.append(f"\n--- Alignment of region REF:{args.start}-{args.end} ---")
    for i in range(0, len(ref_aligned), chunk_size):
        r_chunk = ref_aligned[i:i+chunk_size]
        m_chunk = match_line[i:i+chunk_size]
        q_chunk = qry_aligned[i:i+chunk_size]
        output_lines.append(f"REF ({ref_id}): {r_chunk}")
        output_lines.append(f"               {m_chunk}")
        output_lines.append(f"QRY ({qry_id}): {q_chunk}")
        output_lines.append("")

    output = "\n".join(output_lines)
    print(output)

    if args.outfile:
        with open(args.outfile, "w") as f:
            f.write(output)

if __name__ == "__main__":
    main()

