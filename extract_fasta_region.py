from Bio import SeqIO
import argparse

def extract_region(fasta_path, seq_id, start, end, output_path):
    """Extract a region from a given sequence ID and save as FASTA."""
    with open(fasta_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id == seq_id:
                sub_seq = record.seq[start-1:end]  # convert 1-based to 0-based
                sub_record = record[start-1:end]
                sub_record.id = f"{seq_id}_{start}_{end}"
                sub_record.description = f"{seq_id}:{start}-{end}"
                with open(output_path, "w") as out:
                    SeqIO.write(sub_record, out, "fasta")
                print(f"Extracted: {seq_id}:{start}-{end} to {output_path}")
                return
    raise ValueError(f"Sequence ID '{seq_id}' not found in {fasta_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract region from FASTA by ID and coordinates")
    parser.add_argument("--fasta", required=True, help="FASTA file")
    parser.add_argument("--seq_id", required=True, help="Sequence ID (e.g., NC_009950.1)")
    parser.add_argument("--start", type=int, required=True, help="Start position (1-based)")
    parser.add_argument("--end", type=int, required=True, help="End position (1-based)")
    parser.add_argument("--output", required=True, help="Output FASTA file")

    args = parser.parse_args()
    extract_region(args.fasta, args.seq_id, args.start, args.end, args.output)
