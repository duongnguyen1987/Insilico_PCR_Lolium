import os
import csv
import argparse
from collections import defaultdict
from Bio import SeqIO

def parse_arguments():
    parser = argparse.ArgumentParser(description="Extract sequences by PairID and annotate source.")
    parser.add_argument("--csv", required=True, help="Path to the primer_pairs.csv file")
    parser.add_argument("--fastas", nargs='+', required=True, help="List of input FASTA files")
    parser.add_argument("--outdir", default="grouped_fastas", help="Directory to store grouped FASTA files")
    return parser.parse_args()

def main():
    args = parse_arguments()
    os.makedirs(args.outdir, exist_ok=True)

    # Load PairIDs from CSV
    with open(args.csv, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        pair_ids = sorted(set(row["PairID"] for row in reader))

    # Collect sequences by PairID
    sequences_by_pairid = defaultdict(list)
    for fasta in args.fastas:
        source = os.path.basename(fasta)
        for record in SeqIO.parse(fasta, "fasta"):
            for pid in pair_ids:
                if pid in record.id:
                    record.id = f"{record.id}|{source}"
                    record.description = ""
                    sequences_by_pairid[pid].append(record)
                    break

    # Write grouped FASTA files
    for pid, records in sequences_by_pairid.items():
        output_fasta = os.path.join(args.outdir, f"{pid}_sequences.fasta")
        with open(output_fasta, "w") as out_f:
            SeqIO.write(records, out_f, "fasta")
        print(f"Wrote {len(records)} records for {pid} to {output_fasta}")

if __name__ == "__main__":
    main()
