#!/usr/bin/env python3
import subprocess
import argparse
import os

def parse_args():
    parser = argparse.ArgumentParser(description="Perform MSA using MAFFT and output CLUSTAL format.")
    parser.add_argument("--infile", required=True, help="Input multi-FASTA file")
    parser.add_argument("--outdir", default="msa_output", help="Directory to store the alignment file")
    parser.add_argument("--outfile", default="aligned_output.aln", help="Output alignment filename")
    parser.add_argument("--mafft", default="mafft", help="Path to MAFFT executable")
    return parser.parse_args()

def run_mafft(input_fasta, output_file, mafft_path="mafft"):
    print(f"Running MAFFT (CLUSTAL format) on: {input_fasta}")
    result = subprocess.run(
        [mafft_path, "--auto", "--clustalout", "--adjustdirection", input_fasta],
        capture_output=True,
        text=True
    )

    if result.returncode == 0:
        with open(output_file, "w") as out_f:
            out_f.write(result.stdout)
        print(f"CLUSTAL alignment saved to: {output_file}")
    else:
        print("MAFFT failed:")
        print(result.stderr)

def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    output_path = os.path.join(args.outdir, args.outfile)
    run_mafft(args.infile, output_path, args.mafft)

if __name__ == "__main__":
    main()
