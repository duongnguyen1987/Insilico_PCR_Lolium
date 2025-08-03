import argparse
import csv
from Bio import SeqIO
from Bio.Seq import Seq

def find_amplicons_on_strand(seq_id, seq, f_primer, r_primer_rc, min_len, max_len, strand, seq_len):
    results = []
    start = 0
    while True:
        f_index = seq.find(f_primer, start)
        if f_index == -1:
            break
        search_start = f_index + len(f_primer)
        r_index = seq.find(r_primer_rc, search_start)
        if r_index != -1:
            product_len = r_index + len(r_primer_rc) - f_index
            if min_len <= product_len <= max_len:
                amp_start = f_index + 1
                amp_end = r_index + len(r_primer_rc)
                amplicon_seq = seq[f_index:amp_end]

                if strand == "-":
                    orig_amp_start = seq_len - amp_end + 1
                    orig_amp_end   = seq_len - amp_start + 1
                    fwd_start = seq_len - (f_index + len(f_primer)) + 1
                    fwd_end   = seq_len - f_index
                    rev_start = seq_len - (r_index + len(r_primer_rc)) + 1
                    rev_end   = seq_len - r_index
                else:
                    orig_amp_start = amp_start
                    orig_amp_end   = amp_end
                    fwd_start = amp_start
                    fwd_end   = fwd_start + len(f_primer) - 1
                    rev_start = r_index + 1
                    rev_end   = rev_start + len(r_primer_rc) - 1

                results.append({
                    "seq_id": seq_id,
                    "start": orig_amp_start,
                    "end": orig_amp_end,
                    "length": product_len,
                    "strand": strand,
                    "amplicon_seq": amplicon_seq,
                    "fwd_start": fwd_start,
                    "fwd_end": fwd_end,
                    "rev_start": rev_start,
                    "rev_end": rev_end,
                    "forward_primer": f_primer,
                    "reverse_primer_rc": r_primer_rc
                })
        start = f_index + 1
    return results

def find_amplicons_both_strands(record, f_primer, r_primer_rc, min_len, max_len):
    seq_fwd = str(record.seq).upper()
    seq_rev = str(record.seq.reverse_complement()).upper()
    seq_len = len(record.seq)
    results = []
    results += find_amplicons_on_strand(record.id, seq_fwd, f_primer, r_primer_rc, min_len, max_len, strand="+", seq_len=seq_len)
    results += find_amplicons_on_strand(record.id, seq_rev, f_primer, r_primer_rc, min_len, max_len, strand="-", seq_len=seq_len)
    return results

def main():
    parser = argparse.ArgumentParser(description="Multi-primer in-silico PCR tool")
    parser.add_argument("--fasta", required=True, help="Reference FASTA file")
    parser.add_argument("--primers", required=True, help="CSV file of primer pairs (PairID,Forward,Reverse)")
    parser.add_argument("--minlen", type=int, default=50, help="Minimum product size")
    parser.add_argument("--maxlen", type=int, default=1200, help="Maximum product size")
    parser.add_argument("--outfasta", default="amplicons.fasta", help="Output FASTA file")
    parser.add_argument("--outtable", default="amplicons.tsv", help="Output TSV table")
    args = parser.parse_args()

    with open(args.primers, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        primer_list = [(row['PairID'], row['Forward'].upper(), str(Seq(row['Reverse']).reverse_complement())) for row in reader]

    with open(args.outfasta, "w") as fasta_out, open(args.outtable, "w") as table_out:
        table_out.write("PairID\tSeqID\tStart\tEnd\tLength\tStrand\tFWD_Start\tFWD_End\tREV_Start\tREV_End\tForward_Primer\tReverse_Primer_RC\n")

        for record in SeqIO.parse(args.fasta, "fasta"):
            for pair_id, fwd, rev_rc in primer_list:
                amplicons = find_amplicons_both_strands(record, fwd, rev_rc, args.minlen, args.maxlen)
                for i, amp in enumerate(amplicons, 1):
                    fasta_out.write(f">{pair_id}_{amp['seq_id']}_amplicon_{i}_{amp['start']}_{amp['end']}_strand{amp['strand']}\n{amp['amplicon_seq']}\n")
                    table_out.write(f"{pair_id}\t{amp['seq_id']}\t{amp['start']}\t{amp['end']}\t{amp['length']}\t{amp['strand']}\t{amp['fwd_start']}\t{amp['fwd_end']}\t{amp['rev_start']}\t{amp['rev_end']}\t{amp['forward_primer']}\t{amp['reverse_primer_rc']}\n")

    print(f"\n✅ Done. Amplicons saved to {args.outfasta} and coordinates to {args.outtable}")

if __name__ == "__main__":
    main()
