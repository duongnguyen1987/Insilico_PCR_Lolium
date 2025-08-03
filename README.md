##############################################
#### InSilico development steps of primers
##############################################

# Tools required 

- miniconda (for instalation of relevant modules)
- ncbi-datasets-cli (to download reference genomes)
- unzip
- blast 
- Python and biopython
- mafft 
- clustalo
- pandas
- mummer

## Installation of the required modules

sudo apt install miniconda
conda install -y -c bioconda blast
conda install -c bioconda blast
conda install -c bioconda blast
sudo apt install unzip
pip install biopython
conda install -y -c bioconda blast
conda install -y -c conda-forge biopython pandas tqdm
conda install -c conda-forge ncbi-datasets-cli
conda install -c bioconda mafft
conda install -c bioconda clustalo
conda install -c bioconda mummer


## Dowload the reference genome using "datasets" function

datasets download genome accession GCF_022539505.1 --include gff3,rna,cds,protein,genome,seq-report	# rigidum

datasets download genome accession GCA_030979885.1 --include gff3,rna,cds,protein,genome,seq-report	# multiflorum

datasets download genome accession GCF_019359855.2 --include gff3,rna,cds,protein,genome,seq-report	# perenne

## Python script "inSilicoPCR_multi.py" to search amplicons in the references genomes using a list of primer pairs

python inSilicoPCR_multi.py --fasta GCF_019359855.2_Kyuss_2.0_genomic.fna --primers primer_pairs.csv --outfasta perenne_amplicons.fasta --outtable perenne_amplicons.tsv

python inSilicoPCR_multi.py --fasta GCF_022539505.1_APGP_CSIRO_Lrig_0.1_genomic.fna --primers primer_pairs.csv --outfasta rigidum_amplicons.fasta --outtable rigidum_amplicons.tsv

python inSilicoPCR_multi.py --fasta GCA_030979885.1_Rabiosa_unphased_assembly_genomic.fna --primers primer_pairs.csv --outfasta multiflorum_amplicons.fasta --outtable multiflorum_amplicons.tsv


## Python script to search for common amplicons identified from reference genomes
 
python extract_sequences_by_pairid.py --csv primer_pairs.csv --fastas multiflorum.fasta rigidum.fasta perenne.fasta 

## Python script "msa_mafft.py" to conduct alighment of the common amplicons 

python msa_mafft.py --infile LP_matK-a_sequences.fasta --outdir my_alignments --outfile LP_matK-a_sequences.al 

## Using loop script "run_all_alignments.sh" to go through the list of common amplicon and do alignment for all

bash run_all_alignments.sh 


#### NOTE THE SEQUENCCE FASTA FILE COULD BE ALIGNED AND VIEWED IN https://www.genome.jp/tools-bin/clustalw 


## Align two genome and extract SNPs

nucmer --prefix=plastid NC_009950.1.fasta NC_019651.1.fasta # Align with nucmer

delta-filter -1 plastid.delta > plastid.filtered.delta # Filter best alignments

show-snps -Clr plastid.filtered.delta > plastid.snps.txt  # Extract SNPs and Indels 

## Search location in the plastid.snps.txt and conduct alignment of region with inDel for primer design, focus on the location with high distance to the SNP/variation in either size [BUFF]   [DIST] 

python extract_indels_flanking.py \
  --ref NC_009950.1.fasta \
  --qry NC_019651.1.fasta \
  --ref_start 17992 \  # ← 1 base earlier than biological target
  --ref_end 18010 \
  --qry_center 17944 \
  --flank 150 \
  --outfile alignment_output.txt

## Etract sequence of a region in the genome
python extract_fasta_region.py \ 
--fasta NC_009950.1.fasta \ 
--seq_id NC_009950.1  \
--start 17842 \
--end 18160 \
--output perenne_17842_18160.fasta


python extract_fasta_region.py \
--fasta GCF_020047135_2021.1_ASM2004713v1_genomic.fna \
--seq_id NZ_CP077421.1 \
--start 335673 \
--end 335925 \
--output Pp-QK413-1_335673_335925.fasta
