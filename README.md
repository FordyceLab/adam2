# Adam 2.0 - Hybridization Oligo Design Utility

## How it works

Adam 2.0 is a design utility for hybridization oligos. Given a pool of sequences, provided as a file in FASTA format, the utility outputs putative identification primers that meet the user's design specifications. The program uses k-mer based approach that works in the following way:

1. Given a set of sequences, find the minimum k-mer length that splits all sequences up into unique, unambiguous k-mers
2. Break all sequences into unique k-mers using the identified value of *k*
3. Remove from consideration all k-mers that occur multiple times within the collection of sequences of interest (only consider unique kmers)
4. Generate a matrix of Hamming distances between all unique k-mers within a single sequence and all other k-mers from all other sequences
5. Remove all rows that contain a Hamming distance of 1 for any k-mer pair (unique k-mers must be at least 2 mutations from any other k-mer in the pool)
6. Map Hamming distances to similarity scores S by calculating *S = k−Hamming distance* and minimize the sum of the squared similarity scores
7. Find k-mer within original sequence
8. Perform sliding window primer design of varying size until desired Tm and other properties are achieved (check with primer3) to generate a pool of potential hybridization oligos
9. Prune the oligo pool to a single representative oligo based on pruning properties
10. Report the hybridization oligos in a machine/human readable format
