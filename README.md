# Adam 2.0 - Hybridization Oligo Design Utility

![Adam 2.0](./adam2.png)

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

## Installation

Adam 2.0 requires Python 3.

Installation is as easy as:

```
git clone https://github.com/FordyceLab/adam2.git
cd adam2
pip3 install .
```

This will make the `adam2` command line tool publicly available. You can check this with `which adam2`, which should not return a blank line.

## How to use it

Adam 2.0 (command line tool `adam2`) takes a requires set of IO/ parameters and has three sets of command line options used to alter its behavior:

### I/O flags
- `--input` or `-i` - specifies the input FASTA file containing sequences for oligo design **(required argument)**
- `--output` or `-o` - specifies the output file prefix **(required argument)**

### Design options
- `--size` or `-s` - desired oligo size range, must provide min and max **(required argument)**
- `--tm` or `-t` - desired oligo Tm range in degrees C, must provide min and max **(required argument)**
- `--hairpin` or `-p` - hairpin tolerance in degrees C, default = 10
- `--homodimer` or `-d` - homodimer tolerance in degrees C, default = 10
- `--number` or `-n` - number of oligos to generate per input sequence, default = 10

### Pruning options
- `--length` - length penalty for oligo scoring, default = 0.5
- `--gc_ends` - GC ends reward for oligo scoring, default = 1
- `--gc_comp` - GC composition penalty for oligo scoring, default = 2
- `--tm_mean` - Tm mean penalty for oligo scoring, default = 1
- `--hairpin_tm` - hairpin Tm suppression reward for oligo scoring, default = 0.1
- `--homodimer_tm` - homodimer Tm suppression reward for oligo scoring, default = 0.1

### Tm calculation parameters
- `--dna_conc` - DNA concentration to use for Tm calculation, default = 250
- `--mv_conc` - monovalent cation concentration to use for Tm calculation, default = 50
- `--dv_conc` - divalent cation concentration to use for Tm calculation, default = 0
- `--dntp_conc` - dNTP concentration to use for Tm calculation, default = 0

## Example

Let's say we have the following FASTA file (named `species.fasta` on disk) containing a variable region for two species below:

```
>Species_A
ACGTTGACAGGATTACACAGTAGATCCAGGATTATAGGACCAGGTAGCA

>SpeciesB
ACGTTTGACGTTGTACACAGTAGAGGCTAGGATTAGGTACCAGGTAGCA
```

We want to design a set of oligos to separate these sequences with the following parameters:

- max of 3 potential oligos designed per species
- length of 18-25 nt
- Tm of 55-60C

We can accomplish this using the following command:

```
adam2 -i species.fasta -o species -n 3 -s 18 25 -t 55 60 
```

This command will generate the following YAML file (in this case named `species.yaml`):

```
design_params:
    k: 9
    size_range: (18, 25)
    Tm_range: (55.0, 60.0)
    hairpin_tolerance: 10
    homodimer_tolerance: 10
prune_params:
    GC_ends: 1
    GC_comp: 2
    Tm_mean: 1
    hairpin_Tm: 0.1
    homodimer_Tm: 0.1
Species_A:
    ACGTTGACAGGATTACACAGTAGAT:
        length: 25
        Tm: 55.39
        start: 0
        end: 24
        hairpin_Tm: 0.0
        homodimer_Tm: -46.54
SpeciesB:
    CGTTTGACGTTGTACACAGTAGAG:
        length: 24
        Tm: 55.69
        start: 1
        end: 24
        hairpin_Tm: 32.83
        homodimer_Tm: -11.0
    CACAGTAGAGGCTAGGATTAGGTAC:
        length: 25
        Tm: 55.11
        start: 15
        end: 39
        hairpin_Tm: 0.0
        homodimer_Tm: -31.25
```

