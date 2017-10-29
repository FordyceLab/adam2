# from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from itertools import chain
from collections import Counter, OrderedDict
import numpy as np
from distance import hamming
import primer3
from tqdm import *
import math
import datetime
import yaml

def setup_yaml():
    # https://stackoverflow.com/a/8661021
    represent_dict_order = lambda self, data:  self.represent_mapping('tag:yaml.org,2002:map', data.items())
    yaml.add_representer(OrderedDict, represent_dict_order)

setup_yaml()

class TargetSeq():
    """
    Target sequence against which to make specific hybridization oligo

    Args:
    - seq (Seq) - Biopython Seq object corresponding to (DNA or RNA) sequence of interest
    - size_range (tuple of ints) - minimum and maximum oligo sizes to target
    - tm_range (tuple of floats) - minimum and maximum melting temperatures to target
    - hairpin_tolerance (int) - minimum Tm difference between hairpin and oligo Tm
    - homodimer_tolerance (int) - minimum Tm difference between homodimer and oligo Tm
    """
    def __init__(self, seq, size_range, tm_range, hairpin_tolerance, homodimer_tolerance):
        # Add the sequence record
        self.record = seq
        
        # Get both the DNA and RNA sequences
        self.dna_seq = str(seq.seq.back_transcribe())
        self.rna_seq = str(seq.seq.transcribe())
        
        # Get the length of the total sequence
        self.length = len(self.dna_seq)
        
        # Get the target size range, tm_range, and hairpin tolerance
        self.size_range = size_range
        self.tm_range = tm_range
        self.hairpin_tolerance = hairpin_tolerance
        self.homodimer_tolerance = homodimer_tolerance
    
    def count_kmers(self, k):
        """
        Count all k-mers in the target sequence
        
        Args:
        - k (int) - size of k-mer to use
        
        Return:
        - Counter object corresponding to the counts of each k-mer in the sequence
        """
        counts = Counter([str(self.dna_seq[i:i+k]) for i in range(len(self.dna_seq)-(k-1))])
        return counts
        
    def occ_dist(self, counts):
        """
        Get the distribution of k-mer occurances in the sequence
        
        Args:
        - counts (Counter) - Counter object corresponding to the counts of each k-mer in the sequence
        
        Return:
        - Counter object corresponding to the counts of each k-mer in the sequence
        """
        occ_dist = Counter(counts.values())
        return occ_dist
    
    def choose_k(self):
        """
        Choose the value of k that splits the sequence into all unique k-mers
        
        Return:
        - Value of k for which the sequence is split into only unique k-mers
        """
        # Start with k = 2
        k = 2
        proceed = True
        
        # Keep increasing k until all k-mers are unique
        while proceed:
            kmer_counts = self.count_kmers(k)
            occ_dist = self.occ_dist(kmer_counts)
            occ_dist_keys = [int(key) for key in occ_dist.keys()]
            if any(count > 1 for count in occ_dist_keys):
                k += 1
            else:
                proceed = False
        
        return k
    
    def kmerize(self, k):
        """
        Break sequence into k-mers
        
        Args:
        - k (int) - size of k-mer
        """
        self.counts = self.count_kmers(k)
        self.kmers = set(self.counts.keys())
        
    
    def generate_oligos(self, n, dna_conc=250, mv_conc=50, dv_conc=0, dntp_conc=0):
        """
        Generate oligo set for target sequence
        
        Args:
        - n (int) - max number of seeds to generate oligos for,
        starting from most unique seed and working to least
        """
        self.oligo_sets = []
    
        index = 0
            
        while len(self.oligo_sets) < n:
            try:
                seed = self.unique_kmer_list[index]
                oligo_set = OligoSet(seed, self, dna_conc, mv_conc, dv_conc, dntp_conc)
                if oligo_set.valid_oligos:
                    self.oligo_sets.append(oligo_set)
            except IndexError:
                break
                
            index += 1
            
    def prune_oligos(self, length=0.5, GC_ends=1, GC_comp=2, Tm_mean=1, hairpin_Tm=0.1, homodimer_Tm=0.1):
        """
        Prune an oligo set to only the "best" oligo
        
        Args:
        - length (float) - penalty for additional nucleotides beyond minimum specified length, default = 0.5
        - GC_ends (float) - reward for oligo starting/ending with G/C, default = 1
        - GC_comp (float) - penalty for GC composition deviating from 0.5, default = 2
        - Tm_mean (float) - penalty for Tm per degree of deviation from mean of desired range, default = 1
        - hairpin_Tm (float) - reward for hairpin Tm per degree of deviation minimum threshold, default = 1
        - homodimer_Tm (float) - reward for homodimer Tm per degree of deviation minimum threshold, default = 0.1
        """
        # Set an instance variable for the 
        self.pruned_oligos = []
                
        for oligo_set in self.oligo_sets:
            max_oligo_score = -math.inf
            best_oligo = None
            
            for oligo in oligo_set.valid_oligos:
                length_score = -length * (len(oligo.seq) - self.size_range[0])
                gc_end_score = GC_ends * (float(oligo.seq.startswith(('C', 'G'))) +
                                          float(oligo.seq.endswith(('C', 'G'))))
                gc_comp_score = -GC_comp * abs((oligo.seq.count('G') + oligo.seq.count('G') / len(oligo.seq))- 0.5)
                tm_mean_score = -Tm_mean * abs(oligo.Tm - (self.tm_range[0] + self.tm_range[1] / 2))
                hairpin_tm_score = hairpin_Tm * abs(self.hairpin_tolerance - oligo.hairpin_Tm)
                homodimer_tm_score = homodimer_Tm * abs(self.homodimer_tolerance - oligo.homodimer_Tm)
                
                oligo_score = gc_end_score + gc_comp_score + tm_mean_score + hairpin_tm_score + homodimer_tm_score
                                
                if oligo_score > max_oligo_score:
                    best_oligo = oligo
                    max_oligo_score = oligo_score
                
            self.pruned_oligos.append(best_oligo)
            
    def report(self):
        target_dict = OrderedDict()
        
        for oligo in self.pruned_oligos:
            target_dict[oligo.seq] = oligo.report()
            
        return target_dict
        
        
class TargetSet():
    """
    Set of target sequences to be separated by hybridization oligos

    Args:
    - target_list (list of TargetSeqs) - list of target sequences to be separated
    """
    
    def __init__(self, seq_list, size_range, tm_range, hairpin_tolerance, homodimer_tolerance):
        self.target_list = [TargetSeq(seq, size_range, tm_range, hairpin_tolerance, homodimer_tolerance)
                            for seq in seq_list]
        self.size_range = size_range
        self.tm_range = tm_range
        self.hairpin_tolerance = hairpin_tolerance
        self.homodimer_tolerance = homodimer_tolerance
        self.choose_k()
        self.get_unique_kmers()
        self.assign_kmers()
        self.gen_hamming_matrices()
        self.get_most_unique_kmers()
        
    def choose_k(self):
        """
        kmerize all target sequences using the value of k that breaks all sequences in the set into
        completely unique k-mers
        """
        # Get all values of k and take the max
        ks = [kmer.choose_k() for kmer in self.target_list]
        self.k = max(ks)
        
        # Break all sequences into constituitive k-mer parts
        for target in self.target_list:
            target.kmerize(self.k)
    
    def get_unique_kmers(self):
        """
        Generate sets of all k-mers found and all unique k-mers (k-mers found in only one target species)
        """
        # Figure out how many times each k-mer appears in the master set of all kmers
        super_list = list(chain(*[target.kmers for target in self.target_list]))
        super_counter = Counter(super_list)
        
        # Get the k-mers that only appear once as well as all observed kmers
        self.unique_kmers = {kmer for kmer in super_counter.keys() if super_counter[kmer] == 1}
        self.all_kmers = {kmer for kmer in super_counter.keys()}
        
    def assign_kmers(self):
        """
        Assign the k-mers from the unique pool to the corresponding parental sequence
        """
        # Assign the unique kmers back to the parental sequence
        for target in self.target_list:
            target.unique_kmers = target.kmers.intersection(self.unique_kmers)
            
    def gen_hamming_matrices(self):
        """
        Generate Hamming distance and similarity matrices for all sequences
        """
        # For all target sequences
        for target in tqdm(self.target_list, desc="Generating Hamming distance matrices"):
            
            # Sort the list of k-mers unique to the target and all observed kmers in the set alphabetically
            unique_kmer_list = sorted(list(target.unique_kmers))
            other_kmer_list = sorted(list(self.all_kmers.difference(target.unique_kmers)))
            
            # Get the number of unique k-mers in current target sequence and all observed kmers
            target_unique = len(unique_kmer_list)
            other_kmers = len(other_kmer_list)
            
            # Generate a matrix of appropriate size
            ham_array = np.zeros((target_unique, other_kmers))
            
            # Loop over both lists and fill in the Hamming distance matrix
            for i in range(target_unique):
                for j in range(other_kmers):
                    ham_array[i, j] = hamming(unique_kmer_list[i], other_kmer_list[j])
                 
            # Assign the lists as instance variables to the target sequences in the list
            target.unique_kmer_list = unique_kmer_list
            target.other_kmer_list = other_kmer_list
            
            # Assign the Hamming distance and similarity score matrices to the target sequences in the list 
            target.hamming_dist_matrix = ham_array
            target.similarity_matrix = np.full_like(ham_array, self.k) - ham_array
            
    def get_most_unique_kmers(self):
        """
        Sort the unique k-mers, penalized by squared similarity score
        """
        for target in self.target_list:
            target.dist_sort = np.argsort(np.sum(np.power(target.similarity_matrix, 2), axis = 1))
            
    def generate_oligos(self, n, dna_conc=250, mv_conc=50, dv_conc=0, dntp_conc=0):
        """
        Generate oligo sets for all targets in the set
        
        Args:
        - n (int) - max number of seeds to generate oligos for,
        starting from most unique seed and working to least
        """
        self.n_oligos = n
        
        for target in tqdm(self.target_list, desc="Generating oligos"):
            target.generate_oligos(n, dna_conc=250, mv_conc=50, dv_conc=0, dntp_conc=0)
            
    def prune_oligos(self, length=0.5, GC_ends=1, GC_comp=2, Tm_mean=1, hairpin_Tm=0.1, homodimer_Tm=0.1):
        """
        Prune an oligo set to only the "best" oligo
        
        Args:
        - length (float) - penalty for additional nucleotides beyond minimum specified length, default = 0.5
        - GC_ends (float) - reward for oligo starting/ending with G/C, default = 1
        - GC_comp (float) - penalty for GC composition deviating from 0.5, default = 2
        - Tm_mean (float) - penalty for Tm per degree of deviation from mean of desired range, default = 1
        - hairpin_Tm (float) - reward for hairpin Tm per degree of deviation minimum threshold, default = 1
        - homodimer_Tm (float) - reward for homodimer Tm per degree of deviation minimum threshold, default = 0.1
        """
        
        self.length = length
        self.GC_ends = GC_ends
        self.GC_comp = GC_comp
        self.Tm_mean = Tm_mean
        self.hairpin_Tm = hairpin_Tm
        self.homodimer_Tm = homodimer_Tm
        
        for target in self.target_list:
            target.prune_oligos(length, GC_ends, GC_comp, Tm_mean, hairpin_Tm, homodimer_Tm)
            
    def report(self):
        target_set_dict = OrderedDict()
        
        target_set_dict['design_params'] = OrderedDict({
            'k': self.k,
            'size_range': str(self.size_range),
            'Tm_range': str(self.tm_range),
            'hairpin_tolerance': self.hairpin_tolerance,
            'homodimer_tolerance': self.homodimer_tolerance
        })
        
        target_set_dict['prune_params'] = OrderedDict({
            'GC_ends': self.GC_ends,
            'GC_comp': self.GC_comp,
            'Tm_mean': self.Tm_mean,
            'hairpin_Tm': self.hairpin_Tm,
            'homodimer_Tm': self.homodimer_Tm
        })
        
        for target in self.target_list:
            target_set_dict[target.record.name] = target.report()
        
        return yaml.dump(target_set_dict, default_flow_style=False, indent = 4)
            
                
class Oligo():
    """
    Class to hold all information regarding a specific primer
    
    Args:
    - seq (str) - string containing the primer sequence to be used
    - dna_conc (float) - concentration (nM) of oligo to use for Tm analysis, default is 250
    - mv_conc (float) - concentration (mM) of monovalent cations to use for Tm analysis, default is 50
    - dv_conc (float) - concentration (mM) of divalent cations to use for Tm analysis, default is 0
    - dntp_conc (float) - concentration (mM) of dNTPs to use for Tm analysis, default is 0
    """
    
    def __init__(self, seq, dna_conc=250, mv_conc=50, dv_conc=0, dntp_conc=0):
        # Convert the sequence to uppercase
        seq = seq.upper()
        
        # Set the sequence to an instance variable
        self.seq = seq
        
        # Set the standard IDT OligoAnalyzer settings for analysis by primer3
        self.thermo_settings = primer3.thermoanalysis.ThermoAnalysis(dna_conc, mv_conc, dv_conc, dntp_conc)
        
        # Get the Tm of the oligo
        self.Tm = self.thermo_settings.calcTm(self.seq)
        
        # Get hairpin and homodimer properties
        self.hairpin()
        self.homodimer()
        
    def hairpin(self):
        """
        Calculate hairpin properties of the oligo sequence
        """
        hairpin_result = self.thermo_settings.calcHairpin(self.seq)
        self.hairpin = hairpin_result.structure_found
        self.hairpin_Tm = hairpin_result.tm
        
    def homodimer(self):
        """
        Calculate homodimer properties of the oligo sequence
        """
        homodimer_result = self.thermo_settings.calcHomodimer(self.seq)
        self.homodimer = homodimer_result.structure_found
        self.homodimer_Tm = homodimer_result.tm
        
    def score(self):
        pass
        
    def report(self):
        details = OrderedDict({
            'length': len(self.seq),
            'Tm': round(self.Tm, 2),
            'start': self.start,
            'end': self.end,
            'hairpin_Tm': round(self.hairpin_Tm, 2),
            'homodimer_Tm': round(self.homodimer_Tm, 2)
        })
        
        return details
    
    def __repr__(self):
        return 'Seq: {}; Tm: {:.2f}'.format(self.seq, self.Tm)
    
    def __str__(self):
        return self.__repr__()
            

class OligoSet():
    """
    Generate oligo set for seed and sequence

    Args:
    - seed (str) - seed string to incorporate
    - target (TargetSeq) - target
    """
    
    def __init__(self, seed, target, dna_conc=250, mv_conc=50, dv_conc=0, dntp_conc=0):
        self.seed = seed
        self.target = target
        self.length = len(seed)
        
        self.start = target.dna_seq.find(seed)
        self.end = self.start + self.length
        
        self.valid_oligos = []
        
        self.dna_conc = dna_conc
        self.mv_conc = mv_conc
        self.dv_conc = dv_conc
        self.dntp_conc = dntp_conc
        
        for size in range(self.target.size_range[0], self.target.size_range[1] + 1):
            self.sliding_window(size)
        
    def sliding_window(self, size):
        """
        Perform sliding window oligo optimization
        
        Args:
        - size (int) - size of the sliding window to use
        """
        # Get the maximum overhang length (difference between the given size and the seed)
        max_overhang = size - self.length
        
        # Make sure desired size is larger than the size of the seed sequence
        if max_overhang < 0:
            raise ValueError("Desired oligo size is less than the length of the oligo seed")
        
        # Get the initial starting point of the sliding window
        if self.start - max_overhang >= 0:
            init_start = self.start - max_overhang
        else:
            init_start = 0
            init_start = 0
            
        # Get the final endpoint of the sliding window
        if self.end + max_overhang < self.target.length:
            final_end = self.end + max_overhang
        else:
            final_end = self.target.length
            
        # Slide the window across the sequence
        for i in range(init_start, final_end - size + 1):
            oligo_seq = self.target.dna_seq[i:i + size]
            oligo = Oligo(oligo_seq, self.dna_conc, self.mv_conc, self.dv_conc, self.dntp_conc)
            
            oligo.start = i
            oligo.end = i + (size - 1)
            
            # Save the oligo sequence to the list of valid oligos if it meets selection criteria
            if self.target.tm_range[0] <= oligo.Tm <= self.target.tm_range[1] and \
                oligo.Tm - oligo.hairpin_Tm > self.target.hairpin_tolerance and \
                oligo.Tm - oligo.homodimer_Tm > self.target.homodimer_tolerance:
                
                self.valid_oligos.append(oligo)