from Bio import SeqIO
import numpy as np
from Sequence import *
from Generic_Methods import reverse_complement, calculate_degeneracy, generate_opportunistic_matrix
import math
import itertools
import csv

def generate_sequences(seq_path, meta_path):
    '''
    Function that reads sequences from a fasta file and saves them as Sequence objects with metadata from "metadata.tsv".

    Parameters
    ----------
    seq_path: str
        Absolute path of the sequences.
    meta_path : str
        Absolute path to the folder containing the metadata.tsv file

    Returns
    -------
    sequences : list[Sequence]
        List of sequences contained in seq_path.

    '''
    sequences_temp = {}
    to_delete = []
    
    #Read sequences from fasta file
    sequences_object = SeqIO.parse(open(seq_path), 'fasta')
    for sequence in sequences_object:
        sequences_temp[sequence.id.split('|')[0]] = str(sequence.seq.lower())
        if len(sequence.seq.replace('-','')) < 29000 or sequence.seq.lower().count('n') > 10: #require sequences to have length at least 29k
            to_delete.append(sequence.id.split('|')[0]) #store sequences that should be removed due to too many 'n's or being too short
    #Read metadata unless impossible in which case we assign every sequence its own "lineage"
    skip = -1
    try:
        sequences = []
        for meta in csv.reader(open(meta_path + '/metadata.tsv'), delimiter='\t'):
            if skip == -1: #first line is always the header line
                for cur_meta in range(len(meta)):
                    if 'lineage' in meta[cur_meta].lower():
                        skip = cur_meta
                        break
            else:
                #meta[0] contains id of the sequence
                if meta[0] not in to_delete:
                    sequences.append(Sequence(sequences_temp[meta[0].replace(' ','')], meta[0], lineage=meta[skip]))
    except: #unable to read metadata
        print('Unable to read metadata from file, making up lineages for every sequence')
        sequences = []
        i = 0
        for identifier in sequences_temp:
            sequences.append(Sequence(sequences_temp[identifier], identifier, lineage=str(i)))
            i += 1
            
    return sequences

def read_logfile(filename):
    '''
    Function that reads the logfile of an amplicon finding run and returns the found amplicons.

    Parameters
    ----------
    filename : str
        Absolute path to the log file.

    Returns
    -------
    amplicons : list[ [(int,int), dict] ]
        List of amplicons in the form of a tuple (start, end), and a dictionary containing the keys forward and reverse to store the corresponding
        forward and reverse primers.

    '''
    amplicons = []
    with open(filename, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if 'succesfully' in line:
                amplicons.append( [(int(line.split(')')[0].split('(')[-1].split(',')[0]), int(line.split(')')[0].split('(')[-1].split(',')[1])), {'forward': [], 'reverse': []}] )
    return amplicons

def read_primerfile(filename, amplicons):
    '''
    Function that reads the primerfile of an amplicon finding run and returns the amplicons and their corresponding primersets.

    Parameters
    ----------
    filename : str
        Absolute path to the primer file.
    amplicons : list[ [(int,int), dict] ]
        List of amplicons in the ofrm of a tuple (start, end), and a dictionary containing the keys forward and reverse to store the corresponding
        forward and reverse primers.

    Returns
    -------
    amplicons : list[ [(int,int), dict] ]
        List of amplicons in the form of a tuple (start, end), and a dictionary containing the keys forward and reverse which
        contain for every amplicon the corresponding forward and reverse primers.
    primerlist : dict
        Dictionary with keys forward and reverse which are simply aggregates of the forward and reverse primers of every amplicon in $amplicons.

    '''
    primerlist = {'forward': [], 'reverse': []}
    
    with open(filename, 'r') as fn:
        lines = fn.readlines()
        line_index = 0
        while line_index < len(lines):
            cur_amplicon = lines[line_index]
            if cur_amplicon.split('_')[-1][0] == 'F': #current line corresponds to forward primer
                amplicons[int(cur_amplicon.split('_')[1]) - 1][1]['forward'].append(lines[line_index+1].strip()) #assign forward primer to corresponding amplicon
                primerlist['forward'].append(lines[line_index+1].strip()) #add primer to full set of forward primers
            else: #current line corresponds to a reverse primer
                amplicons[int(cur_amplicon.split('_')[1]) - 1][1]['reverse'].append(lines[line_index+1].strip()) #assign reverse primer to corresponding amplicon
                primerlist['reverse'].append(lines[line_index+1].strip()) #add primer to full set of reverse primers
            line_index += 2 #skip over to next amplicon
        
    return amplicons, primerlist

def locate_primers(sequence, primerlist, comparison_matrix, max_degen=10):
    sequence_reverse = reverse_complement(sequence, rev=True)
    max_degen = min(max_degen, len(sequence)-sequence.count('a')-sequence.count('c')-sequence.count('g')-sequence.count('t'))
    
    def find_hits(sequence, primers, max_degen, comparison_matrix, reverse=False):
        hits = {primer: set() for primer in primers}
        hits_translated = {primer: set() for primer in primers}
        
        #Iterate over primers
        for primer in primers:
            min_stretch = math.ceil((len(primer) - max_degen) / (max_degen + 1)) #worst case scenario stretch of non-degenerate nucleotides to guide exact matching
            cur_kmers = [(primer[kmer:kmer+min_stretch], kmer) for kmer in range(len(primer) - min_stretch + 1)] #split primer into kmers based on worst case stretch
            for kmer, index in cur_kmers: #iterate over kmers and their starting index in the primer
                offset = -1
                cont = True
                while cont:
                    try:
                        offset = sequence.index(kmer, offset+1) #find next occurrence of kmer
                        hit = True
                        if offset >= index and len(sequence) >= offset + len(primer) - index: #check if kmer prefix and suffix can be appended in sequence
                            #Check if prefix in primer corresponds to prefix in sequence
                            for i in range(index):
                                if not comparison_matrix[(primer[i], sequence[offset-index+i])][0]: #not an exact match
                                    hit = False
                                    break
                            #Check if suffix in primer corresponds to suffix in sequence
                            if hit:
                                for i in range(len(primer) - min_stretch - index):
                                    if not comparison_matrix[(primer[index + min_stretch + i], sequence[offset + min_stretch + i])][0]: #not an exact match
                                        hit = False
                                        break
                            if hit: #prefix and suffix agree
                                if not reverse: #store occurrences
                                    hits[primer].add(offset - index)
                                else: #store occurrence translated to forward strand
                                    hits[primer].add(len(sequence) - (offset - index) - len(primer))
                                    
                    except: #no more occurrences of kmer in sequence
                        cont = False
        return hits
                        
    hits_fwd = find_hits(sequence, primerlist['forward'], max_degen, comparison_matrix, reverse=False)
    hits_rev = find_hits(sequence_reverse, primerlist['reverse'], max_degen, comparison_matrix, reverse=True)
    
    return hits_fwd, hits_rev
                            
def locate_amplicons(sequence, amplicons, comparison_matrix, max_degen=10, primer_length=25):
    
    binding_sites = {amplicon[0]: None for amplicon in amplicons}
    #Iterate over amplicons
    for amplicon in amplicons:
        fwd_hits, rev_hits = locate_primers(sequence, amplicon[1], comparison_matrix, max_degen=max_degen)
        amplified = (0, 10**4, False)
        
        fwd_indices = set()
        rev_indices = set()
        for fwd in fwd_hits:
            fwd_indices = fwd_indices.union(fwd_hits[fwd])
        for rev in rev_hits:
            rev_indices = rev_indices.union(rev_hits[rev])
        fwd_indices = list(fwd_indices)
        rev_indices = list(rev_indices)
        
        for fwd, rev in itertools.product(fwd_indices, rev_indices):
            if rev - fwd >= 0 and rev - fwd + primer_length < amplified[1] - amplified[0] + primer_length:
                amplified = (fwd, rev+primer_length, True)
        if amplified[2]:
            binding_sites[amplicon[0]] = amplified
            
    return binding_sites

def generate_simulationfile(sequences_path, metadata_path, logfile_path, primerfile_path, max_degen=10, primer_length=25):
    sequences = generate_sequences(sequences_path, metadata_path) #read sequences
    amplicons = read_logfile(logfile_path) #generate amplicons from logfile
    amplicons, primerlist = read_primerfile(primerfile_path, amplicons) #refine amplicons and determine corresponding primers
    M = generate_opportunistic_matrix() #generate matrix used to determine which nucleotides are identical
    
    fasta_list = []
    S = set()
    for sequence in sequences[:100]:
        print(sequence.id, sequence.lineage)
        S.add(sequence.lineage)
        realized_amplicons = locate_amplicons(sequence.sequence_raw, amplicons, M)
        for amplicon in realized_amplicons:
            if realized_amplicons[amplicon]:
                print(amplicon, realized_amplicons[amplicon][1] - realized_amplicons[amplicon][0])
            else:
                print('Amplicon not amplifiable')
    print(len(sequences))
    print(len(S))

            
    

def main():
    generate_simulationfile('/Users/jaspervanbemmelen/Downloads/Netherlands_October_2022/sequences.fasta', 
                            '/Users/jaspervanbemmelen/Downloads/Netherlands_October_2022',
                            '/Users/jaspervanbemmelen/Documents/Wastewater/Primers_Davida/coverage-1.000/misthresh20_searchwidth50_amps10_all_nseqs2800/logfile_1.txt',
                            '/Users/jaspervanbemmelen/Documents/Wastewater/Primers_Davida/coverage-1.000/misthresh20_searchwidth50_amps10_all_nseqs2800/primers_1.txt')
    sequences = generate_sequences('/Users/jaspervanbemmelen/Documents/Wastewater/source_code/amplivar/testing', '/Users/jaspervanbemmelen/Documents/Wastewater/source_code/amplivar/testing')
    M = generate_opportunistic_matrix()
    return None
    '''
    for s in sequences:
        s.align_to_trim()
    print(s.aligned_to_trim)
    '''
    #A = read_logfile('/Users/jaspervanbemmelen/Documents/Wastewater/Primers_Davida/coverage-0.999/beta-0.05/misthresh20_searchwidth50_amps10_all_nseqs2800/logfile_1.txt')
    A = read_logfile('/Users/jaspervanbemmelen/Documents/Wastewater/Primers_Davida/coverage-1.000/misthresh20_searchwidth50_amps10_all_nseqs2800/logfile_1.txt')
    print(A)
    
    #A, primerlist = read_primerfile('/Users/jaspervanbemmelen/Documents/Wastewater/Primers_Davida/coverage-0.999/beta-0.05/misthresh20_searchwidth50_amps10_all_nseqs2800/primers_1.txt', A)
    A, primerlist = read_primerfile('/Users/jaspervanbemmelen/Documents/Wastewater/Primers_Davida/coverage-1.000/misthresh20_searchwidth50_amps10_all_nseqs2800/primers_1.txt', A)
    print(A)
        
    seqs = [s.sequence_raw for s in sequences]
    lins = [s.lineage for s in sequences]
    return seqs, primerlist, M, A, lins
    '''
    i = 1
    t = time.time()
    for s in seqs:
        print(str(i) + ': ' + str(locate_primers(s, primerlist, M)))
        i += 1
    return sequences, A, primerlist
    print(time.time() - t)
    '''
    
def read_logfile(filename):
    amplicons = []
    with open(filename, 'r') as fn:
        lines = fn.readlines()
        for line in lines:
            if 'succesfully' in line:
                amplicons.append( [(int(line.split(')')[0].split('(')[-1].split(',')[0]), int(line.split(')')[0].split('(')[-1].split(',')[1])), {'forward': [], 'reverse': []}] )
                
    return amplicons

def read_primerfile(filename, amplicons):
    primerlist = {'forward': [], 'reverse': []}
    
    with open(filename, 'r') as fn:
        lines = fn.readlines()
        line_index = 0
        while line_index < len(lines):
            cur_amplicon = lines[line_index]
            if cur_amplicon.split('_')[-1][0] == 'F':
                amplicons[int(cur_amplicon.split('_')[1]) - 1][1]['forward'].append(lines[line_index+1].strip())
                primerlist['forward'].append(lines[line_index+1].strip())
            else:
                amplicons[int(cur_amplicon.split('_')[1]) - 1][1]['reverse'].append(lines[line_index+1].strip())
                primerlist['reverse'].append(lines[line_index+1].strip())
            line_index += 2
            
    return amplicons, primerlist
'''
def locate_primers(sequence, primerlist, comparison_matrix, max_degen=10):
    def all_hits(t, p):
        res = []
        offset = -1
        while True:
            try:
                offset = t.index(p, offset+1)
                res.append(offset)
            except:
                return res
            
    sequence_reverse = reverse_complement(sequence, rev=True)
    max_degen = min(max_degen, len(sequence)-sequence.count('a')-sequence.count('c')-sequence.count('g')-sequence.count('t'))
    #print(max_degen)
    hits_fwd = {primer: set() for primer in primerlist['forward']}
    hits_rev = {primer: set() for primer in primerlist['reverse']}
    hits_rev_translated = {primer: set() for primer in primerlist['reverse']}

    #Forward
    num_hits = 0
    for primer in primerlist['forward']:
        min_stretch = math.ceil((len(primer) - max_degen)/(max_degen + 1))
        cur_kmers = [(primer[kmer:kmer+min_stretch],kmer) for kmer in range(len(primer) - min_stretch + 1)]
        for kmer, index in cur_kmers:
            offset = -1
            cont = True
            while cont:
                try:
                    offset = sequence.index(kmer, offset+1)
                    hit = True
                    if offset >= index and len(sequence) >= offset + len(primer) - index:
                        for i in range(index):
                            if not comparison_matrix[(primer[i], sequence[offset-index+i])][0]:
                                hit = False
                                break
                        if hit:
                            for i in range(len(primer) - min_stretch - index):
                                if not comparison_matrix[(primer[index + min_stretch + i], sequence[offset + min_stretch + i])][0]:
                                    hit = False
                                    break
                        if hit:
                            hits_fwd[primer].add(offset - index)
                            num_hits += 1
                except:
                    cont = False
    #print('number of hits in forward primers:', num_hits)
    #print(hits_fwd)
    #Reverse
    num_hits = 0
    for primer in primerlist['reverse']:
        min_stretch = math.ceil((len(primer) - max_degen)/(max_degen + 1))
        cur_kmers = [(primer[kmer:kmer+min_stretch],kmer) for kmer in range(len(primer) - min_stretch + 1)]
        for kmer, index in cur_kmers:
            offset = -1
            cont = True
            while cont:
                try:
                    offset = sequence_reverse.index(kmer, offset+1)
                    hit = True
                    if offset >= index and len(sequence_reverse) >= offset + len(primer) - index:
                        for i in range(index):
                            if not comparison_matrix[(primer[i], sequence_reverse[offset-index+i])][0]:
                                hit = False
                                break
                        if hit:
                            for i in range(len(primer) - min_stretch - index):
                                if not comparison_matrix[(primer[index + min_stretch + i], sequence_reverse[offset + min_stretch + i])][0]:
                                    hit = False
                                    break
                        if hit:
                            hits_rev[primer].add(offset - index)
                            hits_rev_translated[primer].add(len(sequence_reverse) - (offset - index) - len(primer))
                            #To find starting index of reverse primer in original sequence: len(sequence) - start in rev comp - length : len(sequence) - start in rev comp
                            num_hits += 1
                except:
                    cont = False
    #print('number of hits in reverse primers:', num_hits)
    #print(hits_rev)
    return hits_fwd, hits_rev, hits_rev_translated
'''
        

if __name__ == '__main__':
    #sequences, A, primerlist = main()
    sequences, primerlist, M, A, lins = main()
    sequence_bindings = [{amp[0]:0 for amp in A} for s in sequences]
    
    diff = {}
    index = 0
    fasta_list = []
    for sequence in sequences[:]:
        X = locate_amplicons(sequence, [A[0]], M)
        #print('Sequence:', index+1)
        for key in X:
            break
            print(sequence[X[key][0]:X[key][1]])
        if X[key]:
            fasta_list.append('>' + lins[index])
            fasta_list.append(sequence[X[key][0]:X[key][1]])
            subseq = sequence[X[key][0]:X[key][1]]
            if not subseq in diff:
                diff[subseq] = []
            diff[subseq].append(lins[index])
        index += 1
    '''    
    diff1 = {}
    S = []
    index = 0
    with open('/Users/jaspervanbemmelen/Downloads/Netherlands_September_2022/1673523207660.sequences.fasta', 'r') as f:
        cur = ''
        for line in f:
            if '>' in line:
                if cur != '' and cur.count('N') < 5:
                    S.append(cur.lower())
                cur = ''
            else:
                cur += line.strip()
        S.append(cur.lower())
    for s in S:
        X = locate_amplicons(s, [A[3]], M)
        for key in X:
            break
        if X[key]:
            subseq = s[X[key][0]:X[key][1]]
            if not subseq in diff:
                diff1[subseq] = []
            diff1[subseq].append(index)
        index += 1
        
    
    diff2 = {}
    S = []
    index = 0
    with open('/Users/jaspervanbemmelen/Downloads/Netherlands_October_2022/1673523305918.sequences.fasta', 'r') as f:
        cur = ''
        for line in f:
            if '>' in line:
                if cur != '' and cur.count('N') < 5:
                    print(cur)
                    S.append(cur.lower())
                cur = ''
            else:
                cur += line.strip()
        if cur.count('N') < 5:
            S.append(cur.lower())
    for s in S:
        X = locate_amplicons(s, [A[3]], M)
        for key in X:
            break
        if X[key]:
            subseq = s[X[key][0]:X[key][1]]
            if not subseq in diff:
                diff2[subseq] = []
            diff2[subseq].append(index)
        index += 1
    '''
    '''
    for amplicon in A[:]:
        index = 0
        for s in sequences[:]:
            index += 1
            fwd, rev = locate_primers(s, amplicon[1], M)
            start = 0
            start_found = False
            end = len(s)
            end_found = False
            for primer in fwd:
                if len(fwd[primer]) > 0:
                    start = max(start, max(list(fwd[primer])))
                    start_found = True
            for primer in rev:
                if len(rev[primer]) > 0:
                    end = min(end, min(list(rev[primer])))
                    end_found = True
            if start_found and end_found:
                sequence_bindings[index-1][amplicon[0]] = (start, end)
                continue
                print('Amplicon: ', amplicon)
                print(f'Realized at {start}-{end} with length {end-start}')
            else:
                print(f'For sequence {index} the amplicon {amplicon[0]} cannot be amplified')
    '''
            
    
    
    
    
    
    