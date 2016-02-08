# -*- coding: utf-8 -*-
"""
Last Updated: January 31,2016

This is a gene finding Python program that can accurately determine regions 
of the Salmonella bacterium's DNA that code for proteins.

@author: Erica Lee

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq

from load import load_seq, load_contigs
dna = load_seq("./data/X73525.fa")
contigs = load_contigs()


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))



def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """

    if nucleotide == "A":
        return "T"
    elif nucleotide == "C":
        return "G"
    elif nucleotide == "T":
        return "A"
    elif nucleotide == "G":
        return "C"

    else:
        return "A"



def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """

    string = [] 

    for nucleotide in dna:
            string.append(get_complement(nucleotide))

    return''.join(list(reversed(string))) 




def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    stop = ["TAG", "TAA", "TGA"]

    x = 0

    while x < len(dna):
        if dna[x:x+3] in stop:
            return dna[0:x]
        x+=3
    return dna

def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """

    ORFs = []
    i = 0

    while i < len(dna):
        if dna[i:i+3] == "ATG":
            ORFs.append(rest_of_ORF(dna[i:]))
            i += len(rest_of_ORF(dna[i:]))
        i+=3


    return ORFs


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    # TODO: implement this

    ORFs = []
    ORFs.extend((find_all_ORFs_oneframe(dna[0:])))
    ORFs.extend((find_all_ORFs_oneframe(dna[1:])))
    ORFs.extend((find_all_ORFs_oneframe(dna[2:])))



    return ORFs


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    ORFs = []

    ORFs.extend(find_all_ORFs(dna))

    ORFs.extend(find_all_ORFs(get_reverse_complement(dna)))

    return ORFs



def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    i = 0
    ORFs_both = find_all_ORFs_both_strands(dna)
    longest = ORFs_both[0]
    while i +1 < len(ORFs_both):
        if len(ORFs_both[i+1]) > len(longest):
            longest = ORFs_both[i+1]
        i+=1
    return longest


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    longest = 0
    number = 0
    while number < num_trials:
        print number
        string = shuffle_string(dna)
        if len(str(longest_ORF(string))) > longest:
            longest = len(str(longest_ORF(string)))

        number += 1

    return longest





def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    amino_acid = ''
    i = 0
    while i+3 < len(dna)+1:
        amino_acid += aa_table[dna[i:i+3]]
        i+=3
    return amino_acid


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """

    threshold = longest_ORF_noncoding(dna, 50)
    
    all_both = find_all_ORFs_both_strands(dna)
    amino = []
    for ORFs in all_both:
        if len(ORFs) > threshold:
            amino.append(coding_strand_to_AA(ORFs))
    
    return amino





if __name__ == "__main__":
    import doctest
    doctest.run_docstring_examples(gene_finder, globals())

    print gene_finder(contigs[1][1])



# print find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
# print find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")