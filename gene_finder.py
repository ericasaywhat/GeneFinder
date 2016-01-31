# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Erica Lee

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    # TODO: implement this

    if nucleotide == "A":
        return "T"
    elif nucleotide == "C":
        return "G"
    elif nucleotide == "T":
        return "A"
    elif nucleotide == "G":
        return "C"

    else:
        return "x"

    pass


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
    # TODO: implement this
    # string = ""
    # for i in range(0, len(dna)-1):
    #     string += get_complement(dna[i])
    # print "'" + string[::-1] + "'"
    string = [] 

    for nucleotide in dna:
            string.append(get_complement(nucleotide))

    return''.join(list(reversed(string))) 


    pass


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
    # TODO: implement this
    i = dna.find("ATG")

    index = ["TAG", "TAA", "TGA"]
    for stop in index:
        for x in range(i, len(dna), 3):
            if dna[x:x+3] == stop:
                # print(x)
                return dna[i:x]
    else:
        return dna[i:]


    pass


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
    # TODO: implement this
    i = 0
    ORFs = []
    x = len(rest_of_ORF(dna)) 
    ORFs.append(rest_of_ORF(dna))
    for y in range(x, len(dna)-3,3):
        ORFs.append(rest_of_ORF(dna[x:]))
        y = len(rest_of_ORF(dna[y:])) 
        print x

    # x = 0
    # while x < len(dna) - 2:
    #     ORFs.append(rest_of_ORF(dna[i:]))
    #     x = len(rest_of_ORF(dna[i:])) 
    #     print x
    #     print dna[x:]
    #     if "ATG" in dna[x:x+2]:
    #         i = dna[x:].find("ATG")
    #     print i
    #     if i == -1:
    #         return ORFs





    # while len(dna) > 3:
    #     ORFs.append(rest_of_ORF(dna))
    #     end_index = len(rest_of_ORF(dna))
    #     ORFs.append
    #     # dna = dna.replace(rest_of_ORF(dna),"")
    #     print dna
    return ORFs
    pass
print find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")

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
    print dna
    ORFs = []
    index = 0
    while len(dna) > 3:
        dna = dna[index:]
        ORFs.extend((find_all_ORFs_oneframe(dna)))
        dna = dna[1:]
        index = dna.find("ATG")
        if index == -1:
            break

    

    return ORFs
    pass


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    # TODO: implement this
    ORFs = []

    ORFs.extend(find_all_ORFs(dna))
    print ORFs

    ORFs.extend(find_all_ORFs(get_reverse_complement(dna)))

    return ORFs
    pass

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    pass


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    pass


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
    # TODO: implement this
    pass


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    pass

# if __name__ == "__main__":
#     import doctest
#     doctest.run_docstring_examples(find_all_ORFs_both_strands, globals())


# print find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
# print find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")