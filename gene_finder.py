# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: YOUR NAME HERE

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
    if nucleotide == "A":
        #print("T")
        return "T"
    elif nucleotide == "T":
            #print("A")
        return "A"
    elif nucleotide == "C":
            #print("G")
            return "G"
    elif nucleotide == "G":
            #print("C")
            return "C"


# get_complement("A")
# get_complement("T")
# get_complement("C")
# get_complement("G")
#
#
    # """ Returns the complementary nucleotide
    #
    #     nucleotide: a nucleotide (A, C, G, or T) represented as a string
    #     returns: the complementary nucleotide
    # >>> get_complement('A')
    # 'T'
    # >>> get_complement('C')
    # 'G'
    # """
    # # TODO: implement this
    # pass


def get_reverse_complement(dna):
    dna_list = list(dna)
    n = len(dna_list)
    #print(n)
    complement_dna = [0]*n
    for x in range(0, n):
        y = dna_list[x]
        #print(y)(
        complement_dna[x] = get_complement(y)
    c=complement_dna
    #print(c)
    b = ''.join(complement_dna)
    #print b
    return (b)

#
# get_reverse_complement('ATGC')

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
    pass


def rest_of_ORF(dna):
    # """ Takes a DNA sequence that is assumed to begin with a start
    #         codon and returns the sequence up to but not including the
    #         first in frame stop codon(TAG, TGA, TAA).  If there is no in frame stop codon,
    #         returns the whole string.
    #
    #         dna: a DNA sequence
    #         returns: the open reading frame represented as a string
    # >>> rest_of_ORF("ATGTGAA")"""
    dna_list = list(dna)
    h = int(len(dna_list))
    for i in range(0, h, 3):
        #print(i)
        try:
            codon=dna_list[i]+dna_list[i+1]+dna_list[i+2]
        except IndexError:
           break
        #codon = dna_list[i]+dna_list[i+1]+dna_list[i+2]
        #print(codon)
        if (codon == 'TAG') or (codon == 'TAA') or (codon == 'TGA'):
            b = int(i)
            dna_new=dna[0:b]
            #print("found stop codon")
            #return dna_new
            break
        else:
            dna_n=dna_list[:]
            r= ''.join(dna_n)
            dna_new=r
            #return dna_new
            #print('a')
    #print(dna_new
        r = ''.join(dna_new)
        #print(r)
    print(dna_new)
    return dna_new


rest_of_ORF('ATCTATTAT')
rest_of_ORF('ATCTAAAAA')
    # """ Takes a DNA sequence that is assumed to begin with a start
    #     codon and returns the sequence up to but not including the
    #     first in frame stop codon(TAG, TGA, TAA).  If there is no in frame stop codon,
    #     returns the whole string.
    #
    #     dna: a DNA sequence
    #     returns: the open reading frame represented as a string
    # >>> rest_of_ORF("ATGTGAA")
    # 'ATG'
    # >>> rest_of_ORF("ATGAGATAGG")
    # 'ATGAGA'
    # """
    # # TODO: implement this
    # pass


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
        dna_list = list(dna)
        h = int(len(dna_list))
        proteins = []
        i=0
        while i < h:
            codon=dna_list[i]+dna_list[i+1]+dna_list[i+2]
            if (codon == 'ATG'):
                protein_string = rest_of_ORF(dna[i:])
                proteins.append(protein_string)
                i=i+len(protein_string)
            else:
                i=i+3

        return proteins



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


# def gene_finder(dna):
#     """ Returns the amino acid sequences that are likely
# coded by the specified dna
#
#         dna: a DNA sequence
#         returns: a list of all amino acid sequences
# coded by the sequence dna.
#     """
#     # TODO: implement this
#     pass


if __name__ == "__main__":
    import doctest
    #doctest.testmod()
    doctest.run_docstring_examples(find_all_ORFs_oneframe, globals(), verbose=True)
