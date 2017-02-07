# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE
for next time I should name variables better
@author: Juanita Desouza

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
        The unit tests are sufficient because they test enough cases. I added 2 to do so.
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('T')
    'A'
    >>> get_complement('G')
    'C'
    """
    if nucleotide == "A":
        return "T"
    elif nucleotide == "T":
        return "A"
    elif nucleotide == "C":
        return "G"
    elif nucleotide == "G":
        return "C"


# get_complement("A")
# get_complement("T")
# get_complement("C")
# get_complement("G")
#
#


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
        The unit tests are sufficient because they test enough cases.

    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    dna_list = list(dna)
    n = len(dna_list)
    complement_dna = [0]*n
    for x in range(0, n):
        y = dna_list[x]
        complement_dna[n-x-1] = get_complement(y)
    c = complement_dna
    b = ''.join(complement_dna)
    return (b)

#
# get_reverse_complement('ATGC')



def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.
        dna: a DNA sequence
        returns: the open reading frame represented as a string
        The unit tests are sufficient because they test enough cases. I added 2 to do so.
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    dna_list = list(dna)
    h = int(len(dna_list))
    for i in range(0, h, 3):
        try:
            codon=dna_list[i]+dna_list[i+1]+dna_list[i+2]
        except IndexError:
           break
        #codon = dna_list[i]+dna_list[i+1]+dna_list[i+2]
        if (codon == 'TAG') or (codon == 'TAA') or (codon == 'TGA'):
            b = int(i)
            dna_new = dna[0:b]
            break
        else:
            dna_n = dna_list[:]
            r = ''.join(dna_n)
            dna_new = r
        r = ''.join(dna_new)
        #print(r)
    return dna_new


def find_all_ORFs_oneframe(dna):
        """ Finds all non-nested open reading frames in the given DNA
            sequence and returns them as a list.  This function should
            only find ORFs that are in the default frame of the sequence
            (i.e. they start on indices that are multiples of 3).
            By non-nested we mean that if an ORF occurs entirely within
            another ORF, it should not be included in the returned list of ORFs.

            dna: a DNA sequence
            returns: a list of non-nested ORFs
            The unit tests are sufficient because they test a difficult cases. .
        >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
        ['ATGCATGAATGTAGA', 'ATGTGCCC']
        """
        dna_list = list(dna)
        h = int(len(dna_list))
        proteins = []
        i=0
        while i < h-2:
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
The unit tests are sufficient because they test a difficult case.
    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    allORF = []
    a = find_all_ORFs_oneframe(dna[:])
    b = find_all_ORFs_oneframe(dna[1:])
    c = find_all_ORFs_oneframe(dna[2:])
    allORF += a
    allORF += b
    allORF += c
    return allORF


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
        The unit tests are sufficient because they a difficult case.
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    allORF2_twostrand = []
    dna_rev = get_reverse_complement(dna)
    r = find_all_ORFs(dna)
    b = find_all_ORFs(dna_rev)
    allORF2_twostrand += r
    allORF2_twostrand += b
    return allORF2_twostrand
    # return r


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
        The unit tests are sufficient they test a difficult case
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    full_string = find_all_ORFs_both_strands(dna)
    x = len(full_string)
    r = 'A'
    for i in range(0, x):
        z=full_string[i]
        b=len(z)
        if len(z) > len(r):
            r = z
    return r



def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF
        tested this with a print statement"""
    dna_list = list(dna)
    z=0
    for x in range(0,num_trials):
        random.shuffle(dna_list)
        red= ''.join(dna_list)
        r=len(longest_ORF(red))
        if r > z:
            z=r
    return z

#print(longest_ORF_noncoding("AAATATATGCGAATGTAGCATCAAA", 3))

def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment
        The unit tests are sufficient because they test enough cases.

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    h= len(dna)
    dna_list= list(dna)
    aminoacid=[]
    aminoacid_s =""
    for i in range(0, h, 3):
        try:
            codon=dna_list[i]+dna_list[i+1]+dna_list[i+2]
        except IndexError:
            break
        amino_acid=aa_table[codon]
        aminoacid+=amino_acid
        aminoacid_s=''.join(aminoacid)
    return aminoacid_s
    # i=0
    # while i < h-2:
    # codon=dna_list[i]+dna_list[i+1]+dna_list[i+2]


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely
coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences
coded by the sequence dna.
tested this on 200
    """
    threshold = longest_ORF_noncoding(dna, 1500)
    allORF1= find_all_ORFs_both_strands(dna)
    print(allORF1)
    l= len(allORF1)
    ORfaftert=[]
    for x in range(0,l):
        if len(allORF1[x])>threshold:
            n=allORF1[x]
            ORfaftert.append(n)
    am_ret=[]
    o= len(ORfaftert)
    print(ORfaftert)
    for i in range(0, o):
        z=coding_strand_to_AA(ORfaftert[i])
        print(z)
        am_ret.append(z)
    return am_ret


if __name__ == "__main__":
    import doctest
    doctest.testmod()
    #doctest.run_docstring_examples(coding_strand_to_AA, globals(), verbose=True)
    from load import load_seq
    dna = load_seq("./data/X73525.fa")
    print(gene_finder(dna))
