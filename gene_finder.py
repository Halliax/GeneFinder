# -*- coding: utf-8 -*-
"""
Given a DNA sequence, prints a list of the protein sequences encoded by all
likely coding regions.

@author: Taylor Sheneman

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

    tests: test each use case. All other cases throw an exception.

    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('T')
    'A'
    >>> get_complement('G')
    'C'
    """
    nucleotides = ['A','C','T','G']
    complements = ['T','G','A','C']
    if nucleotide in nucleotides:
        return complements[nucleotides.index(nucleotide)]
    else:
        raise KeyError("There is an invalid nucleotide in that string!")



def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string

    tests: These seem fine. All codons are used, lengths vary.
    Still, added an empty string test.

    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    >>> get_reverse_complement("")
    ''
    """
    # reverse_complements = ''
    # for nucleotide in dna:
    #     current = get_complement(nucleotide)
    #     reverse_complements = current + reverse_complements
    reverse_complements = [get_complement(dna[i]) for i in range(len(dna)-1,-1,-1)]
    return ''.join(reverse_complements)



def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string

    tests: Correct use cases with and without stop codons.

    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATGAGATCAGAC")
    'ATGAGATCAGAC'
    """
    stop_codons = ['TAG','TAA','TGA']
    i = 0
    while i < len(dna):
        if dna[i:i+3] in stop_codons:
            return dna[:i]
        i += 3
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

    tests: Added a nested case.

    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("ATGCATGAAATGAGATAGATGTGCCC")
    ['ATGCATGAAATGAGA', 'ATGTGCCC']
    """
    i = 0
    one_ORF_list = []
    while i < len(dna):
        if dna[i:i+3] == 'ATG':
            ORF = rest_of_ORF(dna[i:])
            one_ORF_list.append(ORF)
            i += len(ORF)
        else:
            i += 3
    return one_ORF_list



def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    tests: This comprehensively tests the reading frame shift. I don't have
    additional tests to run.

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    # ORF_list = []
    # for i in range(3):
    #     ORF_list = ORF_list + find_all_ORFs_oneframe(dna[i:])
    ORF_list = sum([find_all_ORFs_oneframe(dna[i:]) for i in range(3)],[])
    return ORF_list



def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    tests: This tests the use of reverse complement. All other cases are covered
    by previous functions.

    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    ORF_master_list = find_all_ORFs(dna) + find_all_ORFs(get_reverse_complement(dna))
    return ORF_master_list



def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string

    tests: Two sequences with varying numbers of ORFs, and one with none.

    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'

    >>> longest_ORF("ATGCATGAATGTAGATAGATGTGCCC")
    'ATGAATGTAGATAGATGTGCCC'

    >>> longest_ORF("CAGGAACGTAG")
    ''
    """
    ORFs = find_all_ORFs_both_strands(dna)
    if len(ORFs) == 0:
        return ''
    lengths = [len(ORF) for ORF in ORFs]
    return ORFs[lengths.index(max(lengths))]



def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF

    tests: cannot do. The output is randomized. Instead, I printed out a bunch of
    results to see if they seemed reasonable.

    """
    lengths = [len(longest_ORF(shuffle_string(dna))) for i in range(num_trials)]
    # print longest_ORF(shuffle_string(dna))
    # print lengths
    return max(lengths)



def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

    tests: Added a longer case by manually looking up the expected output.
    Included stop codon, for kicks.

    >>> coding_strand_to_AA("ATGCGA")
    'MR'

    >>> coding_strand_to_AA("ATGCCCGCTTT")
    'MPA'

    >>> coding_strand_to_AA("ATGCATGAATGTAGATAGATGTGCCC")
    'MHECR|MC'
    """
    aminos = [aa_table[(dna[i:i+3])] for i in range(0,len(dna),3)
        if len(dna[i:i+3])==3]
    return ''.join(aminos)



def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.

    tests: tricky to run tests for this one too, because the threshold case is
    random. Printed the threshold value to elucidate run-to-run return differences.
    """

    threshold = longest_ORF_noncoding(dna,1500)
    print "Randomly calculated minimum length threshold:",threshold,"base pairs"
    sequences = [coding_strand_to_AA(x)
        for x in find_all_ORFs_both_strands(dna)
        if len(x) > threshold]
    return sequences



if __name__ == "__main__":
    import doctest
    doctest.testmod()
    # doctest.run_docstring_examples(coding_strand_to_AA, globals())
    dna = load_seq("./data/X73525.fa")
    print gene_finder(dna)
