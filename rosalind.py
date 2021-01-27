# Kareem Belle,  2017 

def nuc_count(sequence):
    """
    returns the nucleotide count for a DNA sequence
    """
    a_count = 0
    g_count = 0
    c_count = 0
    t_count = 0

    for nucleotide in sequence:
        if nucleotide == "A":
            a_count += 1
        elif nucleotide == "G":
            g_count += 1
        elif nucleotide == "C":
            c_count += 1
        elif nucleotide == "T":
            t_count += 1
    return a_count, c_count, g_count, t_count


def transcribe(sequence):
    """
    transcribes DNA sequence to
    """
    rna = ""
    for nucleotide in sequence:
        if nucleotide == "T":
            rna += "U"
        else:
            rna += nucleotide
    return rna


def compliment_seq(sequence):
    """
    returns the reverse compliment DNA sequence
    """

    compliment = ""
    for nucleotide in sequence:
        if nucleotide == "A":
            compliment += "T"
        elif nucleotide == "G":
            compliment += "C"
        elif nucleotide == "C":
            compliment += "G"
        elif nucleotide == "T":
            compliment += "A"
    return compliment[::-1]


def fasta_to_dict(path):
    """
    simple fasta converter for rosalind problems
    given a path to a fasta file returns a dictionary where the key
    is the identifier and the sequence is the value
    """
    fast_dict = {}
    fasta_raw = path.readlines()
    for val in fasta_raw:
        ls = ''
        if val[0] == '>':
            fast_dict[val[1:]] = ''
        else:
            for key in fast_dict:
                ls += str(val.replace('\n', ''))
                fast_dict[key] += ls
    return fast_dict


def gc_content(path):
    """given a fasta file retuns a dictionary where each identifier
    is the key and the value is the GC content percentage"""
    fast_dict = fasta_to_dict(path)
    count_dict = {}
    # split fasta sequences into dictionary
    for key in fast_dict:
        count_dict[key] = ((nuc_count(fast_dict[key])[1] +
                            nuc_count(fast_dict[key])[2]) /
                           sum(nuc_count(fast_dict[key])) *
                           100)

    return count_dict


def translate(sequence):
    """
    returns a protein sequence given an rna sequence

    """
    protein = ''
    codon = ''

    table = {"UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
             "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
             "UAU": "Y", "UAC": "Y", "UAA": "STOP", "UAG": "STOP",
             "UGU": "C", "UGC": "C", "UGA": "STOP", "UGG": "W",
             "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
             "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
             "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
             "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
             "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
             "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
             "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
             "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
             "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
             "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
             "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
             "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G", }

    for i in sequence:
        codon += i
        if len(codon) == 3 and table[codon] != "STOP":

            protein += table[codon]
            codon = ''
        else:
            continue
    return protein


def count_point_mutations(str1, str2):
    """
    counts the nucleotide differences given 2 sequences

    """
    mistakes = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            mistakes += 1
    return mistakes


def protein_mass_calc(amino_seq):
    """
    returns the mass of  protein in Daltons
    """
    protein_mass = 0

    aa = {'A': 71.03711, 'R': 156.10111, 'N': 114.04293,
          'D': 115.02694, 'C': 103.00919, 'E': 129.04259,
          'Q': 128.05858, 'G': 57.02146, 'H': 137.05891,
          'I': 113.08406, 'L': 113.08406, 'K': 128.09496,
          'M': 131.04049, 'F': 147.06841, 'P': 97.05276,
          'S': 87.03203, 'T': 101.04768, 'W': 186.07931,
          'Y': 163.06333, 'V': 99.06841, 'X': 0.00000}

    for i in amino_seq:
        protein_mass += aa[i]

    return round(protein_mass, 3)


def find_motifs(sequence, subsequence):
    """
    returns all locations of a subsequence within a sequence
    """
    mot = sequence.find(subsequence, 0)
    motifs = []
    while mot >= 0:
        motifs.append(mot)
        mot = sequence.find(subsequence, mot + 1)
    motifs = [x + 1 for x in motifs]
    return motifs
