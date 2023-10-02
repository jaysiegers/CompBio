def validate_dna(dna_seq):
    #convert uppercase
    seqm = dna_seq.upper()
    print(seqm)

    valid = seqm.count("A") + seqm.count("T") + seqm.count("C") + seqm.count("G")

    if valid == len(seqm):
        print("DNA Valid")
        return True
    else:
        print("DNA Invalid")
        return False

def calcfreq(seq):
    dictionary = {}
    for s in seq.upper():
        if s in dictionary:
            dictionary[s] += 1
        else:
            dictionary[s] = 1
    return dictionary

def count_percentage(freq_data):
    dna_length = sum(freq_data.values())
    print("DNA Length = ", str(dna_length))
    g, c = 0, 0
    if "G" in freq_data:
        g_percentage = (freq_data.get("G")/dna_length) * 100
    if "C" in freq_data:
        c_percentage = (freq_data.get("C")/dna_length) * 100
    return g_percentage, c_percentage

def convert_dna(dna_seq):
    complement = ""
    complement_pairs = {"A": "T", "T": "A", "C": "G", "G": "C"}
    for i in dna_seq.upper():
        complement += complement_pairs[i]
    return complement

def convert_complement(complement_seq):
    mrna = ""
    pairs = {"T": "U"}
    for i in complement_seq:
        if i in pairs:
            mrna += pairs[i]
        else:
            mrna += i
    return mrna

def convertAminoAcid(mrna):
    amino_acid = ''
    genetic_code = {
        'UUU': 'Phe (F)', 'UUC': 'Phe (F)', 'UUA': 'Leu (L)', 'UUG': 'Leu (L)',
        'UCU': 'Ser (S)', 'UCC': 'Ser (S)', 'UCA': 'Ser (S)', 'UCG': 'Ser (S)',
        'UAU': 'Tyr (Y)', 'UAC': 'Tyr (Y)', 'UAA': 'STOP', 'UAG': 'STOP',
        'UGU': 'Cys (C)', 'UGC': 'Cys (C)', 'UGA': 'STOP', 'UGG': 'Trp (W)',
        'CUU': 'Leu (L)', 'CUC': 'Leu (L)', 'CUA': 'Leu (L)', 'CUG': 'Leu (L)',
        'CCU': 'Pro (P)', 'CCC': 'Pro (P)', 'CCA': 'Pro (P)', 'CCG': 'Pro (P)',
        'CAU': 'His (H)', 'CAC': 'His (H)', 'CAA': 'Gln (Q)', 'CAG': 'Gln (Q)',
        'CGU': 'Arg (R)', 'CGC': 'Arg (R)', 'CGA': 'Arg (R)', 'CGG': 'Arg (R)',
        'AUU': 'Ile (I)', 'AUC': 'Ile (I)', 'AUA': 'Ile (I)', 'AUG': 'Met (M)',
        'ACU': 'Thr (T)', 'ACC': 'Thr (T)', 'ACA': 'Thr (T)', 'ACG': 'Thr (T)',
        'AAU': 'Asn (N)', 'AAC': 'Asn (N)', 'AAA': 'Lys (K)', 'AAG': 'Lys (K)',
        'AGU': 'Ser (S)', 'AGC': 'Ser (S)', 'AGA': 'Arg (R)', 'AGG': 'Arg (R)',
        'GUU': 'Val (V)', 'GUC': 'Val (V)', 'GUA': 'Val (V)', 'GUG': 'Val (V)',
        'GCU': 'Ala (A)', 'GCC': 'Ala (A)', 'GCA': 'Ala (A)', 'GCG': 'Ala (A)',
        'GAU': 'Asp (D)', 'GAC': 'Asp (D)', 'GAA': 'Glu (E)', 'GAG': 'Glu (E)',
        'GGU': 'Gly (G)', 'GGC': 'Gly (G)', 'GGA': 'Gly (G)', 'GGG': 'Gly (G)'
    }
    codons = [mrna[i:i + 3] for i in range(0, len(mrna), 3)]

    for codon in codons:
        if codon in genetic_code:
            aminoacid = genetic_code[codon]
            if aminoacid == 'STOP':
                aminoacid = "*"
                break
            amino_acid += aminoacid
            amino_acid += " - "

    return amino_acid

dna = "atagcctga"
validate_result = validate_dna(dna)

if (validate_result == True):
    freq_data = calcfreq(dna)
    print(freq_data)

    result_g, result_c = count_percentage(freq_data)
    print("G = ", round(result_g), "%")
    print("C = ", round(result_c), "%")

    complement_seq = convert_dna(dna)
    print("Complement Sequence:", complement_seq)

    mrna_seq = convert_complement(complement_seq)
    print("mRNA Sequence:", mrna_seq)

    aminoacid_seq = convertAminoAcid(mrna_seq)
    print("Amino Acid Sequence:", aminoacid_seq)

codon_map = {
    'A': ['GCU', 'GCC', 'GCA', 'GCG'],
    'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'N': ['AAU', 'AAC'],
    'D': ['GAU', 'GAC'],
    'C': ['UGU', 'UGC'],
    'Q': ['CAA', 'CAG'],
    'E': ['GAA', 'GAG'],
    'G': ['GGU', 'GGC', 'GGA', 'GGG'],
    'H': ['CAU', 'CAC'],
    'I': ['AUU', 'AUC', 'AUA'],
    'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
    'K': ['AAA', 'AAG'],
    'M': ['AUG'],
    'F': ['UUU', 'UUC'],
    'P': ['CCU', 'CCC', 'CCA', 'CCG'],
    'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'],
    'T': ['ACU', 'ACC', 'ACA', 'ACG'],
    'W': ['UGG'],
    'Y': ['UAU', 'UAC'],
    'V': ['GUU', 'GUC', 'GUA', 'GUG'],
    '*': ['UAA', 'UAG', 'UGA']
}

def translate_mrna(aminoacid_seq):
    seq = ''
    for char in aminoacid_seq:
        if char.isalpha():
            seq += char.upper()
    rna_seq = ""
    for i in seq:
        rna_seq += codon_map[i][0]
    return rna_seq

def codon_freq(seq):
    dictionary = {}

    for i in range(0,len(seq),3):
        codon = seq[i:i +3].upper()
        if codon in dictionary:
            dictionary[codon]+=1
        else:
            dictionary[codon] = 1

    return dictionary

codon_seq = "n_R-L-A-N"

convert_mrna = translate_mrna(codon_seq)
print("mRNA:", convert_mrna)

codonfreq = codon_freq(convert_mrna)
print("Codon Frequency:", codonfreq)
