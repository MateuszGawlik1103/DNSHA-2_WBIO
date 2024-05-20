from data_conversion import nuc_string_to_nuc_list, add_padding, text_to_nuc_list, binary_to_nuc_list, add_padding2
from nucleotide import Nucleotide
from functools import cache
from hashlib import sha512

k_values = [
    "CAAGGAGGAGTTGCGATCCTAGGAGGTGAGAG", "CTACATCTCACAGCACAGATTGTTCGCCTATC", "GTCCTAAATTGTTATTTGTACATCATGTAGTT",
    "TGGCGTCCTCGTGGCCGAACGAGCTCGTGTTA", "ATGCCCCGTAAGCCGTTTATCAGAGTCCATGA", "CCGCTTACACACTTACGTCGAACCTCAAACGC",
    "GCAGATTTGAAGGGCAGGTTACGCCATTGCGT", "GGGTACTACCTGTCCCTCGGCGTCGAACACGA", "TCGAAACTGGGGGCGAGGATAAATAAAGCAAG",
    "ACAGGAATCCGTAAACCACCCTAACGTTGTTG", "AGCAATACGACCGTTGCATGTGCAGTAGGATA", "CCCCAATACTTCTAATTCCCTTTTGTCATGAG",
    "CTAGGTTGCCTCCTCATTAGCTGTGAGCCGTT", "GAAATCTGGTACTTTGATGTACCGGCCGGTAC", "GCGTTCTAAACGGGCTAGCCTACTACAGATCC",
    "TAACGCGTTTACCTCATATTCGGCAGCGGCCA", "TGCAGCGTCGGCTAACGCTGTTACCAGGTCAG", "TGTTGTTGCACTGACGATGACATTAGCCTGAT",
    "AATTTAACGCTCTACGGAGTGATATCCCGTCC", "AGCAAATAGGACTATACTCTGGTAGCTACGCC", "AGTCTGGCAGTACGTTCCGCAGGTAAAGCTCC",
    "CAGGCTCAGACAGGGGCGTGGGCGTGCAGAAT", "CCTAGTAAGGGCTCTAGTTCCAACTTGTTCCA", "CTCGTTGCGAGATCGGGAATACACCCATGTCC",
    "GCGAATTGCCACCCAGTGTGCGCGTCTTGGGT", "GGGAATACTACGCGTCAGTCGTCAATAGACAA", "GTAAAAATAGCTTAGAGCGATTGTAGACATTT",
    "GTTTCCGCCTTTTACTGTTGTGTTAATGTGCA", "TACGTGAAAAGTTTATATTCGGGAGATTTAAG", "TCCCGGCTGCACCACTGCATAAGGGGCTAGCC",
    "AACGTAGGCGATCCACTGAAAAATGAAGCGTT", "ACCAAGGCAGGCCGCTAAGGAATGCGTGCTAA", "AGCTGTCTAAGGGACCCACGTCAGAGTTTTTA",
    "AGTGACGTAGACATGACCTAAGCGTAGCAGCG", "CATCAGTACGTCTTTACCGGTACAAGGGTGTC", "CCATATGAAATCACATGCTCGCCCGTATTCTT",
    "CGCCAAGGCTATCCCAGAGTGGTTCGATTCTG", "CTCGCGGGAAGGGTGTATTACTCTGTAGGGGA", "GAACTAAGTAGCAGTGCACTTGTCGGTGTGCG",
    "GCAGCTAGAGTAGACCACCAGAAGATCCATGT", "GGAGGTTTTGGAGGACCATATTACAAATCGCA", "GGGAACGGCGCGCAGTGTTACAAGATAAAAAC",
    "TAAGCAGTGAGTCTAATCAATTGAGCCTGCAC", "TACTCGTACCACGGATAACGCCCAGTTGATAA", "TCACGCAGTGGAACGCTCCGTGTTCCAGACGA",
    "TCCGGCGCAACGAGCACCCCCGCCGGGCACAA", "TTCAAATGATCCGACCCCCTCTACAGAAAGGG", "ACAACGGGGGAACTAAATAGGTGTTCACGTGA",
    "ACGCGGCATAACACCGGTGATCAGTCAATAGA", "ACTGATCTCGTAAAGACCACCAACGGGTCCAT", "AGCTCAGACTCTCATATCTTGATGTGGTGCGC",
    "ATCAGTAAGTTAGTCCTGACGCGTCAGAGGGA", "ATGCACTAAATAGTATTACCTAGCCCGGCGAT", "CATGTCGAGGGGCAGGTGATCAACGAGGTAGT",
    "CCGTGCTATAGGCATTCTCTCGATTGATCTAT", "CGGAAGTGCGTTTTATTCCGGTAGGTGAGGAT", "CTCAGATTGAAGTGTGCCTCTGTTGTAGTTTA",
    "CTGAGGCCCGATCGTTCAATACCTAGTTCGAA", "GACATAGACTGAACCAGGACTTAAGGGTCTAG", "GATATACTAAAGAAGAACGGCGCAATGCTGTA",
    "GCAAGTTGTTTTTTGGAGATCGATACTGAGGA", "GGCACCAACGTATGGTTCTGGAAGGTTCTGGC", "GTTGTTGCGGATTTCTGTAGTACGCTGCACCC",
    "TACGCTACCTGATTAGTGATCTAGCCATAGGT", "TAGGAGCTATTGTATGTGGGAGCGCGACGCTA", "TCACGACGGTGATACTAGACTAAATAAGAACT",
    "TGGGTCGGCTTCTCCGTATCTGAATGGTACTG", "TTCCCTTCCATTCTTTTGTGCGTGTCACCTGA", "AACGTTAACGCTGGGGCTAGACCTCGTTGTGG",
    "AAGGCGATCTTCTACCGGAGTAGAGCGAGGCG", "ACACATTTGCGAAACAGTTGTTGCAATCGGTG", "ACGTCTACAAGTATCCACATACTACACTACGT",
    "AGGATCGTCTCTTTCCAGATAACACTTCGACA", "ATAGTAGGGGGTCTGTCAAATACTAGCAGCAT", "ATTAGCTGGTTGAAGGACCCTAGCGTTGGTTA",
    "CAATACTCCGCTTACAGCTAACAAAATCCATA", "CATATACCTCCAGTTGTAGTATTGCAAGGTCG", "CCGCCTTTAGGCGCTATTTACGCCCTTGAGGG",
    "CCTTTAGTCGTTGGGTATGGTCCGTTGGTGTA", "CGTACACAACGCGATACAGGCACTCCGAACCT"
]

initial_vector_h0 = [
    nuc_string_to_nuc_list("CGGGAAGCTGCGCGCTTTATGTTATAGCAAGA"),
    nuc_string_to_nuc_list("GTGTCGCTGGTGGACCGACATAGGGGCTATGT"),
    nuc_string_to_nuc_list("ATTACGTGTTATCTAGTTTGGCCATTGAAGGT"),
    nuc_string_to_nuc_list("GGCCCATTTTCCATGGCCTTACTCATCGTTAC"),
    nuc_string_to_nuc_list("CCACAATGCCAGCTTTGGTCTGCGGAAGTCAC"),
    nuc_string_to_nuc_list("GCGTAACCCGGAGATAAGGTATTGCGTAACTT"),
    nuc_string_to_nuc_list("ACTTGAATTCGCGGGTTTGTCAACGTTCCGGT"),
    nuc_string_to_nuc_list("CCGTTGAATATCACGCACATCTTGAGACCTGC")
]


# input artificial DNA (list); output - one bit right shift DNA (list)
def rsob(dna):
    m = len(dna)
    b1 = [Nucleotide.C] * m
    b2 = [Nucleotide.G] * m
    b3 = [Nucleotide.A] + dna[:-1]
    b4 = [x & y for x, y in zip(b1, b3)]
    b5 = [x & y for x, y in zip(dna, b2)]
    return [x.nucleotide_operation(y) for x, y in zip(b4, b5)]


# input: dna - string of nucleotides, k - rotation (in bits)
def r_shift(dna, k):
    if k == 1:
        return rsob(dna)
    elif k % 2 == 0:
        return [Nucleotide.A] * (k // 2) + dna[:-(k // 2)]
    else:
        alfa = [Nucleotide.A] * ((k - 1) // 2) + dna[:-((k - 1) // 2)]
        return rsob(alfa)


def lsob(dna):
    m = len(dna)
    b1 = [Nucleotide.C] * m
    b2 = [Nucleotide.G] * m
    b3 = dna[1:] + [Nucleotide.A]
    b4 = [x & y for x, y in zip(b2, b3)]
    b5 = [x & y for x, y in zip(dna, b1)]
    return [x.nucleotide_operation(y) for x, y in zip(b4, b5)]


# left shift:
def l_shift(dna, k):
    if k % 2 == 0:
        return dna[(k // 2):] + [Nucleotide.A] * (k // 2)
    else:
        alfa = dna[int((k - 1) // 2):] + [Nucleotide.A] * ((k - 1) // 2)
        return lsob(alfa)


# right rotation
def r_rotate(dna, k):
    mm = len(dna)
    b1 = r_shift(dna, k)
    b2 = l_shift(dna, (2 * mm) - k)
    return [x | y for x, y in zip(b1, b2)]


def nuc_lists_xor(*dnas):
    length = len(dnas[0])
    # for dna in dnas[1:]:
    #     if len(dna) != length:
    #         raise ValueError("All DNA sequences must have the same length")
    result = []
    for i in range(length):
        xor_result = dnas[0][i]
        for dna in dnas[1:]:
            xor_result ^= dna[i]
        result.append(xor_result)
    return result


def sum0(alpha):
    return nuc_lists_xor(r_rotate(alpha, 28), r_rotate(alpha, 34), r_rotate(alpha, 39))


def sum1(alpha):
    return nuc_lists_xor(r_rotate(alpha, 14), r_rotate(alpha, 18), r_rotate(alpha, 41))


def dnach(r1: list[Nucleotide], r2: list[Nucleotide], r3: list[Nucleotide]):
    r1_and_r2 = [x & y for (x, y) in zip(r1, r2)]
    not_r1 = [~x for x in r1]
    not_r1_and_r3 = [x & y for (x, y) in zip(not_r1, r3)]
    return [x ^ y for (x, y) in zip(r1_and_r2, not_r1_and_r3)]


def dnamaj(r1: list[Nucleotide], r2: list[Nucleotide], r3: list[Nucleotide]):
    r1_and_r2 = [x & y for (x, y) in zip(r1, r2)]
    r1_and_r3 = [x & y for (x, y) in zip(r1, r3)]
    r2_and_r3 = [x & y for (x, y) in zip(r2, r3)]
    return [x ^ y ^ z for (x, y, z) in zip(r1_and_r2, r1_and_r3, r2_and_r3)]


def nucleotide_addition(word1: list[Nucleotide], word2: list[Nucleotide]):
    word1_copy = word1[::-1]
    word2_copy = word2[::-1]
    result = []
    z0, e = word1_copy[0] + word2_copy[0]
    result.append(z0)
    for i in range(1, 32):
        x, ex = word1_copy[i] + e
        zi, ey = x + word2_copy[i]
        e, a = ex + ey
        result.append(zi)
    result.reverse()
    return result


# input: one 512-nucleotide block M(i)
# output: Wj DNA sequence of 32 nucleotides
@cache
def compute_wj(mi, j):
    if 0 <= j <= 15:
        return tuple(list(mi)[32 * j:32 * (j + 1)])
    else:
        w15 = list(compute_wj(mi, j - 15))
        w2 = list(compute_wj(mi, j - 2))
        sigma0 = nuc_lists_xor(r_rotate(w15, 1), r_rotate(w15, 8), r_shift(w15, 7))
        sigma1 = nuc_lists_xor(r_rotate(w2, 19), r_rotate(w2, 61), r_shift(w2, 6))
        wj = sigma0
        temp = [list(compute_wj(mi, j - 7)), sigma1, list(compute_wj(mi, j - 16))]
        for t in temp:
            wj = nucleotide_addition(wj, t)

        return tuple(wj)


def dnsha512_hash(blocks):
    hi = initial_vector_h0

    for i, block in enumerate(blocks):
        print(i)
        r1, r2, r3, r4, r5, r6, r7, r8 = hi

        for j in range(80):
            k = nuc_string_to_nuc_list(k_values[j])
            t1 = r8
            temp = [sum1(r5), dnach(r5, r6, r7), k, list(compute_wj(tuple(block), j))]
            for t in temp:
                t1 = nucleotide_addition(t1, t)
            t2 = nucleotide_addition(sum0(r1), dnamaj(r1, r2, r3))
            r8, r7, r6, r5, r4, r3, r2 = r7, r6, r5, nucleotide_addition(r4, t1), r3, r2, r1
            r1 = nucleotide_addition(t1, t2)

        hi = [
            nucleotide_addition(r1, hi[0]), nucleotide_addition(r2, hi[1]),
            nucleotide_addition(r3, hi[2]), nucleotide_addition(r4, hi[3]),
            nucleotide_addition(r5, hi[4]), nucleotide_addition(r6, hi[5]),
            nucleotide_addition(r7, hi[6]), nucleotide_addition(r8, hi[7])
        ]

    return hi


def hash_text(text: str):
    print("text to hash " + text)
    print("nucleotide sequence before hash " + add_padding(text_to_nuc_list(text)))
    return dnsha512_hash(add_padding(text_to_nuc_list(text)))


def hash_photo(path: str):
    with open(path, 'rb') as file:
        bits = file.read()
    bit_string = ''.join(format(byte, '08b') for byte in bits)
    return dnsha512_hash(add_padding2(binary_to_nuc_list(bit_string)))



