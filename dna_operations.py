from nucleotide import Nucleotide


def nucleotide_list_from_nucleotide_string(nuc_str: str):
    nuc_list = list()
    for c in nuc_str:
        nuc_list.append(Nucleotide.from_letter(c))

    return nuc_list


initial_vector_h0 = [
    nucleotide_list_from_nucleotide_string("CGGGAAGCTGCGCGCTTTATGTTATAGCAAGA"),
    nucleotide_list_from_nucleotide_string("GTGTCGCTGGTGGACCGACATAGGGGCTATGT"),
    nucleotide_list_from_nucleotide_string("ATTACGTGTTATCTAGTTTGGCCATTGAAGGT"),
    nucleotide_list_from_nucleotide_string("GGCCCATTTTCCATGGCCTTACTCATCGTTAC"),
    nucleotide_list_from_nucleotide_string("CCACAATGCCAGCTTTGGTCTGCGGAAGTCAC"),
    nucleotide_list_from_nucleotide_string("GCGTAACCCGGAGATAAGGTATTGCGTAACTT"),
    nucleotide_list_from_nucleotide_string("ACTTGAATTCGCGGGTTTGTCAACGTTCCGGT"),
    nucleotide_list_from_nucleotide_string("CCGTTGAATATCACGCACATCTTGAGACCTGC")
]


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
    word1.reverse()
    word2.reverse()
    result = list()
    z0, e = word1[0] + word2[0]
    result.append(z0)
    for i in range(1, 32):
        x, ex = word1[i] + e
        zi, ey = x + word2[i]
        e, a = ex + ey
        result.append(zi)
    result.reverse()
    return result


print(nucleotide_addition(
    nucleotide_list_from_nucleotide_string("TCTTTTCAGTACAATTTGCATAACGTGTGGAA"),
    nucleotide_list_from_nucleotide_string("TGATAGCTATTCGATTTACTAAGCATATGTGA")
))

print(dnamaj(
    nucleotide_list_from_nucleotide_string("AGC"),
    nucleotide_list_from_nucleotide_string("ATC"),
    nucleotide_list_from_nucleotide_string("AGG"),
))
