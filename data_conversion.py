from nucleotide import Nucleotide


def nucleotide_list_from_nucleotide_string(nuc_str: str):
    nuc_list = list()
    for c in nuc_str:
        nuc_list.append(Nucleotide.from_letter(c))

    return nuc_list


# dna to binary
def nuc_to_binary(dna):
    nucleotide_to_binary = {'A': '00', 'C': '01', 'G': '10', 'T': '11'}
    binary_sequence = ""
    for nucleotide in dna:
        if nucleotide.name not in nucleotide_to_binary:
            raise ValueError("invalid nucleotide: {}")
        binary_sequence += nucleotide_to_binary[nucleotide.name]
    return binary_sequence


# binary to dna
def binary_to_nuc(binary_sequence):
    binary_to_nucleotide = {'00': Nucleotide.A, '01': Nucleotide.C,
                            '10': Nucleotide.G, '11': Nucleotide.T}
    if len(binary_sequence) % 2 != 0:
        raise ValueError("invalid binary sequence length")
    nucleotide_sequence = []
    for i in range(0, len(binary_sequence), 2):
        nucleotide_binary = binary_sequence[i:i + 2]
        if nucleotide_binary not in binary_to_nucleotide:
            raise ValueError("invalid binary representation")
        nucleotide_sequence.append(binary_to_nucleotide[nucleotide_binary])
    return nucleotide_sequence


def string_to_bits(input_string):
    result_bits = []
    for char in input_string:
        ascii_code = ord(char)
        binary_string = bin(ascii_code)[2:].zfill(8)
        result_bits.extend(map(int, binary_string))
    return result_bits


def bits_to_string(bit_sequence):
    result_string = ""
    for i in range(0, len(bit_sequence), 8):
        bits = bit_sequence[i:i + 8]
        ascii_code = int("".join(map(str, bits)), 2)
        result_string += chr(ascii_code)
    return result_string


# input text; output - list of nucleotides
def dna_encoding(text):
    dna_sequence: list[Nucleotide] = list()
    for char in text:
        for i in range(3, -1, -1):
            two_bits = (ord(char) & (2 ** (2 * i) + 2 ** (2 * i + 1))) >> (2 * i)
            dna_sequence.append(Nucleotide(two_bits))
    return dna_sequence


# input - DNA artificial sequence; output - list of bits
def dna_decoding(alfa):
    res = []
    for letter in alfa:
        binary_nuc = bin(Nucleotide.from_letter(letter).value)[2:].zfill(2)
        res.extend(map(int, binary_nuc))
    return res
