from nucleotide import Nucleotide


def nuc_string_to_nuc_list(nuc_str: str):
    nuc_list = list()
    for c in nuc_str:
        nuc_list.append(Nucleotide.from_letter(c))
    return nuc_list


# dna to binary
def nuc_string_to_binary(dna):
    nucleotide_to_binary = {'A': '00', 'C': '01', 'G': '10', 'T': '11'}
    binary_sequence = ""
    for nucleotide in dna:
        if nucleotide.name not in nucleotide_to_binary:
            raise ValueError("invalid nucleotide: {}")
        binary_sequence += nucleotide_to_binary[nucleotide.name]
    return binary_sequence


# binary to dna
def binary_to_nuc_list(binary_sequence):
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


def text_to_bits_list(input_string):
    result_bits = []
    for char in input_string:
        ascii_code = ord(char)
        binary_string = bin(ascii_code)[2:].zfill(8)
        result_bits.extend(map(int, binary_string))
    return result_bits


def bits_list_to_text(bit_sequence):
    result_string = ""
    for i in range(0, len(bit_sequence), 8):
        bits = bit_sequence[i:i + 8]
        ascii_code = int("".join(map(str, bits)), 2)
        result_string += chr(ascii_code)
    return result_string


# input text; output - list of nucleotides
def text_to_nuc_list(text):
    dna_sequence: list[Nucleotide] = list()
    for char in text:
        for i in range(3, -1, -1):
            two_bits = (ord(char) & (2 ** (2 * i) + 2 ** (2 * i + 1))) >> (2 * i)
            dna_sequence.append(Nucleotide(two_bits))
    return dna_sequence


# input - DNA artificial sequence; output - list of bits
def nuc_string_to_bits(nuc_list):
    res = []
    for letter in nuc_list:
        binary_nuc = bin(Nucleotide.from_letter(letter).value)[2:].zfill(2)
        res.extend(map(int, binary_nuc))
    return res


def nuc_list_to_text(nuc_list):
    result_string = ""
    for i in range(0, len(nuc_list), 4):
        nucleotides = nuc_list[i:i + 4]
        ascii_code = ((nucleotides[0].value << 6) + (nucleotides[1].value << 4)
                      + (nucleotides[2].value << 6) + nucleotides[3].value)
        result_string += chr(ascii_code)
    return result_string


def add_padding(nuc_list: [Nucleotide]) -> [Nucleotide]:
    length_string = bin(len(nuc_list) * 2)[2:]
    length_nuc_list = binary_to_nuc_list(length_string) if len(length_string) % 2 == 0 else binary_to_nuc_list(
        '0' + length_string)
    nuc_list += [Nucleotide.G]
    k = 512 - ((len(nuc_list) + 64) % 512)
    nuc_list += [Nucleotide.A] * k
    nuc_list += [Nucleotide.A] * (64 - len(length_nuc_list))
    nuc_list += length_nuc_list
    return nuc_list

def add_padding2(nuc_list: [Nucleotide]) -> [Nucleotide]:
    length_string = bin(len(nuc_list) * 2)[2:]
    length_nuc_list = binary_to_nuc_list(length_string) if len(length_string) % 2 == 0 else binary_to_nuc_list(
        '0' + length_string)
    nuc_list += [Nucleotide.G]
    k = 512 - ((len(nuc_list) + 64) % 512)
    nuc_list += [Nucleotide.A] * k
    nuc_list += [Nucleotide.A] * (64 - len(length_nuc_list))
    nuc_list += length_nuc_list
    return [nuc_list[i:i + 512] for i in range(0,len(nuc_list),512)]


