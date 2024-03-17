from nucleotide import Nucleotide

def string_to_bits(input_string):
    result_bits = []
    for char in input_string:
        ascii_code = ord(char)
        binary_string = bin(ascii_code)[2:].zfill(8)
        result_bits.extend(map(int, binary_string))
    return result_bits


input_string = "BOB"
result = string_to_bits(input_string)
print(result)


def bits_to_string(bit_sequence):
    result_string = ""
    for i in range(0, len(bit_sequence), 8):
        bits = bit_sequence[i:i+8]
        ascii_code = int("".join(map(str, bits)), 2)
        result_string += chr(ascii_code)
    return result_string


bit_sequence = string_to_bits("BOB") 
result = bits_to_string(bit_sequence)
print(result)



#input text; output - list of nucleotides
def dna_encoding(text):
    dna_sequence: list[Nucleotide] = list()
    for char in text:
        for i in range(3, -1, -1):
            two_bits = (ord(char) & (2 ** (2 * i) + 2 ** (2 * i + 1))) >> (2 * i)
            dna_sequence.append(Nucleotide(two_bits))
    return dna_sequence

print(str(dna_encoding("BOB")))

#input - DNA artificial sequence; output - list of bits
def dna_decoding(alfa):
    res = []
    for letter in alfa:
        binary_nuc = bin(Nucleotide.from_letter(letter).value)[2:].zfill(2)
        res.extend(map(int,binary_nuc))
    return res
print(dna_decoding("GCT"))




