from nucleotide import Nucleotide

rsob_table = {
    ('A', 'A'): Nucleotide.A,
    ('A', 'C'): Nucleotide.G,
    ('A', 'G'): Nucleotide.C,
    ('G', 'C'): Nucleotide.T
}

def rsob_operation(l0,l1):
    if (l0.name,l1.name) in rsob_table:
        return rsob_table[(l0.name,l1.name)]
    elif (l1.name,l0.name) in rsob_table:
        return rsob_table[(l1.name,l0.name)]
    else: return None


#input two lists of nucleotides; output - list of nucleotides after rsob operation
def rsob_nuc(nucl1, nucl2):
    res = []
    for i in range(len(nucl1)):
        res.append(rsob_operation(nucl1[i], nucl2[i]))
    return res

print(rsob_nuc([Nucleotide.A, Nucleotide.A], [Nucleotide.C,Nucleotide.G]))



# input artificial DNA (list); output - one bit right shift DNA (list)
def rsob(dna):
    m = len(dna)
    B1 = [Nucleotide.C] * m
    B2 = [Nucleotide.G] * m
    B3 = [Nucleotide.A] + dna[:-1]
    B4 = [x & y for x, y in zip(B1, B3)]
    B5 = [x & y for x, y in zip(dna, B2)]
    return rsob_nuc(B4,B5)

print(rsob([Nucleotide.T, Nucleotide.A, Nucleotide.C]))

# input: dna - string of nucleotides, k - rotation (in bits)
def r_shift(dna, k):
    if k % 2 == 0:
        return [Nucleotide.A]*int(k/2) + dna[:-int(k/2)]
    else:
        alfa = [Nucleotide.A]*int((k-1)/2) + dna[:-int((k-1)/2)]
        return rsob(alfa)

print(r_shift([Nucleotide.G, Nucleotide.C, Nucleotide.T, Nucleotide.G, Nucleotide.A, Nucleotide.T],3))




# left shift:
def l_shift(dna, k):
    if k % 2 == 0:
        return dna[int(k/2):] + [Nucleotide.A]*int(k/2)
    else:
        alfa = dna[int((k-1)/2):] + [Nucleotide.A]*int((k-1)/2)
        return lsob(alfa)

def lsob(dna):
    m = len(dna)
    B1 = [Nucleotide.C] * m
    B2 = [Nucleotide.G] * m
    B3 = dna[1:] + [Nucleotide.A]
    B4 = [x & y for x, y in zip(B2, B3)]
    B5 = [x & y for x, y in zip(dna, B1)]
    # it uses the same nucleotide operation table as right shift
    return rsob_nuc(B4, B5)


# right rotation
def r_rotate(dna, k):
    m = len(dna)
    B1 = r_shift(dna, k)
    B2 = l_shift(dna, (2*m)-k)
    return nuc_or(B1, B2)


# the nucleotide operation to imitate the bitwise OR
def nuc_or(dna1, dna2):
    binary1 = nuc_to_binary(dna1)
    binary2 = nuc_to_binary(dna2)
    if len(binary1) != len(binary2):
        raise ValueError("not the same length")
    result = ""
    for bit1, bit2 in zip(binary1, binary2):
        if bit1 == '1' or bit2 == '1':
            result += '1'
        else:
            result += '0'
    return binary_to_nuc(result)


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
        nucleotide_binary = binary_sequence[i:i+2]
        if nucleotide_binary not in binary_to_nucleotide:
            raise ValueError("invalid binary representation")
        nucleotide_sequence.append(binary_to_nucleotide[nucleotide_binary])
    return nucleotide_sequence



print("right shift check: (should be AATA)")
print(r_shift([Nucleotide.T, Nucleotide.A, Nucleotide.G, Nucleotide.C],4))

print("left shift check: (should be GCAA)")
print(l_shift([Nucleotide.T, Nucleotide.A, Nucleotide.G, Nucleotide.C],4))

print("right rotate check: even k (should be GTA)")
print(r_rotate([Nucleotide.A, Nucleotide.G, Nucleotide.T],4))

print("right rotate check 2: not even k (should be TAA)")
print(r_rotate([Nucleotide.A, Nucleotide.C, Nucleotide.G],3))