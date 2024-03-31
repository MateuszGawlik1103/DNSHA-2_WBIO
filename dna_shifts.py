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
    elif (l1.name,l0.name in rsob_table):
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
    return lsob_nuc(B4,B5)

def lsob_nuc(nucl1, nucl2):
    res = []
    for i in range(len(nucl1)):
        res.append(lsob_operation(nucl1[i], nucl2[i]))
    return res

# it uses the same nucleotide operation table as right shift
def lsob_operation(l0,l1):
    if (l0.name,l1.name) in rsob_table:
        return rsob_table[(l0.name,l1.name)]
    elif (l1.name,l0.name in rsob_table):
        return rsob_table[(l1.name,l0.name)]
    else: return None



print("right shift check: (should be AATA)")
print(r_shift([Nucleotide.T, Nucleotide.A, Nucleotide.G, Nucleotide.C],4))

print("left shift check: (should be GCAA)")
print(l_shift([Nucleotide.T, Nucleotide.A, Nucleotide.G, Nucleotide.C],4))

# for not even k, both shifts also work (checked)