import nucleotide_operations, coding

rsob_table = {
    ('A', 'A'): 'A',
    ('A', 'C'): 'G',
    ('A', 'G'): 'C',
    ('G', 'C'): 'T'
}

def rsob_operation(l0,l1):
    if (l0,l1) in rsob_table:
        return rsob_table[(l0,l1)]
    elif (l1,l0 in rsob_table):
        return rsob_table[(l1,l0)]
    else: return None


#input two strings of nucleotides
def rsob_nuc_strings(nucl1, nucl2):
    res = ""
    for i in range(len(nucl1)):
        res += rsob_operation(nucl1[i],nucl2[i])
    return res
    


# input artificial DNA (String); output - one bit right shift DNA (String)
def rsob(dna):
    m = len(dna)
    B1 = 'C'*m
    B2 = 'G'*m
    B3 = 'A' + dna[:-1]
    B4 = nucleotide_operations.and_nuc_strings(B1,B3)
    B5 = nucleotide_operations.and_nuc_strings(dna,B2)
    return rsob_nuc_strings(B4,B5)

# input: dna - string of nucleotides, k - rotation (in bits)
def r_shift(dna, k):
    if k % 2 == 0:
        return 'A'*int(k/2) + dna[:-int(k/2)]
    else:
        alfa = 'A'*int((k-1)/2) + dna[:-int((k-1)/2)]
        return rsob(alfa)
    
print(r_shift("GCTGAT",3))


