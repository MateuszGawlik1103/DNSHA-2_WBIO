#AND
and_table = {
    ('A', 'A'): 'A',
    ('A', 'C'): 'A',
    ('A', 'G'): 'A',
    ('A', 'T'): 'A',
    ('C', 'C'): 'C',
    ('C', 'G'): 'A',
    ('C', 'T'): 'C',
    ('G', 'G'): 'G',
    ('G', 'T'): 'G',
    ('T', 'T'): 'T', 
}

or_table = {
    ('A', 'A'): 'A',
    ('A', 'C'): 'C',
    ('A', 'G'): 'G',
    ('A', 'T'): 'T',
    ('C', 'C'): 'C',
    ('C', 'G'): 'T',
    ('C', 'T'): 'T',
    ('G', 'G'): 'G',
    ('G', 'T'): 'T',
    ('T', 'T'): 'T', 
}
xor_table = {
    ('A', 'A'): 'A',
    ('A', 'C'): 'C',
    ('A', 'G'): 'G',
    ('A', 'T'): 'T',
    ('C', 'C'): 'A',
    ('C', 'G'): 'T',
    ('C', 'T'): 'G',
    ('G', 'G'): 'A',
    ('G', 'T'): 'C',
    ('T', 'T'): 'A', 
}

#input two nucleotides
def and_nuc(l0,l1):
    if (l0,l1) in and_table:
        return and_table[(l0,l1)]
    elif (l1,l0 in and_table):
        return and_table[(l1,l0)]
    else: return None

#input two strings of nucleotides
def and_nuc_strings(nucl1, nucl2):
    res = ""
    for i in range(len(nucl1)):
        res += and_nuc(nucl1[i],nucl2[i])
    return res






#input two nucleotides
def or_nuc(l0,l1):
    if (l0,l1) in or_table:
        return or_table[(l0,l1)]
    elif (l1,l0 in or_table):
        return or_table[(l1,l0)]
    else: return None

#input two strings of nucleotides
def or_nuc_strings(nucl1, nucl2):
    res = ""
    for i in range(len(nucl1)):
        res += or_nuc(nucl1[i],nucl2[i])
    return res
    





#input two nucleotides
def xor_nuc(l0,l1):
    if (l0,l1) in xor_table:
        return xor_table[(l0,l1)]
    elif (l1,l0 in xor_table):
        return xor_table[(l1,l0)]
    else: return None


#input two strings of nucleotides
def xor_nuc_strings(nucl1, nucl2):
    res = ""
    for i in range(len(nucl1)):
        res += xor_nuc(nucl1[i],nucl2[i])
    return res





