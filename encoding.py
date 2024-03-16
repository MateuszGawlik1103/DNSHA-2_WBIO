def word_to_binary(word):
    binary_string = ""
    for char in word:
        binary_string += bin(ord(char))[2:].zfill(8)  
    return binary_string

def binary_to_word(binary_string):
    word = ""
    for i in range(0,int(len(binary_string)/8)):
        letter = chr(int(binary_string[0+i*8 : 8+i*8],2))
        word += letter
    return word



word = "BOB"
binary_representation = word_to_binary(word)
#print(binary_representation)

# input 2-bit data, output DNA artificial sequence
def transformation(b1,b0):
    bin_in = int(str(b1)+str(b0),2)

    if bin_in == 0b00:
        return 'A'
    elif bin_in == 0b01:
        return 'C'
    elif bin_in == 0b10:
        return 'G'
    elif bin_in == 0b11:
        return 'T'
    else: return
#print(transformation(0,0))

# input DNA artificial sequence, output binary data

def rev_transformation(alfa):
    if alfa =='A':
        return "00"
    elif alfa == 'C':
        return "01"
    elif alfa == 'G':
        return "10"
    elif alfa == 'T':
        return "11"
    else: return
#print(rev_transformation('A'))



# input binary number (string format) with even number of bits
def DNA_encoding(e):
    res = ""
    m = len(e)
    for i in range(0,int(m/2)):
        x = transformation(e[2*i],e[2*i+1])
        res += x
    return res
#print(DNA_encoding("100111"))


#input DNA artificial sequence
def DNA_decoding(alfa):
    res = ""
    for i in alfa:
        res += rev_transformation(i)
    return res
#print(DNA_decoding("GCT"))


encoded = DNA_encoding(word_to_binary("BOB ANNA"))
print(encoded)
decoded = binary_to_word(DNA_decoding(encoded))
print(decoded)




