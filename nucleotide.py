from enum import Enum

class Nucleotide(Enum):
    A = 0b00
    C = 0b01
    G = 0b10
    T = 0b11

    def __invert__(self):
        return Nucleotide(~self.value & 0b11)

    def __and__(self, other):
        return Nucleotide(self.value & other.value)

    def __or__(self, other):
        return Nucleotide(self.value | other.value)

    def __xor__(self, other):
        return Nucleotide(self.value ^ other.value)

    def __repr__(self):
        return self.name

    @classmethod
    def from_letter(cls, letter):
        letter = letter.upper()
        if letter == 'A':
            return cls.A
        elif letter == 'C':
            return cls.C
        elif letter == 'G':
            return cls.G
        elif letter == 'T':
            return cls.T
        else:
            raise ValueError("Invalid nucleotide letter: {}".format(letter))
