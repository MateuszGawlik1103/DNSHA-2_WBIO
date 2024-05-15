from enum import Enum


class Nucleotide(Enum):
    A = 0b00
    C = 0b01
    G = 0b10
    T = 0b11

    @property
    def _nucleotide_operation_table(self):
        return {
            ('A', 'A'): 0b00,
            ('A', 'C'): 0b10,
            ('A', 'G'): 0b01,
            ('G', 'C'): 0b11,
        }

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

    def __add__(self, other):
        return Nucleotide((self.value + other.value) % 4), Nucleotide(int((self.value + other.value) / 4))

    def __str__(self):
        return self.name

    def nucleotide_operation(self, other):
        key1 = (self.name, other.name)
        key2 = (other.name, self.name)

        if key1 in self._nucleotide_operation_table:
            return Nucleotide(self._nucleotide_operation_table[key1])
        elif key2 in self._nucleotide_operation_table:
            return Nucleotide(self._nucleotide_operation_table[key2])
        else:
            return None

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

