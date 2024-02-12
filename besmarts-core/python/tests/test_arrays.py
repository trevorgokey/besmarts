
"""
besmarts.tests.test_arrays

"""

from besmarts.core.arrays import bitvec, intvec

# the bitvec interface must be constructed with an integer
v = bitvec(4, maxbits=4)

# the bitvec interface must have the on and off members
assert v.on() == [2]

assert v.off() == [0,1,3]

# the bitvec interface must be indexable and take boolean values
v[1] = 1
assert v.v == 6
v[1] = 0
assert v.v == 4

v[1] = True
assert v.v == 6
v[1] = False
assert v.v == 4

# the bitvec interface must support NOT, AND, OR, XOR, EQ operations
assert (v & v) == v
assert (v | v) == v
assert (v ^ v) == bitvec(0) 
# maxbits is 4, so it is 2**2 - 2**3 - 1 = -5
assert ~v == bitvec(-5) 

# 
v = intvec()

# must support indexing
v.v.append(1)
v.v[0] = 2
