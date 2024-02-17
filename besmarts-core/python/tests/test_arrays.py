
"""
besmarts.tests.test_arrays

"""
import unittest

from besmarts.core.arrays import bitvec
from besmarts.core.arrays import intvec

class test_bitvec(unittest.TestCase):

    def setUp(self):
        self.v = bitvec(4, maxbits=4)

    def test_bitvec_on_off(self):
        v = self.v
        # the bitvec interface must have the on and off members
        self.assertEqual(v.on(), [2])
        self.assertEqual(v.off(), [0,1,3])

    def test_bitvec_int_set(self):
        # the bitvec interface must be indexable and take boolean values
        v = self.v
        v[1] = 1
        self.assertEqual(v.v, 6)
        v[1] = 0
        self.assertEqual(v.v, 4)

    def test_bitvec_bool_set(self):
        v = self.v
        v[1] = True
        self.assertEqual(v.v, 6)
        v[1] = False
        self.assertEqual(v.v, 4)

    # the bitvec interface must support NOT, AND, OR, XOR, EQ operations
    def test_bitvec_bitwise_and(self):
        v = self.v
        self.assertEqual((v & v),  v)
    def test_bitvec_bitwise_or(self):
        v = self.v
        self.assertEqual((v | v),  v)

    def test_bitvec_bitwise_xor(self):
        v = self.v
        self.assertEqual((v ^ v),  bitvec(0))

    def test_bitvec_bitwise_not(self):
        # maxbits is 4, so it is 2**2 - 2**3 - 1 = -5
        v = self.v
        self.assertEqual(~v,  bitvec(-5))

    def test_bitvec_bitwise_equal(self):
        # maxbits is 4, so it is 2**2 - 2**3 - 1 = -5
        v = self.v
        self.assertEqual(v, v)

    def test_bitvec_bitwise_in(self):
        # maxbits is 4, so it is 2**2 - 2**3 - 1 = -5
        v = self.v
        self.assertTrue(v in v)
        self.assertTrue(bitvec(0) in v)
        self.assertFalse(bitvec(0) not in v)

    def test_bitvec_bitwise_not_in(self):

        self.assertFalse(bitvec(3) in bitvec(2))
        self.assertTrue(bitvec(2) in bitvec(3))
        self.assertTrue(bitvec(3) not in bitvec(2))

class test_intvec(unittest.TestCase):

    def setUp(self):
        self.v = intvec()

    def test_intvec_indexing(self):
        # must support indexing
        v = self.v
        v.v.append(1)
        self.assertEqual(v.v[0],  1)
        v.v[0] = 2
        self.assertEqual(v.v[0],  2)

if __name__ == "__main__":
    unittest.main(verbosity=2)
