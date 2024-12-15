import unittest
from fractions import Fraction

from . import QuadraticRadical


class MyTestCase(unittest.TestCase):
    def test_init(self):
        self.assertTrue(QuadraticRadical(1))
        self.assertTrue(QuadraticRadical((3, 2)))
        self.assertFalse(QuadraticRadical())
        self.assertFalse(QuadraticRadical((4, 3), (-4, 3)))
        self.assertEqual(QuadraticRadical((1, 2), (2, 3)),
                         QuadraticRadical((2, 3), (1, 2)))

    def test_signs(self):
        self.assertEqual(QuadraticRadical(), -QuadraticRadical(0))
        self.assertNotEqual(QuadraticRadical(3), -QuadraticRadical(3))
        a: QuadraticRadical = QuadraticRadical((6, 7))
        self.assertEqual(a, +a)
        self.assertEqual(-a, ---a)

    def test_casts(self):
        self.assertEqual(int(QuadraticRadical((Fraction(17, 6), 2))), 4)
        self.assertEqual(float(QuadraticRadical((Fraction(25, 16), 1))),1.5625)

    def test_add(self):
        self.assertEqual(QuadraticRadical(1)+1, 2)


if __name__ == '__main__':
    unittest.main()
