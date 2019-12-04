#!/usr/bin/env python3

import unittest
import numpy as np
import huckel


class TestHuckel(unittest.TestCase):

    def test_lin_polyene_n(self):
        for sites in (1, 2, 3, 4, 5, 6, 10, 20, 30, 50, 100):
            m = huckel.lin_polyene_n(sites)
            evals = np.linalg.eigvalsh(m)

            self.assertEqual(len(evals), sites)

            for n, val in enumerate(evals):
                expected = -2*np.cos(np.pi * (n+1) / (sites + 1))  # A4 part 1, page 10 eq. (25)
                self.assertEqual(round(val, 10), round(expected, 10))

        for n in (-10, -3, -1, 0):
            self.assertRaises(ValueError, huckel.lin_polyene_n, n)

    def test_cyc_polyene_n(self):
        for sites in (3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 50, 100):
            m = huckel.cyc_polyene_n(sites)
            evals = np.linalg.eigvalsh(m)

            self.assertEqual(len(evals), sites)

            expected = []
            n = - sites // 2 + 1

            while n <= sites // 2:
                expected.append(-2*np.cos(2 * np.pi * n / sites))  # A4 part 1, page 11 eq. (30)
                n += 1

            expected.sort()

            for n, val in enumerate(evals):
                self.assertEqual(round(val, 10), round(expected[n], 10))

        for n in (-10, -3, -1, 0, 1, 2):
            self.assertRaises(ValueError, huckel.cyc_polyene_n, n)

    def test_get_eigenvalues_with_degeneracies(self):
        self.assertEqual(huckel.get_eigenvalues_with_degeneracies((2, 2, 2, 2, 2)), [(2, 5)])
        self.assertEqual(huckel.get_eigenvalues_with_degeneracies((2, 2, 2, 3, 3)), [(2, 3), (3, 2)])
        self.assertEqual(huckel.get_eigenvalues_with_degeneracies((7, 8, 9)), [(7, 1), (8, 1), (9, 1)])

        # Check the degeneracies for some cyclic poly-enes.
        for sites in (3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 50, 100):
            m = huckel.cyc_polyene_n(sites)
            evals = huckel.get_eigenvalues_with_degeneracies(np.linalg.eigvalsh(m))

            n = 0

            for _, degen in evals:
                if n == 0 or n == sites - 1:
                    self.assertEqual(degen, 1)
                else:
                    self.assertEqual(degen, 2)

                n += degen

    def test_platonic(self):
        for n in (4, 6, 8, 12, 20):
            m = huckel.platonic(n)
            evals = np.linalg.eigvalsh(m)

            self.assertEqual(len(evals), n)

        for n in (3, 5, 9, 10):
            self.assertRaises(ValueError, huckel.platonic, n)


if __name__ == '__main__':
    unittest.main()
