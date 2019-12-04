import numpy as np
import networkx as nx
import sys

from decimal import Decimal
from typing import Iterable, Tuple, List


def lin_polyene_n(n: int) -> np.ndarray:
    """
    Returns the hamiltonian matrix for a linear poly-ene with n sites. The matrix entries are as per the Huckel
    approximation, with alpha = 0 and beta = -1.
    """
    if n < 1:
        raise ValueError("Linear poly-ene must have at least 1 site.")

    m = np.array([[-1 if i == j - 1 or i == j + 1 else 0 for i in range(n)] for j in range(n)])

    return m


def cyc_polyene_n(n: int) -> np.ndarray:
    """
    Returns the hamiltonian matrix for a cyclic poly-ene with n sites. This is just the same as for a linear polyene,
    except now the first and last atom are connected.
    """
    if n < 1:
        raise ValueError("Cyclic poly-ene must have at least 1 site.")

    m = lin_polyene_n(n)
    m[0][n-1] = -1
    m[n-1][0] = -1
    return m


def platonic(n: int) -> np.ndarray:
    """
    Returns the hamiltonian matrix corresponding to the (sp2-hybridised) platonic solid with n vertices.
    Possible values are n = 4, 6, 8, 12, 20.  Use the pre-defined graphs from networkx, as there is no nice way of
    generating them algorithmically (they are just definitions, after all).
    """

    if n == 4:
        g = nx.tetrahedral_graph()
    elif n == 6:
        g = nx.cubical_graph()
    elif n == 8:
        g = nx.octahedral_graph()
    elif n == 12:
        g = nx.dodecahedral_graph()
    elif n == 20:
        g = nx.icosahedral_graph()
    else:
        raise ValueError("Unknown platonic solid")

    mat = nx.adjacency_matrix(g)
    return - mat.todense()


def buckyball() -> np.ndarray:
    """
    Return the hamiltonian matrix (in the Huckel approximation) for buckminsterfullerene, C60. Like for the platonic
    solids, there's no straightforward way to generate this algorithmically, so just return it hard-coded.
    """
    edges = [(0, 2), (0, 48), (0, 59), (1, 3), (1, 9), (1, 58),
             (2, 3), (2, 36), (3, 17), (4, 6), (4, 8), (4, 12),
             (5, 7), (5, 9), (5, 16), (6, 7), (6, 20), (7, 21),
             (8, 9), (8, 56), (10, 11), (10, 12), (10, 20), (11, 27),
             (11, 47), (12, 13), (13, 46), (13, 54), (14, 15), (14, 16),
             (14, 21), (15, 25), (15, 41), (16, 17), (17, 40), (18, 19),
             (18, 20), (18, 26), (19, 21), (19, 24), (22, 23), (22, 31),
             (22, 34), (23, 25), (23, 38), (24, 25), (24, 30), (26, 27),
             (26, 30), (27, 29), (28, 29), (28, 31), (28, 35), (29, 44),
             (30, 31), (32, 34), (32, 39), (32, 50), (33, 35), (33, 45),
             (33, 51), (34, 35), (36, 37), (36, 40), (37, 39), (37, 52),
             (38, 39), (38, 41), (40, 41), (42, 43), (42, 46), (42, 55),
             (43, 45), (43, 53), (44, 45), (44, 47), (46, 47), (48, 49),
             (48, 52), (49, 53), (49, 57), (50, 51), (50, 52), (51, 53),
             (54, 55), (54, 56), (55, 57), (56, 58), (57, 59), (58, 59)
             ]
    g = nx.Graph()
    g.add_edges_from(edges)

    return nx.adjacency_matrix(g).todense()


def get_eigenvalues_with_degeneracies(evals: Iterable[float]) -> Tuple[Decimal, int]:
    """
    Given a set of sorted eigenvalues (possibly including degenerate eigenvalues), return a list of
    (eigenvalue, degeneracy) pairs, with eigenvalues represented as decimals rounded to 3dp.
    """
    cur = None
    count = 0
    result = []

    for e in evals:
        e = Decimal(e)
        e = round(e, 3)

        if e == cur:
            count += 1
        else:
            if cur is not None:
                result.append((cur, count))
            cur = e
            count = 1

    if count > 0:
        result.append((cur, count))

    return result


def print_eigenvalues(evals: List[Tuple[Decimal, int]]) -> None:
    """
    Format and print a sorted (in ascending order of eigenvalue) list of (eigenvalue, degeneracy) pairs as an
    energy-level diagram.
    """
    max_degen = np.amax([e[1] for e in evals])
    line_length = 4 * max_degen - 2
    count = 0

    for val, degen in reversed(evals):
        count += degen
        line = ''
        spacing = (line_length - (4 * degen - 2)) // 2
        line += ' ' * spacing
        line += '――  ' * degen
        line += ' ' * spacing
        if val < 0:
            line += '-'
        else:
            line += ' '
        line += str(abs(val))
        print(line)
        print()

    print("%d orbitals." % count)


def print_usage() -> None:
    print("usage: huckel.py [-l | --linear-polyene] num\n"
          "       huckel.py [-c | --cyclic-polyene] num\n"
          "       huckel.py [-p | --platonic] [4 | 6 | 8 | 12 | 20]\n"
          "       huckel.py [-b | --buckyball]")


def main() -> None:
    if len(sys.argv) < 2:
        print_usage()
        return

    matrix = None

    arguments = {
        '-l': lin_polyene_n, '--linear-polyene': lin_polyene_n,
        '-c': cyc_polyene_n, '--cyclic-polyene': cyc_polyene_n,
        '-p': platonic, '--platonic': platonic
    }

    if sys.argv[1] in arguments.keys():
        try:
            matrix = arguments.get(sys.argv[1])(int(sys.argv[2]))
        except (IndexError, ValueError) as e:
            print_usage()
            return
    elif sys.argv[1] == '-b' or sys.argv[1] == '--buckyball':
        matrix = buckyball()

    if matrix is not None:
        print_eigenvalues(get_eigenvalues_with_degeneracies(np.linalg.eigvalsh(matrix)))
    else:
        print_usage()


if __name__ == '__main__':
    main()
