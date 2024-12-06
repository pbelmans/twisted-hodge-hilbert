r"""
Compute twisted Hodge numbers for Hilbert schemes of points

We implement the _twisted Hodge number conjecture_ as stated in Conjecture E
of [2309.06244] and subsequently proven in [Fu].

It reads

```math
  \sum_{n\geq 0}
  \sum_{p=0}^{2n}
  \sum_{q=0}^{2n}
  \mathrm{h}^{p,q}((\mathop{\rm Hilb}\nolimits^nS,{L}_n)x^py^qt^n
  =
  \prod_{k\ge 1}
  \prod_{p=0}^2
  \prod_{q=0}^2
  \left(
    1-(-1)^{p+q}x^{p+k-1}y^{q+k-1}t^k
  \right)^{-(-1)^{p+q}\mathrm{h}^{p,q}(S,{L}^{\otimes k})}.
```

* [2309.06244] Pieter Belmans, Lie Fu, Andreas Krug:
  Hochschild cohomology of Hilbert schemes of points on surfaces
  [arXiv:2309.06244](https://arxiv.org/abs/2309.06244)
* [Fu] Lie Fu:
  Twisted Hodge numbers and deformation theory of Hilbert schemes of points on surfaces
  via Hodge modules

AUTHORS:

- Pieter Belmans (2024-12-12): initial version
"""

import twisted_ci
from sage.categories.cartesian_product import cartesian_product
from sage.matrix.constructor import matrix
from sage.misc.table import table
from sage.rings.integer_ring import ZZ
from sage.rings.power_series_ring import PowerSeriesRing


def twisted_hodge_diamond(S, n):
    r"""
    Computes the twisted Hodge diamond of a Hilbert scheme of points

    This implements the twisted Hodge number formula describes above, to compute
    the entries of the twisted Hodge diamond for a given number of points,
    i.e., the matrix describing

    ```math
      \mathrm{h}^{p,q}(\mathop{\rm Hilb}\nolimits^nS,L_n)
    ```

    INPUT:

    - ``S`` -- sequence of twisted Hodge diamonds of the surface

    - ``n`` -- the number of points on the Hilbert scheme

    EXAMPLES::

        sage: from twisted_hilbert import *
        sage: P2 = [None, None, None]
        sage: P2[0] = matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        sage: P2[1] = matrix([[10, 0, 0], [8, 0, 0], [1, 0, 0]])
        sage: P2[2] = matrix([[28, 0, 0], [35, 0, 0], [10, 0, 0]])
        sage: twisted_hodge_diamond(P2, 2)
        [55  0  0  0  0]
        [80 28  0  0  0]
        [38 35  0  0  0]
        [ 8 10  0  0  0]
        [ 1  0  0  0  0]

    """
    assert len(S) >= n + 1, "need enough entries in list"
    A = PowerSeriesRing(ZZ, ("x", "y"), default_prec=(2 * n + 1))
    B = PowerSeriesRing(A, "t", default_prec=n + 1)

    x, y = A.gens()
    t = B.gen()

    H = B(1)

    for m in range(1, n + 1):
        for p, q in cartesian_product([range(3), range(3)]):
            H *= (
                1 - (-1) ** (p + q) * x ** (p + m - 1) * y ** (q + m - 1) * t**m
            ) ** (-((-1) ** (p + q)) * S[m][p, q])

    M = matrix.zero(ZZ, 2 * n + 1)

    for p, q in cartesian_product([range(2 * n + 1), range(2 * n + 1)]):
        try:
            M[p, q] = H.dict()[n].dict()[A(x**p * y**q).exponents()[0]]
        except KeyError:
            pass

    return M


class TwistedHodgeDiamond:
    r"""
    Container structure for twisted Hodge diamonds.

    EXAMPLES:

    """

    __M = []

    def __init__(self, M):
        """
        INPUT:

        - ``M`` -- matrix encoding twisted Hodge diamond
        """
        assert M.is_square()

        self.__M = matrix(M)

    def pprint(self):
        T = []
        d = self.__M.nrows()

        for i in reversed(range(2 * d + 1)):
            row = [""] * abs(d - i)

            for j in range(max(0, i - d), min(i, d) + 1):
                row.extend([self.__M[i - j, j], ""])

            T.append(row)

        # padding all rows to full length
        for i in range(len(T)):
            T[i].extend([""] * (2 * d - len(T[i]) + 1))

        return table(T, align="center")

    def __str__(self):
        return str(self.pprint())

    def __getitem__(self, key):
        r"""Return the `(p,q)`th entry.

        INPUT:

        - ``key``: tuple of indices for the twisted Hodge diamond
        """
        p, q = key

        # TODO return TwistedHodgeDiamond?
        return self.__M[p, q]


class TwistedSurfaceDiamonds:
    r"""
    Encodes twisted Hodge diamonds of surface and powers of a line bundle
    """

    __L = None

    def __init__(self):
        pass

    @classmethod
    def from_list(cls, L):
        r"""
        Construct a TwistedSurfaceDiamonds object from a list of matrices

        This is the basic approach, and limits the calculation of twisted Hodge
        numbers to however many entries are provided.

        INPUT:

        - ``L`` -- list of 3x3 matrices with twisted Hodge numbers

        EXAMPLES:

        Twisted Hodge numbers of the projective plane::

            sage: from twisted_hilbert import *
            sage: P2 = [None, None, None]
            sage: P2[0] = matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
            sage: P2[1] = matrix([[10, 0, 0], [8, 0, 0], [1, 0, 0]])
            sage: P2[2] = matrix([[28, 0, 0], [35, 0, 0], [10, 0, 0]])
            sage: S = TwistedSurfaceDiamonds.from_list(P2)

        """
        pair = TwistedSurfaceDiamonds()
        pair.__L = list(map(matrix, L))

        return pair

    def __getitem__(self, k):
        r"""Get the ``k``-twisted Hodge diamond for the line bundle ``L``"""
        # TODO return TwistedHodgeDiamond?
        return self.__L[k]


class CompleteIntersectionSurface(TwistedSurfaceDiamonds):
    r"""
    TwistedSurfaceDiamonds for a complete intersection and the line bundle O(1)
    """

    __d = None
    __i = 1

    def __init__(self, d, i=1):
        try:
            len(d)
        except TypeError:
            d = [d]

        self.__d = d
        self.__i = i

    def __getitem__(self, k):
        return twisted_ci.TwistedHodgeDiamond((len(self.__d) + 2, self.__d), k * self.__i)


class TwistedHilbertSchemeDiamond:
    __diamond = None

    def __init__(self, S: TwistedSurfaceDiamonds, n):
        self.__diamond = twisted_hodge_diamond([S[k] for k in range(n + 1)], n)

    def pprint(self):
        return self.__diamond
