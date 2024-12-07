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
from sage.categories.rings import Rings
from sage.matrix.constructor import matrix
from sage.misc.fast_methods import Singleton
from sage.misc.table import table
from sage.rings.integer_ring import ZZ
from sage.rings.polynomial.multi_polynomial import MPolynomial
from sage.rings.power_series_ring import PowerSeriesRing
from sage.structure.element import Element
from sage.structure.parent import Parent


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

    for k in range(1, n + 1):
        for p, q in cartesian_product([range(3), range(3)]):
            H *= (
                1 - (-1) ** (p + q) * x ** (p + k - 1) * y ** (q + k - 1) * t**k
            ) ** (-((-1) ** (p + q)) * S[k][p, q])

    M = matrix.zero(ZZ, 2 * n + 1)

    for p, q in cartesian_product([range(2 * n + 1), range(2 * n + 1)]):
        try:
            M[p, q] = H.dict()[n].dict()[A(x**p * y**q).exponents()[0]]
        except KeyError:
            pass

    return M


class TwistedHodgeDiamond(Element):
    r"""
    Container structure for twisted Hodge diamonds.

    EXAMPLES:

    The twisted Hodge diamond for the projective plane and anticanonical bundle::

        sage: from twisted_hilbert import *
        sage: H = TwistedHodgeDiamond.from_matrix([[10, 0, 0], [8, 0, 0], [1, 0, 0]])
        sage: H
        twisted Hodge diamond
        sage: H.pprint()
                      0
                  0        0
              1       0        0
                  8        0
                      10

    Notice how (twisted) Hodge diamond are printed in a funny way, with `h^{0,0}` at
    the bottom.

    """

    def __init__(self, parent, M):
        """
        INPUT:

        - ``M`` -- matrix encoding twisted Hodge diamond
        """
        self.__M = M

        Element.__init__(self, parent)

    @classmethod
    def from_matrix(cls, M):
        r"""
        Construct a twisted Hodge diamond from a matrix

        INPUT:

        - ``M`` -- square matrix encoding twisted Hodge diamond

        EXAMPLES:

        The twisted Hodge diamond for the projective plane and anticanonical bundle::

            sage: from twisted_hilbert import *
            sage: H = TwistedHodgeDiamond.from_matrix([[10, 0, 0], [8, 0, 0], [1, 0, 0]])
            sage: H.pprint()
                          0
                      0        0
                  1       0        0
                      8        0
                          10
        """

        M = matrix(M)
        assert M.is_square()

        return TwistedHodgeDiamondRing()(M)

    def pprint(self):
        r"""Pretty print the twisted Hodge diamond

        EXAMPLES:

        The twisted Hodge diamond for the projective plane and anticanonical bundle::

            sage: from twisted_hilbert import *
            sage: H = TwistedHodgeDiamond.from_matrix([[10, 0, 0], [8, 0, 0], [1, 0, 0]])
            sage: H.pprint()
                          0
                      0        0
                  1       0        0
                      8        0
                          10
        """
        T = []
        d = self.__M.nrows() - 1

        for i in reversed(range(2 * d + 1)):
            row = [""] * abs(d - i)

            for j in range(max(0, i - d), min(i, d) + 1):
                row.extend([self.__M[i - j, j], ""])

            T.append(row)

        # padding all rows to full length
        for i in range(len(T)):
            T[i].extend([""] * (2 * d - len(T[i]) + 1))

        return table(T, align="center")

    def _repr_(self):
        r"""Output diagnostic information

        This is a one-line string giving some basic information about the twisted Hodge
        diamond. You'll see this when you just evaluate something which returns
        a twisted Hodge diamond. To see something more useful, you'll likely want to use

        * :meth:`TwistedHodgeDiamond.__str__` via `print`
        * :meth:`TwistedHodgeDiamond.pprint`
        * :meth:`TwistedHodgeDiamond.polynomial`

        It is also possible to override this output by using the built-in
        functionality for parents and renaming.

        EXAMPLES:

        The default behavior::

            sage: from twisted_hilbert import *
            sage: TwistedHodgeDiamond.from_matrix([[10, 0, 0], [8, 0, 0], [1, 0, 0]])
            twisted Hodge diamond

        Special constructors give additional information::

            sage: CompleteIntersectionSurface(3, 2)[2]
            twisted Hodge diamond for complete intersection of degree [3] and L=O(4)

        """
        return "twisted Hodge diamond"

    def __str__(self):
        r"""Pretty print a twisted Hodge diamond

        This gets called when you specifically print the object.

        EXAMPLES:

        The twisted Hodge diamond for the projective plane and anticanonical bundle::

            sage: from twisted_hilbert import *
            sage: H = TwistedHodgeDiamond.from_matrix([[10, 0, 0], [8, 0, 0], [1, 0, 0]])
            sage: print(H)
                          0
                      0        0
                  1       0        0
                      8        0
                          10

        """
        return str(self.pprint())

    def __getitem__(self, key):
        r"""Return the `(p,q)`th entry.

        This is

        ```math
        \dim\mathrm{H}^q(X,\Omega_X^p\otimes L)
        ```

        corresponding to the entry indexed by `p` and `q` in the matrix.

        INPUT:

        - ``key``: tuple of indices for the twisted Hodge diamond

        EXAMPLES:

        The twisted Hodge diamond for the projective plane and anticanonical bundle::

            sage: from twisted_hilbert import *
            sage: H = TwistedHodgeDiamond.from_matrix([[10, 0, 0], [8, 0, 0], [1, 0, 0]])
            sage: H[0, 0]
            10
            sage: H[1, 0]
            8

        """
        return self.__M[key[0], key[1]]


class TwistedHodgeDiamondRing(Singleton, Parent):
    def __init__(self):
        Parent.__init__(self, category=Rings().Commutative())

    def _element_constructor_(self, *args, **keywords):
        m = args[0]
        if m in ZZ:
            m = matrix(ZZ, 1, 1, [m])
        elif isinstance(m, MPolynomial):
            raise NotImplementedError()
        elif isinstance(m, (list, tuple)):
            m = matrix(m)
        elt = self.element_class(self, m)
        return elt

    def from_matrix(self, M):
        return self.element_class(self, matrix(M))

    def _repr_(self) -> str:
        return "Ring of twisted Hodge diamonds"

    Element = TwistedHodgeDiamond


class TwistedSurfaceDiamonds:
    r"""Encodes twisted Hodge diamonds of surface and powers of a line bundle"""

    def __init__(self):
        pass

    @classmethod
    def from_list(cls, diamonds):
        r"""
        Construct a TwistedSurfaceDiamonds object from a list of matrices

        This is the basic approach, and limits the calculation of twisted Hodge
        numbers to however many entries are provided.

        INPUT:

        - ``diamonds`` -- list of 3x3 matrices with twisted Hodge numbers

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
        pair.__diamonds = list(map(TwistedHodgeDiamond.from_matrix, diamonds))

        return pair

    def __getitem__(self, k):
        r"""Get the twisted Hodge diamond for the line bundle ``L^k``"""
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
        # TODO make this cleaner...
        H = twisted_ci.TwistedHodgeDiamond((len(self.__d) + 2, self.__d), k * self.__i)
        H = TwistedHodgeDiamond.from_matrix(H._TwistedHodgeDiamond__M)
        H.rename(
            "twisted Hodge diamond for complete intersection of "
            + f"degree {self.__d} and L=O({self.__i * k})"
        )

        return H


class BiellipticSurface(TwistedSurfaceDiamonds):
    r"""
    TwistedSurfaceDiamonds for a bielliptic surface and the line bundle O(1)

    Here O(1) refers to the fraction of the canonical bundle.
    """

    def __init__(self, order):
        assert order in [2, 3, 4, 6]

        self.__order = order

    def __getitem__(self, k):
        if self.__order == 2:
            diamonds = [
                matrix([[1, 1, 0], [1, 2, 1], [0, 1, 1]]),
                matrix([[0, 1, 1], [1, 2, 1], [1, 1, 0]]),
            ]
            return TwistedHodgeDiamond.from_matrix(diamonds[k % 2])
        if self.__order == 3:
            diamonds = [
                matrix([[1, 1, 0], [1, 2, 1], [0, 1, 1]]),
                matrix([[0, 0, 0], [1, 1, 0], [1, 1, 0]]),
                matrix([[0, 1, 1], [0, 1, 1], [0, 0, 0]]),
            ]
            return TwistedHodgeDiamond.from_matrix(diamonds[k % 3])
        if self.__order in [4, 6]:
            raise NotImplementedError()


class EnriquesSurface(TwistedSurfaceDiamonds):
    r"""
    TwistedSurfaceDiamonds for an Enriques surface and the (anti)canonical line bundle

    EXAMPLES:

    The following is Appendix B of [AJM.2017.v21.n6.a4]

        sage: from twisted_hilbert import *
        sage: TwistedHilbertSchemeDiamond(EnriquesSurface(), 2)[3, 1]
        10

    * [AJM.2017.v21.n6.a4] Taro Hayashi
      Universal covering Calabi–Yau manifolds of the Hilbert schemes of points of
      Enriques surfaces
      https://dx.doi.org/10.4310/AJM.2017.v21.n6.a4
    """

    def __getitem__(self, k):
        if k % 2 == 0:
            return TwistedHodgeDiamond.from_matrix([[1, 0, 0], [0, 10, 0], [0, 0, 1]])
        if k % 2 == 1:
            return TwistedHodgeDiamond.from_matrix([[0, 0, 1], [0, 10, 0], [1, 0, 0]])


def TwistedHilbertSchemeDiamond(S: TwistedSurfaceDiamonds, n):
    r"""
    Construct twisted Hodge diamond Hilbert schemes of $n$ points on $S$

    EXAMPLES:

    Anticanonically twisted Hodge diamond of Hilb^2 P^2::

        sage: from twisted_hilbert import *
        sage: S = CompleteIntersectionSurface([], 3)
        sage: H = TwistedHilbertSchemeDiamond(S, 2)
        sage: H
        twisted Hodge diamond for Hilb^2 S
        sage: H.pprint()
                            0
                       0         0
                  0         0        0
              0        0         0       0
          1       10        0        0       0
              8        35        0       0
                  38        28       0
                       80        0
                            55

    Anticanonically twisted Hodge diamond of Hilb^2 S where S is bielliptic order 2::

        sage: TwistedHilbertSchemeDiamond(BiellipticSurface(2), 2).pprint()
                          0
                      0       0
                  0       2       0
              1       4       4       1
          1       3       8       3       1
              1       4       4       1
                  0       2       0
                      0       0
                          0

    Anticanonically twisted Hodge diamond of Hilb^2 S where S is bielliptic order 3::

        sage: TwistedHilbertSchemeDiamond(BiellipticSurface(3), 2).pprint()
                          0
                      0       0
                  0       0       0
              1       1       1       0
          1       2       2       1       0
              1       1       1       0
                  0       0       0
                      0       0
                          0

    """
    H = TwistedHodgeDiamond.from_matrix(
        twisted_hodge_diamond([S[k] for k in range(n + 1)], n)
    )
    H.rename(f"twisted Hodge diamond for Hilb^{n} S")

    return H
