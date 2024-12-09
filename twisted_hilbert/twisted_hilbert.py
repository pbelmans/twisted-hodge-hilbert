r"""
Compute twisted Hodge numbers for Hilbert schemes of points

We implement the twisted Hodge number formula for Hilbert schemes of points on surfaces,
as

* stated in Conjecture E of [2309.06244],
* proven in [Fu].

It reads

.. MATH::

  \sum_{n\geq 0}
  \sum_{p=0}^{2n}
  \sum_{q=0}^{2n}
  \mathrm{h}^{p,q}(\mathop{\rm Hilb}\nolimits^nS,{L}_n)x^py^qt^n
  =
  \prod_{k\ge 1}
  \prod_{p=0}^2
  \prod_{q=0}^2
  \left(
    1-(-1)^{p+q}x^{p+k-1}y^{q+k-1}t^k
  \right)^{-(-1)^{p+q}\mathrm{h}^{p,q}(S,{L}^{\otimes k})}.

Here, $L$ is a line bundle on a smooth projective surface $S$ (or a compact complex
surface), $\\operatorname{Hilb}^nS$ is the Hilbert scheme of $n$ points on $S$,
and $L_n$ is the induced line bundle on the Hilbert scheme. One is referred to [Fu]
for more details.

This formula computes the twisted Hodge numbers of $L_n$ on the Hilbert scheme, i.e.,

.. MATH::

  \mathrm{h}^{p,q}(\mathop{\rm Hilb}\nolimits^nS,L_n)

An interesting example where things can be computed also using an explicit description
of the Hilbert scheme is $\\operatorname{Hilb}^2\\mathbb{P}^2$, where the Hochschild
cohomology (or rather, Hochschild--Kostant--Rosenberg decomposition) is computed,
which corresponds to the twisted Hodge diamond using the anticanonical line bundle,
cf. [Section 4.2, 2309.06244]::

    sage: from twisted_hilbert import *
    sage: S = CompleteIntersectionSurface([], 3) # anticanonical twist
    sage: TwistedHilbertSchemeDiamond(S, 2).as_parallelogram()
      1
      0   8
      0   10   38
      0   0    35   80
      0   0    0    28   55
          0    0    0    0
               0    0    0
                    0    0
                         0

* [2309.06244] Pieter Belmans, Lie Fu, Andreas Krug,
  Hochschild cohomology of Hilbert schemes of points on surfaces
  `arXiv:2309.06244 <https://arxiv.org/abs/2309.06244>`_
* [Fu] Lie Fu:
  Twisted Hodge numbers and deformation theory of Hilbert schemes of points on surfaces
  via Hodge modules

AUTHORS:

- Pieter Belmans (2024-12-12): initial version
"""

import twisted_ci
from sage.arith.misc import binomial
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

    This implements the twisted Hodge number formula described above, to compute
    the entries of the twisted Hodge diamond for a given number of points,
    i.e., the matrix describing

    .. MATH::

      \mathrm{h}^{p,q}(\mathop{\rm Hilb}\nolimits^nS,L_n)

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

    Notice how (twisted) Hodge diamond are printed in a funny way, with
    $\\mathrm{h}^{0,0}$ at the bottom.

    """

    def __init__(self, parent, M):
        """
        Constructor for a TwistedHodgeDiamond (not to be called directly)

        INPUT:

        - ``M`` -- matrix encoding twisted Hodge diamond

        This function should not be called directly, use the class method
        :meth:`TwistedHodgeDiamond.from_matrix` instead.
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

    def dimension(self):
        r"""Dimension of the variety underlying the twisted Hodge diamond

        EXAMPLES::

            sage: from twisted_hilbert import *
            sage: EnriquesSurface()[0].dimension()
            2
            sage: TwistedHilbertSchemeDiamond(EnriquesSurface(), 3).dimension()
            6

        """
        return self.__M.nrows() - 1

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
        d = self.dimension()

        for i in reversed(range(2 * d + 1)):
            row = [""] * abs(d - i)

            for j in range(max(0, i - d), min(i, d) + 1):
                row.extend([self.__M[i - j, j], ""])

            T.append(row)

        # padding all rows to full length
        for i in range(len(T)):
            T[i].extend([""] * (2 * d - len(T[i]) + 1))

        return table(T, align="center")

    def as_parallelogram(self):
        r"""Return the twisted Hodge diamond as polyvector parallelogram

        It is up to the user to make sure that the twisted Hodge diamond is computed
        using the anticanonical bundle.

        EXAMPLES:

        Our favourite example is still $\\operatorname{Hilb}^2\\mathbb{P}^2$::

            sage: from twisted_hilbert import *
            sage: H = TwistedHilbertSchemeDiamond(CompleteIntersectionSurface([], 3), 2)
            sage: H.as_parallelogram()
              1
              0   8
              0   10   38
              0   0    35   80
              0   0    0    28   55
                  0    0    0    0
                       0    0    0
                            0    0
                                 0
        """
        T = []
        d = self.dimension()

        for n in range(2 * d + 1):
            T.append([])

            for p in range(d + 1):
                q = n - p
                if q in range(d + 1) and p in range(d + 1):
                    T[-1].append(self[d - p, q])
                else:
                    T[-1].append("")

        return table(T, align="center")

    def _repr_(self):
        r"""Output diagnostic information

        This is a one-line string giving some basic information about the twisted Hodge
        diamond. You'll see this when you just evaluate something which returns
        a twisted Hodge diamond. To see something more useful, you'll likely want to use

        * :meth:`TwistedHodgeDiamond.__str__` via `print`
        * :meth:`TwistedHodgeDiamond.pprint`

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

    def __eq__(self, other):
        r"""Compare two twisted Hodge diamonds

        INPUT:

        - ``other`` -- the other twisted Hodge diamond

        EXAMPLES:

        Twisted Hodge diamonds for bielliptic surfaces are (not) the same::

            sage: from twisted_hilbert import *
            sage: BiellipticSurface(2)[0] == BiellipticSurface(3)[0]
            True
            sage: BiellipticSurface(2)[1] == BiellipticSurface(3)[1]
            False

        """

        return self.__M == other.__M

    def __getitem__(self, key):
        r"""Return $\\mathrm{h}^{p,q}(X,L)$

        This is

        .. MATH::

            \dim\mathrm{H}^q(X,\Omega_X^p\otimes L)

        corresponding to the entry indexed by $p$ and $q$ in the matrix, where we
        take `key=(p,q)`.

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
    r"""Encodes twisted Hodge diamonds of surface and powers of a line bundle

    This makes it possible to implement both a class that knows all about a surface
    and one that only contains a finite amount of data.
    """

    @classmethod
    def from_list(cls, diamonds):
        r"""
        Construct a :class:`TwistedSurfaceDiamonds` object from a list of matrices

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
    :class:`TwistedSurfaceDiamonds` for a complete intersection

    It is possible to vary the line bundle being used to any power of $\\mathcal{O}(1)$
    """

    __d = None
    __i = 1

    def __init__(self, d, i=1):
        r"""
        Construct a complete intersection

        The twisted Hodge numbers for complete intersection surfaces are computed
        using [twisted-hodge-ci].

        * [twisted-hodge-ci] Twisted Hodge numbers for complete intersections
          `twisted-hodge-ci <https://github.com/pbelmans/twisted-hodge-ci>`_

        INPUT:

        - ``d`` -- degree, or list of degrees

        - ``i`` (default: 1) -- power of $\\mathcal{O}(1)$ to be used

        EXAMPLES:

        Anticanonically twisted projective plane::

            sage: from twisted_hilbert import *
            sage: CompleteIntersectionSurface([], 3)[1].pprint()
                      0
                  0        0
              1       0        0
                  8        0
                      10

        Quadric surface with default twist::

            sage: CompleteIntersectionSurface(2)[1].pprint()
                      0
                  0       0
              0       0       0
                  0       0
                      4

        del Pezzo surface of degree 4 with default (and anticanonical) twist::

            sage: CompleteIntersectionSurface([2, 2])[1].pprint()
                      0
                  0       0
              1       2       0
                  0       0
                      5

        """

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


class ProductSurface(TwistedSurfaceDiamonds):
    r"""
    :class:`TwistedSurfaceDiamonds` for the product of two curves

    The line bundle is the anticanonical bundle.
    """

    def __init__(self, g, h):
        assert g >= 0 and h >= 0, "genera need to be positive"

        self.__g = g
        self.__h = h

    @classmethod
    def __H(cls, g, i):
        r"""Cohomology of $i$th power of anticanonical line bundle on genus $g$ curve

        Helper function, not to be called directly.

        EXAMPLES:

            The projective line:

                sage: from twisted_hilbert import *
                sage: ProductSurface._ProductSurface__H(0, 0)
                [1, 0]
                sage: ProductSurface._ProductSurface__H(0, 3)
                [7, 0]
                sage: ProductSurface._ProductSurface__H(0, -3)
                [0, 5]

            On an elliptic curve:

                sage: ProductSurface._ProductSurface__H(1, 0)
                [1, 1]
                sage: ProductSurface._ProductSurface__H(1, 5)
                [1, 1]

            On a curve of genus 3:

                sage: from twisted_hilbert import *
                sage: ProductSurface._ProductSurface__H(3, 0)
                [1, 3]
                sage: ProductSurface._ProductSurface__H(3, -1)
                [3, 1]
                sage: ProductSurface._ProductSurface__H(3, 3)
                [0, 14]
                sage: ProductSurface._ProductSurface__H(3, -3)
                [10, 0]

        """
        if g == 0:
            if i <= -1:
                return [0, -2 * i - 1]
            if i >= 0:
                return [2 * i + 1, 0]
        if g == 1:
            return [1, 1]
        if g >= 2:
            if i <= -2:
                return [(-2 * i - 1) * (g - 1), 0]
            if i == -1:
                return [g, 1]
            if i == 0:
                return [1, g]
            if i >= 1:
                return [0, (2 * i + 1) * (g - 1)]

    def __getitem__(self, k):
        r"""Get the twisted Hodge diamond for the kth power of the anticanonical bundle

        EXAMPLES:

        The quadric is the product of two curves of genus 0::

            sage: from twisted_hilbert import *
            sage: S = ProductSurface(0, 0)
            sage: T = CompleteIntersectionSurface(2, 2)
            sage: T[0] == S[0]
            True
            sage: T[3] == S[3]
            True
            sage: T[-5] == S[-5]
            True

        """
        # introduce shorthands
        g = self.__g
        h = self.__h
        H = ProductSurface.__H

        return TwistedHodgeDiamond.from_matrix(
            [
                [
                    H(g, k)[0] * H(h, k)[0],
                    H(g, k)[0] * H(h, k)[1] + H(g, k)[1] * H(h, k)[0],
                    H(g, k)[1] * H(h, k)[1],
                ],
                [
                    H(g, k)[0] * H(h, k - 1)[0] + H(g, k - 1)[0] * H(h, k)[0],
                    H(g, k)[0] * H(h, k - 1)[1]
                    + H(g, k)[1] * H(h, k - 1)[0]
                    + H(g, k - 1)[0] * H(h, k)[1]
                    + H(g, k - 1)[1] * H(h, k)[0],
                    H(g, k)[1] * H(h, k - 1)[1] + H(g, k - 1)[1] * H(h, k)[1],
                ],
                [
                    H(g, k - 1)[0] * H(h, k - 1)[0],
                    H(g, k - 1)[0] * H(h, k - 1)[1] + H(g, k - 1)[1] * H(h, k - 1)[0],
                    H(g, k - 1)[1] * H(h, k - 1)[1],
                ],
            ]
        )


class BiellipticSurface(TwistedSurfaceDiamonds):
    r"""
    :class:`TwistedSurfaceDiamonds` for a bielliptic surface

    The line bundle is the anticanonical line bundle
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
    :class:`TwistedSurfaceDiamonds` for an Enriques surface

    The line bundle is the (anti)canonical line bundle

    EXAMPLES:

    The following is Appendix B of [AJM.2017.v21.n6.a4]::

        sage: from twisted_hilbert import *
        sage: TwistedHilbertSchemeDiamond(EnriquesSurface(), 2)[3, 1]
        10
        sage: TwistedHilbertSchemeDiamond(EnriquesSurface(), 2).pprint()
                           0
                       0        0
                  0        1        0
              0        0        0        0
          1       10       66       10       1
              0        0        0        0
                  0        1        0
                       0        0
                           0

    * [AJM.2017.v21.n6.a4] Taro Hayashi,
      Universal covering Calabi–Yau manifolds of the Hilbert schemes of points of
      Enriques surfaces
      `AJM.2017.v21.n6.a4 <https://dx.doi.org/10.4310/AJM.2017.v21.n6.a4>`_
    """

    def __getitem__(self, k):
        if k % 2 == 0:
            return TwistedHodgeDiamond.from_matrix([[1, 0, 0], [0, 10, 0], [0, 0, 1]])
        if k % 2 == 1:
            return TwistedHodgeDiamond.from_matrix([[0, 0, 1], [0, 10, 0], [1, 0, 0]])


def HilbertSchemeDeformations(HKR: TwistedHodgeDiamond):
    r"""
    The degree 0, 1 and 2 cohomology of the tangent bundle of a Hilbert scheme

    This implements Theorem 1.5 of [Fu]. It it expected that the user inputs the
    anticanonical Hodge diamond, describing the Hochschild cohomology of the surface.

    For this one doesn't have to be able to compute all the twisted Hodge diamonds
    required to describe the twisted Hodge diamond giving all of the Hochschild
    cohomology of the Hilbert scheme, a single twisted Hodge diamond, for the
    anticanonical twist, suffices.

    EXAMPLES:

    Our favourite example is still $\\operatorname{Hilb}^2\\mathbb{P}^2$::

        sage: from twisted_hilbert import *
        sage: HilbertSchemeDeformations(CompleteIntersectionSurface([], 3)[1])
        [8, 10, 0]

    """
    return [
        HKR[1, 0],
        HKR[1, 1] + HKR[1, 0] * HKR[2, 1] + HKR[0, 0],
        HKR[1, 2]
        + HKR[1, 1] * HKR[2, 1]
        + HKR[1, 0] * HKR[2, 2]
        + HKR[1, 0] * binomial(HKR[2, 1], 2)
        + HKR[2, 1] * HKR[0, 0]
        + HKR[0, 1],
    ]


def TwistedHilbertSchemeDiamond(S: TwistedSurfaceDiamonds, n):
    r"""
    Construct twisted Hodge diamond Hilbert schemes of $n$ points on $S$

    INPUT:

    - ``S`` -- twisted Hodge diamonds for powers of a line bundle

    - ``n`` -- number of points

    EXAMPLES:

    Anticanonically twisted Hodge diamond of $\\operatorname{Hilb}^2\\mathbb{P}^2$::

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

    Anticanonically twisted Hodge diamond of $\\operatorname{Hilb}^2S$ where $S$ is a
    bielliptic surface with canonical bundle of order 2::

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

    Anticanonically twisted Hodge diamond of $\\operatorname{Hilb}^2S$ where $S$ is a
    bielliptic surface with canonical bundle of order 3::

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
