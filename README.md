[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14334379.svg)](https://doi.org/10.5281/zenodo.14334379)
[![license](https://badgen.net/github/license/pbelmans/twisted-hodge-hilbert)](https://github.com/pbelmans/twisted-hodge-hilbert/blob/master/LICENSE)

[![tests](https://github.com/pbelmans/twisted-hodge-hilbert/actions/workflows/tests.yml/badge.svg)](https://github.com/pbelmans/twisted-hodge-hilbert/actions)

# Twisted Hodge numbers for Hilbert schemes of points

We implement the twisted Hodge number formula for Hilbert schemes of points on surfaces,
as

* stated in Conjecture E of [2309.06244],
* proven in Theorem 1.1 of [2412.09975]

It reads

```math
  \sum_{n\geq 0}\sum_{p=0}^{2n}\sum_{q=0}^{2n}\mathrm{h}^{p,q}(\mathop{\rm Hilb}\nolimits^nS,{L}_n)x^py^qt^n
  =
  \prod_{k\ge 1}\prod_{p=0}^2\prod_{q=0}^2\left( 1-(-1)^{p+q}x^{p+k-1}y^{q+k-1}t^k\right)^{-(-1)^{p+q}\mathrm{h}^{p,q}(S,{L}^{\otimes k})}.
```

Here, $L$ is a line bundle on a smooth projective surface $S$ (or a compact complex
surface), $\mathop{\rm Hilb}^nS$ is the Hilbert scheme of $n$ points on $S$,
and $L_n$ is the induced line bundle on the Hilbert scheme. One is referred to [2412.09975]
for more details.

This formula computes the twisted Hodge numbers of $L_n$ on the Hilbert scheme, i.e.,

```math
  \mathrm{h}^{p,q}(\mathop{\rm Hilb}\nolimits^nS,L_n)
```

An interesting example where things can be computed also using an explicit description
of the Hilbert scheme is $\mathop{\rm Hilb}^2\mathbb{P}^2$, where the Hochschild
cohomology (or rather, Hochschild--Kostant--Rosenberg decomposition) is computed,
which corresponds to the twisted Hodge diamond using the anticanonical line bundle,
cf. [Section 4.2, 2309.06244]:

```sage
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
```

* [2309.06244] Pieter Belmans, Lie Fu, Andreas Krug: Hochschild cohomology of Hilbert schemes of points on surfaces
  [arXiv:2309.06244](https://arxiv.org/abs/2309.06244)
* [2412.09975] Lie Fu: Twisted Hodge numbers and deformation theory of Hilbert schemes of points on surfaces via Hodge modules
  [arXiv:2412.09975](https://arxiv.org/abs/2412.09975)


## Getting started

You can install it as follows:

``sage --pip install git+https://github.com/pbelmans/twisted-hodge-hilbert.git``

and then you can use

``from twisted_hilbert import *``

to use it.

The documentation with examples can be [read online](https://twisted-hilbert.ncag.info)
or as [a pdf](https://twisted-hilbert.ncag.info/documentation.pdf).


## How to cite

If you have used this code in any way,
please consider citing it as explained on [Zenodo](https://doi.org/10.5281/zenodo.14334379).
https://doi.org/10.5281/zenodo.14334379
You can choose to cite a specific version, or always the latest version.

The following BibTeX entry is a good starting point:

```bibtex
@software{hodge-diamond-cutter,
  author = {Belmans, Pieter},
  title = {Twisted Hodge numbers for Hilbert schemes of points},
  url = {https://github.com/pbelmans/twisted-hodge-hilbert},
  doi = {10.5281/zenodo.14334379},
}
```

which leads to something like

> Pieter Belmans. _Twisted Hodge numbers for Hilbert schemes of points_. doi:10.5281/zenodo.14334379. url: ht<span>tps://github.com/pbelmans/twisted-hodge-hilbert.
