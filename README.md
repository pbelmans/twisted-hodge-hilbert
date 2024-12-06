[![tests](https://github.com/pbelmans/twisted-hodge-hilbert/actions/workflows/tests.yml/badge.svg)](https://github.com/pbelmans/twisted-hodge-hilbert/actions)

# Twisted Hodge numbers for Hilbert schemes of points

This package computes twisted Hodge numbers for Hilbert schemes of points on surfaces
in terms of the cohomology of powers of line bundles on surfaces.

We implement the _twisted Hodge number conjecture_ as stated in Conjecture E of [2309.06244]
and subsequently proven in [Fu].

It reads

```math
  \sum_{n\geq 0}\sum_{p=0}^{2n}\sum_{q=0}^{2n}\mathrm{h}^{p,q}((\mathop{\rm Hilb}\nolimits^nS,{L}_n)x^py^qt^n
  =
  \prod_{k\ge 1}\prod_{p=0}^2\prod_{q=0}^2\left( 1-(-1)^{p+q}x^{p+k-1}y^{q+k-1}t^k\right)^{-(-1)^{p+q}\mathrm{h}^{p,q}(S,{L}^{\otimes k})}.
```

* [2309.06244] Pieter Belmans, Lie Fu, Andreas Krug: Hochschild cohomology of Hilbert schemes of points on surfaces
  [arXiv:2309.06244](https://arxiv.org/abs/2309.06244)
* [Fu] Lie Fu: Twisted Hodge numbers and deformation theory of Hilbert schemes of points on surfaces via Hodge modules
