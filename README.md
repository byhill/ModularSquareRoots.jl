# ModularSquareRoots

[![Build Status](https://ci.appveyor.com/api/projects/status/github/byhill/ModularSquareRoots.jl?svg=true)](https://ci.appveyor.com/project/byhill/ModularSquareRoots-jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

## Overview
This module provides support for finding modular square-roots.
In particular, for a given integer $n$ and modulus $m$,
this module provides support for solving the congruence $x^2 \equiv n \pmod m$.

To do this, you can use the function `sqrtmod(n, m)`.

```julia
julia> using ModularSquareRoots

julia> sqrtmod(4, 5)
2-element Vector{Int64}:
 3
 2

julia> all(powermod(x, 2, 5) == 4 for x in sqrtmod(4, 5))
true

julia> sqrtmod(1240, 289032)
8-element Vector{Int64}:
 107056
 251572
  10712
 155228
 278320
 133804
 181976
  37460

julia> all(powermod(x, 2, 289032) == 1240 for x in sqrtmod(1240, 289032))
true

julia> sqrtmod(23, 200)
Int64[]

julia> !any(powermod(x, 2, 200) == 23 for x in sqrtmod(23, 200))
true
```
