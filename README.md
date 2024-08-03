# StaticBitVectors.jl
Bit vectors backed by `StaticVector`s. By exploiting local storage and whole-word CPU operations on dense-packed bits, bitwise operations on static bit vectors can often be performed significantly faster than on other vector-of-Bool types.

## Installation
```
pkg> add https://github.com/benninkrs/StaticBitVectors.jl
```
## Usage

This package provides the abstract type `StaticBitVector` and two concrete subtypes:
 * `SBitVector` for an immutable vector of `Bool`s 
 * `MBitVector` for a mutable vector of `Bool`s

`SBitVector` and `MBitVector` have nearly identical interfaces, supporting most vector operations that are non length-changing. (notably, insertion/deletion of elements is currently not supported). The main difference between the two is that `MBitVector` supports `setindex!` whereas `SBitVector` does not.

`StaticBitVector`s can be constructed in the standard way from just about any array or iterator whose elements are convertable to `Bool`.  The functions `trues` and `falses` are extended to take `SBitVector` or `MBitVector` as the first argument, which allows one to specify the output type.

Bitwise operations on `StaticBitVector`s (`~`, `|`, `&`, `xor`, `nor`, and `nand`) are implemented efficiently, leveraging whole-word CPU operations to act on 64 bits at a time. An even larger set of operations (including the bitwise operators, comparisons, `min`, and, `max`) are implemented efficiently for `map` and `map!`.  Broadcasting can also be used, though it is not quite as fast as using `map` or the inherently parallel bitwise operations.  The results of such operations are `StaticBitVectors` of the same type as the input.  In the case of mixed `SBitVector` and `MBitVector` inputs, the output is an `SBitVector` (which follows the convention of `StaticArrays`).  When one of the arguments is not a `StaticBitVector`, `map` and broadcasting fall back to standard element-by-element operations.

Efficient implementations of reducing operations, including `count`, `parity`, `dot`, and `hamming` (Hamming weight and Hamming distance) are also provided.
## Example
```
julia> s = SBitVector([false, true, true, false, true])
5-element SBitVector{1}:
 0
 1
 1
 0
 1

 julia> m = SBitVector((isodd(i) for i in 0:4))
5-element MBitVector{1}:
 0
 1
 0
 1
 0

julia> s[2]
true

julia> s[4:-1:2]
3-element SBitVector{1}:
 0
 1
 1

 julia> s & m
 5-element SBitVector{1}:
 0
 1
 0
 0
 0

julia> map(<=, s, m)
 5-element SBitVector{1}:
 1
 1
 0
 1
 0
```
Vertical concatenation of `StaticBitVector`s yields another `StaticBitVector`, however, horizontal concatenation yields a (non-static) `BitMatrix`.


## Performance
