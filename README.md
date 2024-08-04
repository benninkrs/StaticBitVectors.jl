# StaticBitVectors.jl
Bit vectors backed by `StaticVector`s. By exploiting local storage and whole-word CPU operations on dense-packed bits, bitwise operations on static bit vectors can often be performed significantly faster than on other vector-of-Bool types.

## Installation
```
pkg> add https://github.com/benninkrs/StaticBitVectors.jl
```
## Usage

This package provides four concrete types: 
 * `SBitCol` for an immutable column vector of `Bool`s 
 * `MBitCol` for a mutable column vector of `Bool`s
 * `SBitRow` for an immutable row vector of `Bool`s 
 * `MBitRow` for a mutable row vector of `Bool`s

 The two column types are subtypes of `AbstractVector{Bool}`, while the row types are aliases for adjoints of column types, and are subtypes of `AbstractMatrix{Bool}`.  These four types have nearly identical interfaces, supporting most array operations that are non length-changing. (notably, insertion/deletion of elements is currently not supported). The main difference is that the two mutable types (`MBitCol` and `MBitRow`) support `setindex!`, whereas the other two do not. Collectively, these four types are called `BitVec`s.

`BitVec`s can be constructed and converted from just about any array or iterator whose elements are convertable to `Bool`.  The functions `trues` and `falses` are extended to take a `BitVec` type as the first argument, which allows one to specify the output type.

`BitVecs` can be indexed by integers, Cartesian indices, and iterables. **Logical indexing is not yet supported.**

Bitwise operations on `BitVec`s (`~`, `|`, `&`, `xor`, `nor`, and `nand`) are implemented efficiently, leveraging whole-word CPU operations to act on 64 bits at a time. An even larger set of operations (including the bitwise operators, comparisons, `min`, and, `max`) are implemented efficiently for `map` and `map!`.  Broadcasting on arrays of the same shape can also be used, with speed comparable to that of `map` and the inherently parallel bitwise operations.  The results of such operations are `BitVecs` of the same type as the input (but see the next section for what happens in the case of mixed input types).

Efficient implementations of reducing operations, including `count`, `parity`, `dot`, and `hamming` (Hamming weight and Hamming distance) are also provided.

## Example
```
julia> s = SBitCol([false, true, true, false, true])
5-element SBitCol{1}:
 0
 1
 1
 0
 1

 julia> m = SBitCol((isodd(i) for i in 0:4))
5-element MBitCol{1}:
 0
 1
 0
 1
 0

julia> s[2]
true

julia> s[4:-1:2]
3-element SBitCol{1}:
 0
 1
 1

 julia> s & m
 5-element SBitCol{1}:
 0
 1
 0
 0
 0

julia> map(<=, s, m)
 5-element SBitCol{1}:
 1
 1
 0
 1
 0
```
Vertical concatenation of `BitCol`s yields another `BitCol`; likewise, , however, horizontal concatenation yields a (non-static) `BitMatrix`.

## Corner Cases

For most ways of using `BitVec`s, the expected behavior is obvious.  In a few circumstances, however, it is less clear what behavior should be expected.  Some of these ambiguous circumstances and the chosen behaviors are:

* Bitwise operations involving a combination of mutable and immutable `BitVec`s return an immutable `BitVec`. (This follows the convention of `StaticArray`s.)
* Bitwise operations involving a combination of column and row `BitVec`s return a column `BitVec`.
* Bitwise operations require `BitVecs` to have the same shape (length and orientation). 
* In contrast, `map` only requires `BitVec` arguments to have the same lenghth, not the same orientation. (This follows the behavior of `map` for standard vector types.)


## Performance
