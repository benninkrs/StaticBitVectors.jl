# StaticBitVectors.jl
#
# TODO:
#	- Implement broadcasting for row vectors and mixed row-column vectors
#	- Support logical indexing
#	- Consider incorporating this into the BitArray framework. This would
#		require a reworking BitArray; possibly parameterize by backing type?
#		e.g.	BitArray{N} = AbstractBitArray{N, Array{N,UInt64}}
#				SBitArray{S,N} = AbstractBitArray{N, SArray{S,N,UInt64}}
#				MBitArray{S,N} = AbstractBitArray{N, MArray{S,N,UInt64}}
#		then SBitCol{L} = SBitArray{Tuple{L}, 1}, etc.
# - Support for other binary element types?
# - Support multidimensional arrays?
# - Support for Int64[] backing type?

"""
	StaticBitVectors (module)

Dense-packed bit vectors backed by `StaticArray`s.
"""
module StaticBitVectors

using StaticArrays
using LinearAlgebra: Adjoint
using Base: tail

import Base: convert, promote_type, promote_rule
import Base: string, print, show, display
import Base: length, size, checkbounds, ndims, axes, getindex, setindex!, iterate
import Base: trues, falses, vcat
import Base: map, map!		#, bit_map!
import Base: count, sum
import Base: +, -, *, /, ==, ~, &, |, xor, nor, nand, cmp
import LinearAlgebra: dot
import Base.Broadcast: broadcastable, BroadcastStyle, Broadcasted, broadcasted

export AbstractBitVector, SBitCol, MBitCol, SBitRow, MBitRow, BitCol, BitRow, BitVec
export dot, hamming, parity, ⪯, ⪰


include("utils.jl")
include("types.jl")
include("indexing.jl")
include("ops.jl")
include("broadcast.jl")


end	# module
