# TODO
#	- Implment broadcasting
#
#	- Consider generalizing this to static bit arrays of any dimension.
#	- Also, consider incorporating this into the BitArray framework. This would
#		require a reworking BitArray; possibly parameterize by backing type?
#		e.g.	BitArray{N} = AbstractBitArray{N, Array{N,UInt64}}
#				SBitArray{S,N} = AbstractBitArray{N, SArray{S,N,UInt64}}
#				MBitArray{S,N} = AbstractBitArray{N, MArray{S,N,UInt64}}
#		then SBitVector{L} = SBitArray{Tuple{L}, 1}, etc.
"""
	StaticBitVectors (module)

Stack-allocated dense-packed bit vectors. (BitVectors with SVector or MVector for storage.)
"""
module StaticBitVectors

using StaticArrays
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

export AbstractBitVector, StaticBitVector, SBitVector, MBitVector
export hamming, parity, ⪯, ⪰

# TODO:	Support for other binary element types
# TODO:	Support multidimensional arrays?


## Preliminaries

const Iterable = Union{Tuple, AbstractArray, UnitRange, Base.Generator}

function allequal(it::Iterable)
	if length(it) <= 1
		return true
	end
	v = first(it)
	for v_ in it
		v_ === v || return false
	end
	return true
end


# Bitwise utilities (copied from bitarray.jl)
const _msk64 = ~UInt64(0)
@inline _div64(l) = l >> 6
@inline _mod64(l) = l & 63
@inline _msk_end(l::Int) = _msk64 >>> _mod64(-l)

# the number of chunks needed to store a given number of bits
@inline nchunks(n::Int) = _div64(n+63)

# return the chunk index and bitmask for a given bit position
@inline function chunk_index(i::Integer)
	i1 = _div64(i-1)+1
	i2 = _mod64(i-1)
	msk = UInt64(1) << i2
	return (i1, msk)
end


# Compute chunks from an iterable.
# Returns a tuple instead of an S/MVector to be acceptable by both SBitVector and MBitVector
function compute_chunks(A)
	len = length(A)
	nch = nchunks(len)

	# empty vector
	nch == 0 && return MVector{0,UInt64}()

	chunks = MVector{nch,UInt64}(undef)
	itr = iterate(A)
	@inbounds begin
		for ich = 1:nch-1
			ch = UInt64(0)
			for j = 0:63
				ch |= (UInt64(convert(Bool, itr[1])) << j)
				itr = iterate(A, itr[2])
			end
			chunks[ich] = ch
		end
		ch = UInt64(0)
		for j = 0:_mod64(len-1)
			ch |= (UInt64(convert(Bool, itr[1])) << j)
			itr = iterate(A, itr[2])
		end
		chunks[nch] = ch
	end
	return Tuple(chunks)
end


## The Types and Constructors

abstract type StaticBitVector{C} <: AbstractVector{Bool} end
const AbstractBitVector = Union{BitVector, StaticBitVector}

nchunks(::StaticBitVector{C}) where {C} = C

"""
	SBitVector{C}

Immutable dense-packed bit vector.  Uses an SVector{C,UInt64} for storage.
Constructors:

	SBitVector(data::SVector{C,UInt64}, len::Int)
	SBitVector(data::NTuple{C,UInt64}, len::Int)
	SBitVector(bits::AbstractArray{Bool})
"""
struct SBitVector{C} <: StaticBitVector{C}
	chunks::SVector{C,UInt64}
	len::Int
	SBitVector(chunks::SVector{C,UInt64}, len) where {C} = SBitVector{C}(chunks, len)
	function SBitVector{C}(chunks::SVector{C,UInt64}, len::Int) where {C}
		nchunks(len) == C || error("Require length(chunks) == C >> 6")
		new{C}(chunks, len)
	end
end


"""
	MBitVector{C}

Mutable dense-packed bit vector.  Uses an MVector{C,UInt64} for storage.
Constructors:

	MBitVector(data::SVector{C,UInt64}, len::Int)
	MBitVector{C}(undef, len::Int)
	MBitVector(data::NTuple{C,UInt64}, len::Int)
	MBitVector(bits::AbstractArray{Bool})
"""
struct MBitVector{C} <: StaticBitVector{C}
	chunks::MVector{C,UInt64}
	len::Int
	MBitVector(chunks::MVector{C,UInt64}, len) where {C} = MBitVector{C}(chunks, len)
	function MBitVector{C}(chunks::MVector{C,UInt64}, len::Int) where {C}
		nchunks(len) == C || error("Require length(chunks) == C >> 6")
		new{C}(chunks, len)
	end
	MBitVector{C}(::UndefInitializer, n) where {C} = new{C}(MVector{C,UInt64}(undef), len)
end


# construct from tuples
SBitVector(t::Tuple{Vararg{UInt64}}, len) = SBitVector(SVector(t), len)
MBitVector(t::Tuple{Vararg{UInt64}}, len) = MBitVector(MVector(t), len)

# consruct from StaticBitVectors
SBitVector(b::StaticBitVector) = SBitVector(SVector(b.chunks), b.len)
MBitVector(b::StaticBitVector) = MBitVector(MVector(b.chunks), b.len)

# Construction from BitVector
SBitVector(b::BitVector) = SBitVector(tuple(b.chunks...), length(b))
MBitVector(b::BitVector) = MBitVector(tuple(b.chunks...), length(b))

# Construction from Bool
SBitVector(b::Bool) = SBitVector((UInt64(b),), 1)
MBitVector(b::Bool) = MBitVector((UInt64(b),), 1)

# Fallback constructor (from iterable)
function (::Type{T})(itr) where T<:StaticBitVector
	chunks = compute_chunks(itr)
	T(chunks, length(itr))
end

# Construct all true StaticBitVector
trues(::Type{T}, len) where {T<:StaticBitVector} = trues(T{nchunks(len)}, len)
function trues(::Type{T}, len) where {T<:StaticBitVector{C}} where {C}
	chunks = ntuple(i-> i<C ? _msk64 : _msk_end(len), Val(C))
	return constructor_type(T)(chunks, len)
end


# Construct all false StaticBitVector
falses(::Type{T}, len) where {T<:StaticBitVector} = falses(T{nchunks(len)}, len)
function falses(::Type{T}, len) where {T<:StaticBitVector{C}} where {C}
	chunks = ntuple(i-> UInt64(0), Val(C))
	return constructor_type(T)(chunks, len)
end


# sbitmask(len, I) = bitmask(SBitVector, len, I)
# mbitmask(len, I) = bitmask(MBitVector, len, I)

# function bitmask(::Type{T}, I) where {T}
# 	T()
#  # do something
# end


## Functions to determine constructor type for generic methods

# Strip type parameters from a StaticBitVector
constructor_type(::Type{<:SBitVector}) = SBitVector
constructor_type(::Type{<:MBitVector}) = MBitVector
constructor_type(b::StaticBitVector) = constructor_type(typeof(b))

promote_type_args(arg) = constructor_type(arg)
promote_type_args(args...) = promote_type(typeof(args[1]), promote_type_args(tail(args)...))

# We follow the convention of StaticArrays:  the result is immutable unless both args are mutable
promote_rule(::Type{<:StaticBitVector}, ::Type{<:StaticBitVector}) = SBitVector
promote_rule(::Type{<:MBitVector}, ::Type{<:MBitVector}) = MBitVector


## Conversions
convert(::Type{SBitVector}, x) = SBitVector(x)
convert(::Type{MBitVector}, x) = MBitVector(x)



length(b::StaticBitVector) = b.len
ndims(b::StaticBitVector) = 1
size(b::StaticBitVector) = (length(b),)
axes(b::StaticBitVector) = (Base.OneTo(length(b)),)

# # Conversion to BitVector
# convert(::Type{BitVector}, bv::StaticBitVector) = BitVector(bv)
# convert(::Type{StaticBitVector}, bv::BitVector) = StaticBitVector(bv)
#
# function BitVector(bv::StaticBitVector)
# 	bv = BitVector(undef, length(bv))
# 	bv.chunks = Vector(bv.chunks)
# 	bv
# end


##  Indexing

checkbounds(b::StaticBitVector, I...) = Base.checkbounds_indices(Bool, axes(b), I) || Base.throw_boundserror(b, I)

# getindex entry point
getindex(b::StaticBitVector, ::Colon) = constructor_type(b)(b)

@inline function getindex(b::StaticBitVector, i)
   @boundscheck checkbounds(b, i)
	_getindex(b, i)
end

# Implementation -- assumes bounds have already been checked
@inline function _getindex(b::StaticBitVector, i::Integer)
	ich, msk = chunk_index(i)
	@inbounds r = (b.chunks[ich] & msk) != 0
	return r
end

@inline _getindex(b::StaticBitVector, i::CartesianIndex{1}) = _getindex(b, i[1])

# fallback / index by an iterable
function _getindex(b::StaticBitVector, itr)
	 v = [_getindex(b, i) for i in itr]
	 constructor_type(b)(v)
end


# setindex! entry point
function setindex!(bv::MBitVector, val, i)
	@boundscheck checkbounds(bv, i)
	_setindex!(bv, val, i)
end


# Implementation -- assumes bounds have already been checked
@inline function _setindex!(bv::MBitVector, val, i::Integer)
	i1, i2 = chunk_index(i)
	msk = ~(UInt64(1) << i2)
	@inbounds bv.chunks[i1] = (bv.chunks[i1] & msk) | (Bool(val) << i2)
end


@inline _setindex!(b::StaticBitVector, val, i::CartesianIndex{1}) = _setindex!(b, val, i[1])

function _setindex!(bv::MBitVector, val, itr)
	for (iv,ib) in enumerate(itr)
		_setindex!(bv, val[iv], ib)
	end
end


# function setindex(bv::StaticBitVector, val, i)
# 	@boundscheck checkbounds(bv, i)
# 	_setindex(bv, val, i)
# end
#
#
# @inline function _setindex(bv::StaticBitVector{C}, val::Bool, i::Integer) where {C}
# 	temp = MVector{C,UInt64}(bv.chunks)
# 	i1, i2 = Base.get_chunks_id(i)
# 	msk = ~(UInt64(1) << i2)
# 	@inbounds temp[i1] = (temp[i1] & msk) | (val << i2)
# 	StaticBitVector{C}(SVector{C,UInt64}(temp))
# end
#
#
#
# @inline function _setindex(bv::StaticBitVector{C}, val::AbstractVector{Bool}, itr) where {C}
# 	temp = MVector{C,UInt64}(bv.chunks)
# 	for (iv,ib) in enumerate(itr)
# 		i1, i2 = Base.get_chunks_id(ib)
# 		msk = ~(UInt64(1) << i2)
# 		@inbounds temp[i1] = (temp[i1] & msk) | (val[iv] << i2)
# 	end
# 	StaticBitVector{C}(SVector{C,UInt64}(temp))
# end
#
#
#



##  Iteration

function iterate(bv::StaticBitVector, i::Int=0)
    i >= length(bv) && return nothing
    (bv.chunks[_div64(i)+1] & (UInt64(1)<<_mod64(i)) != 0, i+1)
end


# concatenation
vcat(a::StaticBitVector) = constructor_type(a)(a)

# Without inlining this is SLOW
# For SBitVector, would it be faster to construct using ntuple(f, ...) instead of temporary MVector?
@inline function vcat(args::StaticBitVector...)
	totlen = sum(map(a->length(a), args))
	nch = nchunks(totlen)
	chunks = MVector{nch, UInt64}(undef)
	a = args[1]
	for i = 1:length(a.chunks)
		chunks[i] = a.chunks[i]
	end
	len = length(a)

	for a in tail(args)
		ich = nchunks(len)
		off = _mod64(len)
		for j = 1:length(a.chunks)
			chunks[ich - 1 + j]  |= a.chunks[j] << off
			if ich + j <= nch
				chunks[ich + j] = a.chunks[j] >> (64-off)
			end
		end
		len += length(a)
	end
	return promote_type_args(args...)(Tuple(chunks), totlen)
end


## Basic arithmetic operations -- result in a non-dense array
-(b::StaticBitVector) = (-).(b)
*(b::StaticBitVector, x::Number) = b .* x
*(x::Number, b::StaticBitVector) = x .* b
/(b::StaticBitVector, x::Number) = b ./ x



## Whole-array operationsl
==(a::StaticBitVector, b::StaticBitVector) = (length(a) == length(b)) && (a.chunks == b.chunks)
⪯(a::StaticBitVector, b::StaticBitVector) = all_chunks(_preceq, a, b)
⪰(a::StaticBitVector, b::StaticBitVector) = all_chunks(_succeq, a, b)

function all_chunks(op, a::StaticBitVector, b::StaticBitVector)
	length(a) == length(b) || throw(DimensionMismatch("sizes of A and B must match"))
	@inbounds for i in nchunks(a)
		op(a.chunks[i], b.chunks[i]) || return false
	end
	return true
end

# (define internally to avoid type piracy)
_preceq(x::UInt64, y::UInt64) = x == (x & y)
_succeq(x::UInt64, y::UInt64) = y == (x & y)


function count(b::AbstractBitVector)
  s = 0
  chk = b.chunks
  @inbounds for i = 1:length(chk)
		s += count_ones(chk[i])
  end
  s
end

sum(b::AbstractBitVector) = count(b)

parity(b::AbstractBitVector) = isodd(count(b))


@inline function dot(x::StaticBitVector, y::StaticBitVector)
	# simplest way to mimic Array dot behavior
	length(x) == length(y) || throw(DimensionMismatch())
	s = 0
	xc = x.chunks
	yc = y.chunks
	@inbounds for i = 1:length(xc)
		 s += count_ones(xc[i] & yc[i])
	end
	s
end


"""
  hamming(x::StaticBitVector)

Hamming weight of `x`.
"""
hamming(b::AbstractBitVector) = count(b)

"""
  hamming(x::StaticBitVector, y::StaticBitVector)

Hamming distance from `x` to `y`.
"""
function hamming(x::AbstractBitVector, y::AbstractBitVector)
  length(x) == length(y) || throw(DimensionMismatch())
  s = 0
  xc = x.chunks
  yc = y.chunks
  @inbounds for i = 1:length(xc)
		s += count_ones(xc[i] ⊻ yc[i])
  end
  s
end




## bitwise operations - return a StaticBitVector
# Should we use broadcasting instead?

checklengths(a) = nothing

@inline function checklengths(a, b)
	length(a) == length(b) || throw(DimensionMismatch("sizes of A and B must match"))
	nothing
end

@inline function checklengths(args...)
	allequal(length(arg) for arg in args) || throw(DimensionMismatch("all arguments must have the same length"))
	nothing
end

# These are redundant with broadcasting, but are faster.
(~)(a::StaticBitVector) = bit_map(~, a)
(&)(a::StaticBitVector, b::StaticBitVector) = bit_map(&, a, b)
(|)(a::StaticBitVector, b::StaticBitVector) = bit_map(|, a, b)
xor(a::StaticBitVector, b::StaticBitVector) = bit_map(xor, a, b)
nor(a::StaticBitVector, b::StaticBitVector) = bit_map(nor, a, b)
nand(a::StaticBitVector, b::StaticBitVector) = bit_map(nand, a, b)


# Efficient versions of map and map! for boolean functions.
map(::Union{typeof(~), typeof(!)}, a::StaticBitVector) = bit_map(~, a)
map(::Union{typeof(&), typeof(min)}, a::StaticBitVector, b::StaticBitVector) = bit_map(&, a, b)
map(::Union{typeof(|), typeof(max)}, a::StaticBitVector, b::StaticBitVector) = bit_map(|, a, b)
map(::Union{typeof(xor), typeof(!=)}, a::StaticBitVector, b::StaticBitVector) = bit_map(xor, a, b)
map(::typeof(nor), a::StaticBitVector, b::StaticBitVector) = bit_map(nor, a, b)
map(::typeof(nand), a::StaticBitVector, b::StaticBitVector) = bit_map(nand, a, b)
map(::typeof(*), a::StaticBitVector, b::StaticBitVector) = bit_map(*, a, b)
map(::typeof(==), a::StaticBitVector, b::StaticBitVector) = bit_map((x,y) -> ~xor(x,y), a, b)
map(::typeof(^), a::StaticBitVector, b::StaticBitVector) = bit_map((x,y) -> x | ~y, a, b)
map(::typeof(>), a::StaticBitVector, b::StaticBitVector) = bit_map((x,y) -> x & ~y, a, b)
map(::typeof(>=), a::StaticBitVector, b::StaticBitVector) = bit_map((x,y) -> x | ~y, a, b)
map(::typeof(<), a::StaticBitVector, b::StaticBitVector) = bit_map((x,y) -> y & ~x, a, b)
map(::typeof(<=), a::StaticBitVector, b::StaticBitVector) = bit_map((x,y) -> y | ~x, a, b)
map(::typeof(min), args...) where T = bit_map(&, args...)
map(::typeof(max), args...) = bit_map(|, args...)

masked(x::UInt64, ich, nch, msk) = ich < nch ? x : x & msk

# For some reason, doing length-checking first and dispatching to non-checking functions was slower

# 1-ary functions
@inline function bit_map(f::F, a::StaticBitVector{C}) where {C,F}
	C==0 && return promote_type_args(a)((), 0)

	len = length(a)
	chunkf = i -> masked(f(@inbounds a.chunks[i]), i, C, _msk_end(len))
	chunks = ntuple(chunkf, Val(C))
   promote_type_args(a)(chunks, len)
end


# 2-ary functions
@inline function bit_map(f::F, a::StaticBitVector{C},  b::StaticBitVector{C}) where {F,C}
	checklengths(a, b)
	len = length(a)

	C==0 && return promote_type_args(a,b)((), 0)

	chunkf = i -> masked((@inbounds f(a.chunks[i], b.chunks[i])), i, C, _msk_end(len))
	chunks = ntuple(chunkf, Val(C))
   promote_type_args(a, b)(chunks, len)
end


# n-ary functions
@inline function bit_map(f::F, args::StaticBitVector{C}...) where {F,C}
	checklengths(args...)
	len = length(args[1])

	C==0 && return promote_type_args(a,b)((), 0)

	chunkf = i -> masked(f( (@inbounds arg.chunks[i] for arg in args)...), i, C, _msk_end(len))
	chunks = ntuple(chunkf, Val(C))
   promote_type_args(args...)(chunks, len)
end


# # These are a slightly slower for big, inferrable arrays
# # 1-ary functions
# @inline function bit_map_old(f::F, a::StaticBitVector{C}) where {C,F}
# 	temp = MVector{C,UInt64}(undef)
# #	isempty(A) && return StaticBitVector{0,1}()
# 	@inbounds for i = 1:C
#    	temp[i] = f(a.chunks[i])
#    end
#    temp[C] &= _msk_end(length(a))
#    constructor_type(a)(Tuple(temp), length(a))
# end


# # 2-ary functions
# function bit_map_old(f::F, a::StaticBitVector{C},  b::StaticBitVector{C}) where {F,C}
# 	checklengths(a, b)
# 	len = length(a)

# 	C==0 && return promote_type_args(a,b)((), 0)

# 	chunks = MVector{C,UInt64}(undef)
# 	for i = 1:C
#    	chunks[i] = f(a.chunks[i], b.chunks[i])
#    end
#    chunks[C] &= _msk_end(length(a))
#    promote_type_args(a, b)(Tuple(chunks), len)
# end

# # n-ary functions
# @inline function bit_map_old(f::F, args::StaticBitVector{C}...) where {F,C}
# 	checklengths(args...)
# 	len = length(args[1])

# 	temp = MVector{C,UInt64}(undef)
# 	for i = 1:C
#    	temp[i] = f((arg.chunks[i] for arg in args)...)
#    end
#    temp[C] &= _msk_end(len)
#    promote_type_args(args...)(Tuple(temp), len)
# end



map!(::Union{typeof(~), typeof(!)}, dest::MBitVector, a::StaticBitVector) = bit_map!(~, dest, a)
map!(::Union{typeof(&), typeof(*), typeof(min)}, dest::MBitVector, a::StaticBitVector, b::StaticBitVector) = bit_map!(&, dest, a, b)
map!(::Union{typeof(|), typeof(max)}, dest::MBitVector, a::StaticBitVector, b::StaticBitVector) = bit_map!(|, dest, a, b)
map!(::typeof(xor), dest::MBitVector, a::StaticBitVector, b::StaticBitVector) = bit_map!(xor, dest, a, b)
map!(::typeof(nor), dest::MBitVector, a::StaticBitVector, b::StaticBitVector) = bit_map!(nor, dest, a, b)
map!(::typeof(nand), dest::MBitVector, a::StaticBitVector, b::StaticBitVector) = bit_map!(nand, dest, a, b)
map!(::typeof(*), dest::MBitVector, a::StaticBitVector, b::StaticBitVector) = bit_map!(*, dest,a, b)
map!(::typeof(!=), dest::MBitVector, a::StaticBitVector, b::StaticBitVector) = bit_map!(xor, dest, a, b)
map!(::typeof(==), dest::MBitVector, a::StaticBitVector, b::StaticBitVector) = bit_map!((x,y) -> ~xor(x,y), dest, a, b)
map!(::typeof(^), dest::MBitVector, a::StaticBitVector, b::StaticBitVector) = bit_map!((x,y) -> x | ~y, dest, a, b)
map!(::typeof(>), dest::MBitVector, a::StaticBitVector, b::StaticBitVector) = bit_map!((x,y) -> x & ~y, dest, a, b)
map!(::typeof(>=), dest::MBitVector, a::StaticBitVector, b::StaticBitVector) = bit_map!((x,y) -> x | ~y, dest, a, b)
map!(::typeof(<), dest::MBitVector, a::StaticBitVector, b::StaticBitVector) = bit_map!((x,y) -> y & ~x, dest, a, b)
map!(::typeof(<=), dest::MBitVector, a::StaticBitVector, b::StaticBitVector) = bit_map!((x,y) -> y | ~x, dest, a, b)
map!(::typeof(min), dest::MBitVector, args::StaticBitVector...) where T = bit_map!(&, dest, args...)
map!(::typeof(max), dest::MBitVector, args::StaticBitVector...) = bit_map!(|, dest, args...)



# 1-ary functions
@inline function bit_map!(f::F, dest::MBitVector{C}, a::StaticBitVector{C}) where {F,C}
	checklengths(dest, a)
	chunks = dest.chunks
	@inbounds for i = 1:C
   	chunks[i] = @inbounds f(a.chunks[i])
   end
	@inbounds chunks[C] &= _msk_end(length(dest))
	dest
end


# 2-ary functions
@inline function bit_map!(f::F, dest::MBitVector{C}, a::StaticBitVector{C}, b::StaticBitVector{C}) where {F,C}
	checklengths(dest, a, b)
	chunks = dest.chunks
	@inbounds for i = 1:C
   	chunks[i] =  @inbounds f(a.chunks[i], b.chunks[i])
   end
	@inbounds chunks[C] &= _msk_end(length(dest))
	dest
end


# n-ary functions
function bit_map!(f::F, dest::MBitVector{C}, args::StaticBitVector{C}...) where {F,C}
	checklengths(dest, args...)
	chunks = dest.chunks
	@inbounds for i = 1:C
   	chunks[i] =  @inbounds f(a.chunks[i], b.chunks[i])
   end
	@inbounds chunks[C] &= _msk_end(length(dest))
	dest
end


## Broadcasting
#
# We bypass most of the Broadcasting infrastructure because it slows things down terribly.
# In general, using map is more performant than using broadcasting.

# # Create a custom style so that we can dispatch to our custom implementation
# struct BitVectorStyle <: Broadcast.AbstractArrayStyle{1} end
# BroadcastStyle(::Type{<:StaticBitVector}) = BitVectorStyle()

# Functions for which broadcasting can be done bitwise
const BitwiseFun = Union{typeof(&), typeof(|), typeof(xor), typeof(~), typeof(nor), typeof(nand)}


broadcasted(f, a::StaticBitVector) = bit_map(f, a)

#  functions that can be broadcasting by recasting in terms of bitwise functions
broadcasted(::typeof(&), a::StaticBitVector, b::StaticBitVector) = bit_broadcast(&, a, b)
broadcasted(::typeof(|), a::StaticBitVector, b::StaticBitVector) = bit_broadcast(|, a, b)
broadcasted(::typeof(xor), a::StaticBitVector, b::StaticBitVector) = bit_broadcast(xor, a, b)
broadcasted(::typeof(nor), a::StaticBitVector, b::StaticBitVector) = bit_broadcast(nor, a, b)
broadcasted(::typeof(nand), a::StaticBitVector, b::StaticBitVector) = bit_broadcast(nand, a, b)
broadcasted(::typeof(min), a::StaticBitVector, b::StaticBitVector) = bit_broadcast(&, a, b)
broadcasted(::typeof(max), a::StaticBitVector, b::StaticBitVector) = bit_broadcast(|, a, b)
broadcasted(::typeof(*), a::StaticBitVector, b::StaticBitVector) = bit_broadcast(*, a, b)
broadcasted(::typeof(==), a::StaticBitVector, b::StaticBitVector) = bit_broadcast((x,y) -> ~xor(x,y), a, b)
broadcasted(::typeof(^), a::StaticBitVector, b::StaticBitVector) = bit_broadcast((x,y) -> x | ~y, a, b)
broadcasted(::typeof(>), a::StaticBitVector, b::StaticBitVector) = bit_broadcast((x,y) -> x & ~y, a, b)
broadcasted(::typeof(>=), a::StaticBitVector, b::StaticBitVector) = bit_broadcast((x,y) -> x | ~y, a, b)
broadcasted(::typeof(<), a::StaticBitVector, b::StaticBitVector) = bit_broadcast((x,y) -> y & ~x, a, b)
broadcasted(::typeof(<=), a::StaticBitVector, b::StaticBitVector) = bit_broadcast((x,y) -> y | ~x, a, b)

function bit_broadcast(f::F, a::StaticBitVector, b::StaticBitVector) where {F<:BitwiseFun}
	if length(a) == length(b)
		bit_map(f, a, b)
	elseif length(a) == 1
		a_ = expand_bitvec(a, b)
		bit_map(f, a_, b)
	elseif length(b) == 1
		b_ = expand_bitvec(b, a)
		bit_map(f, a, b_)
	else
		error("invalid combination of lengths")
	end
end

# Expand a length-1 StaticBitVector to the size of another StaticBitVector 
function expand_bitvec(a::StaticBitVector, b::StaticBitVector{C}) where {C}
	a[1] ? trues(SBitVector{C}, length(b)) : falses(SBitVector{C}, length(b))
end



end	# module
