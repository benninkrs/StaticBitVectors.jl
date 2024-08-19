# StaticBitVectors: ops.jl
#
#  Operations on bit vectors
#



## Basic arithmetic operations -- result in a non-dense array
-(b::BitVec) = (-).(b)
*(b::BitVec, x::Number) = b .* x
*(x::Number, b::BitVec) = x .* b
/(b::BitCol, x::Number) = b ./ x



## Whole-array operationsl
==(a::BitVec, b::BitVec) = (size(a) == size(b)) && (chunks(a) == chunks(b))
⪯(a::BitVec, b::BitVec) = all_chunks(_preceq, a, b)
⪰(a::BitVec, b::BitVec) = all_chunks(_succeq, a, b)


# helper functions
function all_chunks(op, a::BitVec, b::BitVec)
	length(a) == length(b) || throw(DimensionMismatch("sizes of A and B must match"))
	@inbounds for i in nchunks(a)
		op(chunks(a)[i], chunks(b)[i]) || return false
	end
	return true
end

# (define internally to avoid type piracy)
_preceq(x::UInt64, y::UInt64) = x == (x & y)
_succeq(x::UInt64, y::UInt64) = y == (x & y)


function count(b::AnyBitVector)
  s = 0
  chk = chunks(b)
  @inbounds for i = 1:length(chk)
		s += count_ones(chk[i])
  end
  s
end

sum(b::AnyBitVector) = count(b)

parity(b::AnyBitVector) = isodd(count(b))

# Computed over the integers, not GF(2)
@inline function dot(x::BitVec, y::BitVec)
	# simplest way to mimic Array dot behavior
	length(x) == length(y) || throw(DimensionMismatch())
	s = 0
	xc = chunks(x)
	yc = chunks(y)
	@inbounds for i = 1:length(xc)
		 s += count_ones(xc[i] & yc[i])
	end
	s
end

# Leave this for a separate algebra package?
# *(a::BitRow, b::BitCol) = dot(a,b)


"""
  hamming(x::BitVec)

Hamming weight of `x`.
"""
hamming(b::AnyBitVector) = count(b)

"""
  hamming(x::BitVec, y::BitVec)

Hamming distance from `x` to `y`.
"""
function hamming(x::AnyBitVector, y::AnyBitVector)
  length(x) == length(y) || throw(DimensionMismatch())
  s = 0
  xc = chunks(x)
  yc = chunks(y)
  @inbounds for i = 1:length(xc)
		s += count_ones(xc[i] ⊻ yc[i])
  end
  s
end




## bitwise operations - return a BitCol/BitVe

checklengths(a) = nothing

@inline function checklengths(a, b)
	length(a) == length(b) || throw(DimensionMismatch("sizes of A and B must match"))
	nothing
end

@inline function checklengths(args...)
	allequal(length(arg) for arg in args) || throw(DimensionMismatch("all arguments must have the same length"))
	nothing
end

checkorientation(a::BitCol, b::BitCol) = nothing
checkorientation(a::BitRow, b::BitRow) = nothing
checkorientation(a::BitVec, b::BitVec) = throw(DimensionMismatch("Arguments must have the same ortientation"))

# These are redundant with bit_map and broadcasting, but are faster.
# Arguments must have the same size (length and orientation)
(~)(a::BitVec) = bit_map(~, a)
 (&)(a::BitVec, b::BitVec) = begin checkorientation(a,b); bit_map(&, a, b); end
(|)(a::BitVec, b::BitVec) = begin checkorientation(a,b); bit_map(|, a, b); end
xor(a::BitVec, b::BitVec) = begin checkorientation(a,b); bit_map(xor, a, b); end
nor(a::BitVec, b::BitVec) = begin checkorientation(a,b); bit_map(nor, a, b); end
nand(a::BitVec, b::BitVec) = begin checkorientation(a,b); bit_map(nand, a, b); end


# Efficient versions of map and map! for boolean functions.
# Inlining the 2-argument functions (and 2-arg bit_map below) makes them about
# 30% faster for nchunks >= 66, but could conceivably cause memory bloat.
# Inlining seems to have no effect on 1-arg map
map(::Union{typeof(~), typeof(!)}, a::BitVec) = bit_map(~, a)
map(::Union{typeof(&), typeof(min)}, a::BitVec, b::BitVec) = bit_map(&, a, b)
map(::Union{typeof(|), typeof(max)}, a::BitVec, b::BitVec) = bit_map(|, a, b)
map(::Union{typeof(xor), typeof(!=)}, a::BitVec, b::BitVec) = bit_map(xor, a, b)
map(::typeof(nor), a::BitVec, b::BitVec) = bit_map(nor, a, b)
map(::typeof(nand), a::BitVec, b::BitVec) = bit_map(nand, a, b)
map(::typeof(*), a::BitVec, b::BitVec) = bit_map(*, a, b)
map(::typeof(==), a::BitVec, b::BitVec) = bit_map((x,y) -> ~xor(x,y), a, b)
map(::typeof(^), a::BitVec, b::BitVec) = bit_map((x,y) -> x | ~y, a, b)
map(::typeof(>), a::BitVec, b::BitVec) = bit_map((x,y) -> x & ~y, a, b)
map(::typeof(>=), a::BitVec, b::BitVec) = bit_map((x,y) -> x | ~y, a, b)
map(::typeof(<), a::BitVec, b::BitVec) = bit_map((x,y) -> y & ~x, a, b)
map(::typeof(<=), a::BitVec, b::BitVec) = bit_map((x,y) -> y | ~x, a, b)
map(::typeof(min), args...) = bit_map(&, args...)
map(::typeof(max), args...) = bit_map(|, args...)


# For some reason, doing length-checking first and dispatching to non-checking functions was slower



## 1-ary functions

# Default - VBitCol, VBitRow
function bit_map(f::F, a::BitVec) where {F<:Function}
	nch = nchunks(a)
	nch==0 && return typeof(a)(UInt64[], 0)

	len = length(a)

	chunks_ = Vector{UInt64}(undef, nch)
	for i = 1:nch
	 	@inbounds chunks_[i] = f(chunks(a)[i])
	end
	@inbounds chunks_[nch] &= _msk_end(len)

	typeof(a)(chunks_, len)
end


# Specialized for static-based 
function bit_map(f::F, a::SMBitVec) where {F<:Function}
	nch = nchunks(a)
	nch==0 && return typeof(a)((), 0)

	len = length(a)
	msk = _msk_end(len)

	masked(v,i) = i<nch ? v : v & msk
	chunkfun = Base.@constprop :none @inbounds i -> masked(f(chunks(a)[i]), i)

	chunks_ = ntuple(chunkfun, Val(nch))
	typeof(a)(chunks_, len)
	# SBitCol(chunks_, len)
end



## 2-ary functions

# Default - VBitCol, VBitRow
function bit_map(f::F, a::BitVec,  b::BitVec) where {F<:Function}
	checklengths(a, b)
	nch = nchunks(a)
	outtype = promote_typeof(a,b)

	nch==0 && return outtype(UInt64[], 0)

	len = length(a)

	chunks_ = Vector{UInt64}(undef, nch)
	for i = 1:nch
		@inbounds chunks_[i] = f(chunks(a)[i], chunks(b)[i])
	end
	@inbounds chunks_[nch] &= _msk_end(len)

	outtype(chunks_, len)
end


# Specialized for static-based 
function bit_map(f::F, a::SMBitVec,  b::SMBitVec) where {F<:Function}
	checklengths(a, b)
	nch = nchunks(a)
	outtype = promote_typeof(a,b)

	nch==0 && return outtype((), 0)

	len = length(a)
	msk = _msk_end(len)

	masked(v,i) = i<nch ? v : v & msk
	chunkfun = Base.@constprop :none @inbounds i -> masked(f(chunks(a)[i], chunks(b)[i]), i)
	chunks_ = ntuple(chunkfun, Val(nch))
	outtype(chunks_, len)
end



##  n-ary functions

# Default - VBitCol, VBitRow
function bit_map(f::F, args::BitVec...) where {F<:Function}
	checklengths(args...)
	nch = nchunks(args[1])
	outtype = promote_typeof(args...)
	
	nch==0 && return outtype(UInt64[], 0)

	len = length(args[1])

	chunks_ = Vector{UInt64}(undef, nch)
	for i = 1:nch
		@inbounds chunks_[i] = f((chunks(arg)[i] for arg in args)...)
	end
	@inbounds chunks_[nch] &= _msk_end(len)
	outtype(chunks_, len)
end


# Specialized for static-based 
function bit_map(f::F, args::SMBitVec...) where {F<:Function}
	checklengths(args...)
	nch = nchunks(args[1])
	outtype = promote_typeof(args...)
	nch==0 && return outtype((), 0)

	len = length(args[1])
	msk = _msk_end(len)

	masked(v,i) = i<nch ? v : v & msk
	chunkfun = Base.@constprop :none @inbounds i -> masked(f((chunks(arg)[i] for arg in args)...), i)
   chunks_ = ntuple(chunkfun, Val(nch))
	outtype(chunks_, len)
end


# # These are a slightly slower for big, inferrable arrays
# # 1-ary functions
# @inline function bit_map_old(f::F, a::BitCol{C}) where {C,F}
# 	temp = MVector{C,UInt64}(undef)
# #	isempty(A) && return BitCol{0,1}()
# 	@inbounds for i = 1:C
#    	temp[i] = f(chunks(a)[i])
#    end
#    temp[C] &= _msk_end(length(a))
#    constructor_type(a)(Tuple(temp), length(a))
# end


# # 2-ary functions
# function bit_map_old(f::F, a::BitCol{C},  b::BitCol{C}) where {F,C}
# 	# checklengths(a, b)
# 	len = length(a)

# 	C==0 && return promote_type_args(a,b)((), 0)

# 	chunks = MVector{C,UInt64}(undef)
# 	@inbounds for i = 1:C
#    	chunks[i] = f(chunks(a)[i], chunks(b)[i])
#    end
#    chunks[C] &= _msk_end(len)
#    promote_type_args(a, b)(Tuple(chunks), len)
# end

# # n-ary functions
# @inline function bit_map_old(f::F, args::BitCol{C}...) where {F,C}
# 	checklengths(args...)
# 	len = length(args[1])

# 	temp = MVector{C,UInt64}(undef)
# 	for i = 1:C
#    	temp[i] = f((chunks(arg)[i] for arg in args)...)
#    end
#    temp[C] &= _msk_end(len)
#    promote_type_args(args...)(Tuple(temp), len)
# end


# It is ok for dest to alias one of the input args.
map!(::Union{typeof(~), typeof(!)}, dest::MBitVec, a::BitVec) = bit_map!(~, dest, a)
map!(::Union{typeof(&), typeof(*), typeof(min)}, dest::MBitVec, a::BitVec, b::BitVec) = bit_map!(&, dest, a, b)
map!(::Union{typeof(|), typeof(max)}, dest::MBitCol, a::BitVec, b::BitVec) = bit_map!(|, dest, a, b)
map!(::typeof(xor), dest::MBitVec, a::BitVec, b::BitVec) = bit_map!(xor, dest, a, b)
map!(::typeof(nor), dest::MBitVec, a::BitVec, b::BitVec) = bit_map!(nor, dest, a, b)
map!(::typeof(nand), dest::MBitVec, a::BitVec, b::BitVec) = bit_map!(nand, dest, a, b)
map!(::typeof(*), dest::MBitVec, a::BitVec, b::BitVec) = bit_map!(*, dest,a, b)
map!(::typeof(!=), dest::MBitVec, a::BitVec, b::BitVec) = bit_map!(xor, dest, a, b)
map!(::typeof(==), dest::MBitVec, a::BitVec, b::BitVec) = bit_map!((x,y) -> ~xor(x,y), dest, a, b)
map!(::typeof(^), dest::MBitVec, a::BitVec, b::BitVec) = bit_map!((x,y) -> x | ~y, dest, a, b)
map!(::typeof(>), dest::MBitVec, a::BitVec, b::BitVec) = bit_map!((x,y) -> x & ~y, dest, a, b)
map!(::typeof(>=), dest::MBitVec, a::BitVec, b::BitVec) = bit_map!((x,y) -> x | ~y, dest, a, b)
map!(::typeof(<), dest::MBitVec, a::BitVec, b::BitVec) = bit_map!((x,y) -> y & ~x, dest, a, b)
map!(::typeof(<=), dest::MBitVec, a::BitVec, b::BitVec) = bit_map!((x,y) -> y | ~x, dest, a, b)
map!(::typeof(min), dest::MBitVec, args::BitVec...) = bit_map!(&, dest, args...)
map!(::typeof(max), dest::MBitVec, args::BitVec...) = bit_map!(|, dest, args...)


# TODO:  Are these slower than bit_map?

# 1-ary functions
function bit_map!(f::F, dest::Union{MBitVec, VBitVec}, a::BitVec) where {F<:Function}
	checklengths(dest, a)
	nch = nchunks(dest)
	@inbounds for i = 1:nch
   	chunks(dest)[i] = @inbounds f(chunks(a)[i])
   end
	@inbounds chunks(dest)[nch] &= _msk_end(length(dest))
	dest
end


# 2-ary functions
function bit_map!(f::F, dest::Union{MBitVec, VBitVec}, a::BitVec, b::BitVec) where {F<:Function}
	checklengths(dest, a, b)
	nch = nchunks(dest)
	@inbounds for i = 1:nch
   	chunks(dest)[i] = f(chunks(a)[i], chunks(b)[i])
   end
	@inbounds chunks(dest)[nch] &= _msk_end(length(dest))
	dest
end




# n-ary functions
function bit_map!(f::F, dest::Union{MBitVec, VBitVec}, args::BitVec...) where {F<:Function}
	checklengths(dest, args...)
	nch = nchunks(dest)
	@inbounds for i = 1:nch
   	chunks(dest)[i] =  @inbounds f((chunks(arg)[i] for arg in args)...)
   end
	@inbounds chunks(dest)[nch] &= _msk_end(length(dest))
	dest
end

