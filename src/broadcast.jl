# StaticBitVectors: broadcast.jl
#
# Broadcasting for bit vectors.
#
# We bypass most of the Broadcasting infrastructure because it slows things down terribly.
# In general, using map is more performant than using broadcasting.


# Create a custom style so that we can dispatch to our custom implementation
# struct BitVectorStyle <: Broadcast.AbstractArrayStyle{1} end
# BroadcastStyle(::Type{<:BitCol}) = BitVectorStyle()

# Functions for which broadcasting can be done bitwise
# const BitwiseFun = Union{typeof(&), typeof(|), typeof(xor), typeof(nor), typeof(nand)}


#  functions that can be broadcasting by recasting in terms of bitwise functions
broadcasted(::typeof(~), a::BitCol) = bit_map(~, a)
broadcasted(::typeof(&), a::BitCol, b::BitCol) = bit_broadcast(&, a, b)
broadcasted(::typeof(|), a::BitCol, b::BitCol) = bit_broadcast(|, a, b)
broadcasted(::typeof(xor), a::BitCol, b::BitCol) = bit_broadcast(xor, a, b)
broadcasted(::typeof(nor), a::BitCol, b::BitCol) = bit_broadcast(nor, a, b)
broadcasted(::typeof(nand), a::BitCol, b::BitCol) = bit_broadcast(nand, a, b)
broadcasted(::typeof(min), a::BitCol, b::BitCol) = bit_broadcast(&, a, b)
broadcasted(::typeof(max), a::BitCol, b::BitCol) = bit_broadcast(|, a, b)
broadcasted(::typeof(*), a::BitCol, b::BitCol) = bit_broadcast(*, a, b)
broadcasted(::typeof(==), a::BitCol, b::BitCol) = bit_broadcast((x,y) -> ~xor(x,y), a, b)
broadcasted(::typeof(^), a::BitCol, b::BitCol) = bit_broadcast((x,y) -> x | ~y, a, b)
broadcasted(::typeof(>), a::BitCol, b::BitCol) = bit_broadcast((x,y) -> x & ~y, a, b)
broadcasted(::typeof(>=), a::BitCol, b::BitCol) = bit_broadcast((x,y) -> x | ~y, a, b)
broadcasted(::typeof(<), a::BitCol, b::BitCol) = bit_broadcast((x,y) -> y & ~x, a, b)
broadcasted(::typeof(<=), a::BitCol, b::BitCol) = bit_broadcast((x,y) -> y | ~x, a, b)

function bit_broadcast(f::F, a::BitCol, b::BitCol) where {F}
	if length(a) == length(b)
		bit_map(f, a, b)
	elseif length(a) == 1
		a_ = expand_bitvec(a, b)	# expanding and mapping may not be the fastest approach
		bit_map(f, a_, b)
	elseif length(b) == 1
		b_ = expand_bitvec(b, a)	# expanding and mapping may not be the fastest approach
		bit_map(f, a, b_)
	else
		error("invalid combination of lengths")
	end
end

# Expand a length-1 BitCol to the size of another BitCol 
function expand_bitvec(a::BitCol, b::BitCol{C}) where {C}
	@assert length(a) == 1
	a[1] ? trues(SBitCol{C}, length(b)) : falses(SBitCol{C}, length(b))
end


