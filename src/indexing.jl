# StaticBitVectors: indexing.jl
#
#  Indexing 
#


#
#  BitCol indexing
#
#	v[scalar]						::Bool
#	v[collection]					::BitCol 
#  v[x, 1, 1, ...] = v[x]		::BitCol
# Anyything else should fall back to generic array methods

checkbounds(b::BitCol, I...) = Base.checkbounds_indices(Bool, axes(b), I) || Base.throw_boundserror(b, I)

# getindex 

# v[scalar]
# v[collection]
# v[x, 1, 1, ...]
@inline function getindex(b::BitCol, i, j::Integer...)
   @boundscheck checkbounds(b, i, j...)
	_getindex(b, i)
end

#  (implementation - bounds have already been checked)

@inline function _getindex(b::BitCol, i::Integer)
	ich, msk = chunk_index(i)
	@inbounds r = (chunks(b)[ich] & msk) != 0
	return r
end

@inline _getindex(b::BitCol, i::CartesianIndex{1}) = _getindex(b, i[1])

@inline _getindex(b::BitCol, ::Colon) = b

function _getindex(b::BitCol, itr)
	basetype(b)(_getindex(b, i) for i in itr)
end


# setindex!

# v[scalar]
# v[collection]
# v[x, 1, 1, ...]
@inline function setindex!(b::MBitCol, val, i, j::Integer...)
	@boundscheck checkbounds(b, i, j...)
	_setindex!(b, val, i)
end


#  (implementation - bounds have already been checked)

@inline function _setindex!(b::MBitCol, val, i::Integer)
	i1, msk, mskb = chunk_index(i, convert(Bool, val))
	@inbounds chunks(b)[i1] = (chunks(b)[i1] & ~msk) | mskb
end


@inline _setindex!(b::MBitCol, val, i::CartesianIndex{1}) = _setindex!(b, val, i[1])

@inline _setindex!(b::MBitCol, val, ::Colon) = _setindex!(b, val, ntuple(length(b)))

# v[:] = bitvector
@inline function _setindex!(b::MBitCol, v::BitVec, ::Colon)
	if length(v) == length(b)
		for i in 1:nchunks(b)
			chunks(b)[i] = chunks(v)[i]
		end
	else
		throw(DimensionMismatch)
	end
end

# fallback - iterable
@inline function _setindex!(bv::MBitCol, val, itr)
	for (iv,ib) in enumerate(itr)
		_setindex!(bv, val[iv], ib)
	end
end


# function setindex(bv::BitCol, val, i)
# 	@boundscheck checkbounds(bv, i)
# 	_setindex(bv, val, i)
# end
#
#
# @inline function _setindex(bv::BitCol{C}, val::Bool, i::Integer) where {C}
# 	temp = MVector{C,UInt64}(chunks(bv))
# 	i1, i2 = Base.get_chunks_id(i)
# 	msk = ~(UInt64(1) << i2)
# 	@inbounds temp[i1] = (temp[i1] & msk) | (val << i2)
# 	BitCol{C}(SVector{C,UInt64}(temp))
# end
#
#
#
# @inline function _setindex(bv::BitCol{C}, val::AbstractVector{Bool}, itr) where {C}
# 	temp = MVector{C,UInt64}(chunks(bv))
# 	for (iv,ib) in enumerate(itr)
# 		i1, i2 = Base.get_chunks_id(ib)
# 		msk = ~(UInt64(1) << i2)
# 		@inbounds temp[i1] = (temp[i1] & msk) | (val[iv] << i2)
# 	end
# 	BitCol{C}(SVector{C,UInt64}(temp))
# end
#
#
#



# BitRow indexing
#
#	v[scalar]							::Bool
#	v[collection]						::BitRow
#	v[1, scalar, 1, 1, ...]			::Bool
#  v[1, collection, 1, 1, ...]	::BitRow
# Anyything else should fall back to generic array methods
#
# Unlike other adjoints, indexing with a single argument returns an adjoint

checkbounds(b::BitRow, I) = Base.checkbounds_indices(Bool, axes(b.parent), (I,)) || Base.throw_boundserror(b, I)
function checkbounds(b::BitRow, I, J...)
	i_ok = Base.checkbounds_indices(Bool, (Base.OneTo(1),), (I,))
	j_ok = Base.checkbounds_indices(Bool, axes(b.parent), J)
	(i_ok && j_ok) || Base.throw_boundserror(b, (I, J...))
end



# getindex - entry

# v[scalar]
# v[collection]
@inline function getindex(b::BitRow, i)
   @boundscheck  checkbounds(b, i)
	_getindex(b, i)
end

# needed to avoid ambiguity with getindex(::Adjoint, i::Int64) defined in Base
@inline function getindex(b::BitRow, i::Int64)
   @boundscheck checkbounds(b, i)
	_getindex(b, i)
end

# v[1, scalar, 1, 1, ...]
# v[1, collection, 1, 1, ...]
@inline function getindex(b::BitRow, i::Integer, j, k::Integer...)
   @boundscheck checkbounds(b, i, j, k...)
	_getindex(b.parent, j)
end

# getindex - implementation

@inline function _getindex(b::BitRow, i::Integer)
	ich, msk = chunk_index(i)
	@inbounds r = (chunks(b)[ich] & msk) != 0
	return r
end

@inline _getindex(b::BitRow, i::CartesianIndex{1}) = _getindex(b, i[1])

@inline _getindex(b::BitRow, ::Colon) = b

# fallback - iterable
function _getindex(b::BitRow, itr)
	basetype(b)(_getindex(b, i) for i in itr)
end





#
# Row vector setindex!
#


# setindex! entry point
@inline function setindex!(bv::MBitRow, val, i)
	@boundscheck checkbounds(bv, i)
	_setindex!(bv, val, i)
end

# needed to avoid method ambiguity with Base.Adjoint
@inline function setindex!(bv::MBitRow, val, i::Int64)
	@boundscheck checkbounds(bv.parent, i)
	_setindex!(bv, val, i)
end


# Implementation -- assumes bounds have already been checked
@inline function _setindex!(b::MBitRow, val, i::Integer)
	i1, msk, mskb = chunk_index(i, convert(Bool, val))
	@inbounds chunks(b)[i1] = (chunks(b)[i1] & ~msk) | mskb
end


@inline _setindex!(b::MBitRow, val, i::CartesianIndex{1}) = _setindex!(b, val, i[1])

function _setindex!(bv::MBitRow, val, itr)
	for (iv,ib) in enumerate(itr)
		_setindex!(bv, val[iv], ib)
	end
end


##  Iteration

function iterate(b::BitVec, i::Int=0)
    i >= length(b) && return nothing
    (chunks(b)[_div64(i)+1] & (UInt64(1)<<_mod64(i)) != 0, i+1)
end


# concatenation
vcat(a::BitCol) = a

function vcat(args::BitCol...)
	(chunks, len) = catchunks(args...)
	resulttype = basetype(promote_typeof(args...))
	return resulttype(chunkstype(resulttype)(chunks), len)
end


# Base already defines hcat for adjoint (row) vectors, so we don't have to
