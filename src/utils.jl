# StaticBitVectors: utils.jl
#
# Primitive utilities for all BitVec types
#


# Bitwise utilities (copied from bitarray.jl)
const _msk64 = ~UInt64(0)
@inline _div64(k) = k >> 6
@inline _mod64(k) = k & 63
@inline _msk_end(k::Int) = _msk64 >>> unsigned(_mod64(-k))  # mask to keep only the lowest k bits

# the number of chunks needed to store a given number of bits
@inline nchunks(n::Int) = _div64(n+63)


# return the chunk index and bitmask for a given bit position
@inline function chunk_index(i::Integer)
	i1 = _div64(i-1)+1
	i2 = _mod64(i-1)
	msk = UInt64(1) << i2
	return (i1, msk)
end

@inline function chunk_index(i::Integer, b::Bool)
	i1 = _div64(i-1)+1
	i2 = _mod64(i-1)
	msk = UInt64(1) << i2
	mskb = b << i2
	return (i1, msk, mskb)
end


# Compute chunks from an iterable.
function compute_chunks(A)
	len = length(A)
	nch = nchunks(len)

	# empty vector
	nch == 0 && return UInt64[]

	chunks = Vector{UInt64}(undef, nch)
	itr = iterate(A)
	@inbounds for ich = 1:nch
		ch = UInt64(0)
		maxj = ich < nch ? 63 : _mod64(len-1)
		for j = 0:maxj
			ch |= (UInt64(convert(Bool, itr[1])) << j)
			itr = iterate(A, itr[2])
		end
		chunks[ich] = ch
	end
	return chunks
end



# Uses MVector for chunks, instead of Vector.
# This is much slower. It has 11 allocations, whereas with Vector it has only 1.  Why?
# Maybe because the length is not a compiler constant
function compute_chunks_m(A)
	len = length(A)
	nch = nchunks(len)

	# empty vector
	nch == 0 && return ()

	chunks = MVector{nch,UInt64}(undef)
	itr = iterate(A)
	@inbounds for ich = 1:nch
		ch = UInt64(0)
		maxj = ich < nch ? 63 : _mod64(len-1)
		for j = 0:maxj
			ch |= (UInt64(convert(Bool, itr[1])) << j)
			itr = iterate(A, itr[2])
		end
		chunks[ich] = ch
	end
	return Tuple(chunks)
end



# this is slowest
function compute_chunks_s(A)
	len = length(A)
	nch = nchunks(len)

	# empty vector
	nch == 0 && return ()

	chunks = fill(UInt64(0), SVector{nch,UInt64})
	itr = iterate(A)
	@inbounds for ich = 1:nch
		ch = UInt64(0)
		maxj = ich < nch ? 63 : _mod64(len-1)
		for j = 0:maxj
			ch |= (UInt64(convert(Bool, itr[1])) << j)
			itr = iterate(A, itr[2])
		end
		chunks = setindex(chunks, ch, ich)
	end
	return chunks
end



# Without inlining this is SLOW
# For SBitCol, would it be faster to construct using ntuple(f, ...) instead of temporary MVector?
@inline function catchunks(args...)
	totlen = sum(map(a->length(a), args))
	nch = nchunks(totlen)
	chunks = MVector{nch, UInt64}(undef)
	a = args[1]
	for i = 1:length(a.chunks)
		chunks[i] = a.chunks[i]
	end
	cumlen = length(a)

	for a in tail(args)
		ich = nchunks(cumlen)
		off = _mod64(cumlen)
		for j = 1:length(a.chunks)
			chunks[ich - 1 + j]  |= a.chunks[j] << off
			if ich + j <= nch
				chunks[ich + j] = a.chunks[j] >> (64-off)
			end
		end
		cumlen += length(a)
	end
	return (Tuple(chunks), totlen)
end
