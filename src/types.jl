# StaticBitVectors: types.jl
#
# The Types and Constructors
#

# Copied from TypeTools.jl
# basetype(dt::DataType) = dt.name.wrapper
# basetype(ut::UnionAll) = basetype(ut.body)

const Iterable = Union{AbstractArray, Tuple, Base.Generator}

"""
	BitCol{VT <: AbstractVector{UInt64}}

Column bit vector stored as dense-packed bits in an array of type `VT`.
Three particular subtypes are defined for common use:

   SBitCol{C} = BitCol{SVector{C,UInt64}}		# Static length, immutable
   MBitCol{C} = BitCol{MVector{C,UInt64}}		# Static length, mutable
   VBitCol    = BitCol{Vector{UInt64}}			# Variable length, mutable

Constructors:

	BitCol(chunks::VT, len::Int)
	BitCol{VT}(chunks::VT, len::Int)
	BitCol{VT}(::Iterable{Bool})
	BitCol{VT}(::Union(BitCol, BitRow, BitVector)

See also [`BitRow`](@ref).
"""
struct BitCol{VT<:AbstractVector{UInt64}} <: AbstractVector{Bool}
	chunks::VT
	len::Int
	# default constructor
	function BitCol{VT}(chunks::VT, len) where {VT<:AbstractVector{UInt64}}
		nchunks(len) == length(chunks) || error("Require length(chunks) == len >> 6")
		new{VT}(chunks, len)
	end
	# type-inferred constructor
	BitCol(chunks::VT, len) where {VT<:AbstractVector{UInt64}} = BitCol{VT}(chunks, len)

end



# For interop with Base's BitVector 
const AnyBitVector = Union{BitVector, BitCol}

# Standard subtypes
const SBitCol{C} = BitCol{SVector{C,UInt64}}
const MBitCol{C} = BitCol{MVector{C,UInt64}}
const SMBitCol{C} = Union{SBitCol{C}, MBitCol{C}}
const VBitCol = BitCol{Vector{UInt64}}

# type-inferred
SBitCol(ch::SVector{C,UInt64}, len) where {C} = SBitCol{C}(ch, len)
MBitCol(ch::MVector{C,UInt64}, len) where {C} = MBitCol{C}(ch, len)

"""
	BitRow{VT <: AbstractVector{Bool}}

Row bit vector stored as dense-packed bits in an array of type `VT`.
Alias for Adjoint{Bool, BitCol{VT}}.

BitRow objects can be constructed by taking the adjoint of a BitCol, or
using any of the same signatures as a BitCol constructor.

See also [[`BitCol`](@ref).
"""
const BitRow{VT} = Adjoint{Bool, BitCol{VT}}
const SBitRow{C} = BitRow{SVector{C,UInt64}}
const MBitRow{C} = BitRow{MVector{C,UInt64}}
const SMBitRow{C} = Union{SBitRow{C}, MBitRow{C}}
const VBitRow = BitRow{Vector{UInt64}}

# Core constructors
Adjoint{Bool, BT}(b::BT) where {BT<:BitCol} = invoke(Adjoint{Bool, BT}, Tuple{Any}, b)

# Construct directly from chunks
BitRow{VT}(ch::VT, len) where {VT <: AbstractVector{UInt64}} = BitCol{VT}(ch, len)'
# type inferred
BitRow(ch::VT, len) where {VT <: AbstractVector{UInt64}} = BitCol{VT}(ch, len)'
SBitRow(ch::SVector{C,UInt64}, len) where {C} = SBitCol{C}(ch, len)'
MBitRow(ch::MVector{C,UInt64}, len) where {C} = MBitCol{C}(ch, len)'



# Various useful type groups
const BitVec{VT} = Union{BitCol{VT}, BitRow{VT}}
const SBitVec{C} = Union{SBitCol{C}, SBitRow{C}}
const MBitVec{C} = Union{MBitCol{C}, MBitRow{C}}
const SMBitVec{C} = Union{SBitVec{C}, MBitVec{C}}
const VBitVec = Union{VBitCol, VBitRow}

# const StaticBitVec{C} = BitVec{StaticVector{C}}

# Type helpers

# Extract the BitCol type underlying a BitRow
coltype(x) = coltype(typeof(x))
coltype(::Type{BitRow{VT}}) where {VT<:AbstractVector{UInt64}} = BitCol{VT}
# handle UnionAlls over BitRow subtypes
coltype(::Type{BR}) where {BR<:BitRow} = UnionAll(BR.var, BR.body.parameters[2])

# Extrac the type of chunks
chunkstype(t::BitVec) = chunkstype(typeof(t))
chunkstype(::Type{BitCol{VT}}) where {VT<:AbstractVector{UInt64}} = VT
chunkstype(::Type{BitRow{VT}}) where {VT<:AbstractVector{UInt64}} = VT
# handle UnionAlls
chunkstype(::Type{BT}) where {BT<:BitCol} = UnionAll(BT.var, BT.body.parameters[1])
chunkstype(::Type{BT}) where {BT<:BitRow} = chunkstype(UnionAll(BT.var, BT.body.parameters[2]))

# Strip type parameters from a BitCol
basetype(b::BitVec) = basetype(typeof(b))
basetype(::Type{<:SBitCol}) = SBitCol
basetype(::Type{<:MBitCol}) = MBitCol
basetype(::Type{<:SBitRow}) = SBitRow
basetype(::Type{<:MBitRow}) = MBitRow
basetype(::Type{<:VBitCol}) = VBitCol
basetype(::Type{<:VBitRow}) = VBitRow



# Property access
chunks(bv::BitCol) = bv.chunks
chunks(av::BitRow) = av.parent.chunks
chunks(bv::BitVector) = bv.chunks	# for ease of interop with BitVector

nchunks(::SBitVec{C}) where {C} = C
nchunks(::MBitVec{C}) where {C} = C
nchunks(bv::BitVec) = length(chunks(bv))
nchunks(bv::BitVector) = length(chunks(bv))

length(bv::BitCol) = bv.len
length(av::BitRow) = av.parent.len

size(bv::BitCol) = (length(bv),)
size(bv::BitRow) = (1,length(bv))

axes(b::BitCol) = (Base.OneTo(length(b)),)
axes(b::BitRow) = (Base.OneTo(1),Base.OneTo(length(b)),)


#
# Constructors
#

# Helper constructor - construct from function that returns chunk elements
# @inline function SBitCol{C}(f::F, len) where {C,F<:Function}
# 	chunks_ = SVector(ntuple(f, Val(C)))
# 	BitCol(chunks_, len)
# end

# @inline function MBitCol{C}(f::F, len) where {C,F<:Function}
# 	chunks_ = MVector(ntuple(f, Val(C)))
# 	BitCol(chunks_, len)
# end

# @inline function BitCol{VT}(f::F, len) where {VT,F<:Function}
# 	chunks_ = UInt64[f(i) for i in 1:nchunks(len)]
# 	# chunks_ = Vector{UInt64}(undef, nchunks(len))
# 	# for i in 1:nchunks(len)
# 	# 	@inbounds chunks_[i] = f(i)
# 	# end
# 	BitCol(chunks_, len)
# end

# @inline function SBitRow{C}(f::F, len) where {C,F<:Function}
# 	chunks_ = SVector(ntuple(f, Val(C)))
# 	b = BitCol(chunks_, len)'
# 	# Adjoint{Bool, typeof(b)}(b)
# end

# @inline function MBitRow{C}(f::F, len) where {C,F<:Function}
# 	chunks_ = MVector(ntuple(f, Val(C)))
# 	BitCol(chunks_, len)'
# end

# @inline function BitRow{VT}(f::F, len) where {VT,F<:Function}
# 	chunks_ = UInt64[f(i) for i in 1:nchunks(len)]
# 	BitCol(chunks_, len)'
# end


# MBitCol(t::Tuple{Vararg{UInt64}}, len) = BitCol(MVector(t), len)
# (bt::Type{<:BitRow})(t::Tuple{Vararg{UInt64}}, len) = coltype(bt)(t, len)'


# failed attempts
# SBitCol(t::Tuple{Vararg{UInt64}}, len) = BitCol(SVector(t), len)
# MBitCol(t::Tuple{Vararg{UInt64}}, len) = BitCol(MVector(t), len)
# VBitCol(t::Tuple{Vararg{UInt64}}, len) = BitCol(collect(t), len)

#(::Type{BitCol{VT{C,UInt}} where {C}})() where {VT<:SVector} = println("ok")
# (::Type{BT})() where {BT <: BitCol{SVector{C,UInt}} where {C}} = println("ok")
# (::Type{BT})() where {BT <: BitCol{VT}} where {VT} = println("It works! $VT")	 # doesn't match SBitCol
# (::Type{BT})() where {BT <: BitCol{VT} where {VT}} = println("It works! $VT")	 # VT not available in function body
# (::Type{BT})() where {BT <: BitCol{SV{C}} where {C}} where {SV} = println("It works! $SV")	


# construct from Tuples
SBitCol(t::NTuple{C,UInt64}, len) where {C} = BitCol(SVector{C,UInt64}(t), len) 
MBitCol(t::NTuple{C,UInt64}, len) where {C} = BitCol(MVector{C,UInt64}(t), len) 
SBitCol{C}(t::NTuple{C,UInt64}, len) where {C} = BitCol(SVector{C,UInt64}(t), len) 
MBitCol{C}(t::NTuple{C,UInt64}, len) where {C} = BitCol(MVector{C,UInt64}(t), len) 

SBitRow(t::NTuple{C,UInt64}, len) where {C} = BitCol(SVector{C,UInt64}(t), len)'
MBitRow(t::NTuple{C,UInt64}, len) where {C} = BitCol(MVector{C,UInt64}(t), len)'
SBitRow{C}(t::NTuple{C,UInt64}, len) where {C} = BitCol(SVector{C,UInt64}(t), len)'
MBitRow{C}(t::NTuple{C,UInt64}, len) where {C} = BitCol(MVector{C,UInt64}(t), len)'



# consruct from BitVectors
SBitCol(bv::AnyBitVector) = SBitCol(SVector{nchunks(bv)}(chunks(bv)), length(bv))
MBitCol(bv::AnyBitVector) = MBitCol(MVector{nchunks(bv)}(chunks(bv)), length(bv))
VBitCol(bv::AnyBitVector) = VBitCol(collect(chunks(bv)), length(bv))
(BT::Type{<:BitRow})(bv::AnyBitVector, len) = coltype(BT)(bv, len)'


# # Construction from a Bool
(BT::Type{<:BitCol})(b::Bool) = BitCol(chunkstype(BT)(UInt64(b)), 1)
VBitCol(b::Bool) = BitCol(UInt64[b], 1)
(BT::Type{<:BitRow})(b::Bool) = coltype(BT)(b)'


# Construct from iterable
function (BT::Type{<:SMBitCol})(itr::Iterable)
	chunks = chunkstype(BT)(Tuple(compute_chunks(itr)))
	BitCol(chunks, length(itr))
end

#
function (BT::Type{<:BitCol})(itr::Iterable)
	chunks = chunkstype(BT)(compute_chunks(itr))
	BitCol(chunks, length(itr))
end

# associated BitRow constructor
(BT::Type{<:BitRow})(itr::Iterable) = coltype(BT)(itr)'



# Construct all true BitCol
trues(::Type{T}, len) where {T <: SMBitCol} = trues(T{nchunks(len)}, len)
function trues(::Type{T}, len) where {T<:SMBitCol{C}} where {C}
	chunkfun = i-> i<nchunks(len) ? _msk64 : _msk_end(len)
	chunks_ = chunkstype(T)(ntuple(chunkfun, C))
	T(chunks_, len)
end

function trues(::Type{T}, len) where {T<:BitCol}
	chunks_ = fill(_msk64, nchunks(len))
	T(chunks_, len)
end

# BitRow
trues(::Type{T}, len) where {T<:BitRow} = trues(coltype(T), len)'



# Construct all false BitCol
falses(::Type{T}, len) where {T <: SMBitCol} = falses(T{nchunks(len)}, len)
function falses(::Type{T}, len) where {T<:SMBitCol{C}} where {C}
	chunkfun = i-> UInt64(0)
	chunks_ = chunkstype(T)(ntuple(chunkfun, C))
	T(chunks_, len)
end

function falses(::Type{T}, len) where {T<:BitCol}
	chunks_ = fill(UInt64(0), nchunks(len))
	T(chunks_, len)
end

# BitRow
falses(::Type{T}, len) where {T<:BitRow} = falses(coltype(T), len)'




#
# Promotion rules
#

promote_typeof(arg) = typeof(arg)
promote_typeof(args...) = promote_type(typeof(args[1]), promote_typeof(tail(args)...))

# We follow the convention of StaticArrays:  the result is immutable unless both args are mutable
# SBitxxx > MBitxxx > VBitxxx
# xBitCol > xBitRow

# Identical types are automatically promoted to the same type.
# These handle mixed-type arguments.  Only one order needs to be defined.
promote_rule(::Type{<:SBitCol}, ::Type{<:MBitCol}) = SBitCol
promote_rule(::Type{<:MBitCol}, ::Type{<:VBitCol}) = MBitCol
promote_rule(::Type{<:SBitCol{C}}, ::Type{<:MBitCol{C}}) where {C} = SBitCol{C}
promote_rule(::Type{<:MBitCol{C}}, ::Type{<:VBitCol}) where {C} = MBitCol{C}

promote_rule(::Type{<:SBitRow}, ::Type{<:MBitRow}) = SBitRow
promote_rule(::Type{<:MBitRow}, ::Type{<:VBitRow}) = MBitRow
promote_rule(::Type{<:SBitRow{C}}, ::Type{<:MBitRow{C}}) where {C} = SBitRow{C}
promote_rule(::Type{<:MBitRow{C}}, ::Type{<:VBitRow}) where {C} = MBitRow{C}

promote_rule(t1::Type{<:BitCol}, t2::Type{<:BitRow}) = promote_type(t1, coltype(t2))

# promote_rule(::Type{<:SMBitCol{C}}, ::Type{<:SMBitCol{C}}) where {C} = SBitCol{C}
# promote_rule(::Type{<:SMBitRow{C}}, ::Type{<:SMBitRow{C}}) where {C} = SBitRow{C}
# promote_rule(::Type{<:MBitVec{C}}, ::Type{<:MBitVec{C}}) where {C} = MBitCol{C}
# promote_rule(::Type{<:BitVec{C}}, ::Type{<:BitVec{C}}) where {C} = SBitCol{C}


# Conversions
convert(::Type{BT}, x) where {BT<:BitCol} = BT(x)
convert(::Type{BT}, x) where {BT<:BitRow} = BT(x)



# # Conversion to BitVector
# convert(::Type{BitVector}, bv::BitCol) = BitVector(bv)

# # Construct BitVector by copying chunks.
# # (Faster than BitArray's default method that copies bits individually.)
# function BitVector(bv::BitCol)
# 	bv = BitVector(undef, length(bv))
# 	chunks(bv) = Vector(chunks(bv))
# 	bv
# end

