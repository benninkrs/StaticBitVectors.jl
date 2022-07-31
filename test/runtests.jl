using StaticBitVectors
using Test


@info "Testing Constructors"
# Construction from explicit Bool array
s4 = SBitVector([false, false, true, true])
m4 = MBitVector([false, true, false, true])
@test typeof(s4) == SBitVector{1}
@test typeof(m4) == MBitVector{1}

# Construction from 0,1-integer array
@test SBitVector([0,0,1,1]) == s4
@test MBitVector([0,1,0,1]) == m4

# Construction from other BitVectors
b4 = BitVector([true,false,false,true])
@test typeof(SBitVector(s4)) == SBitVector{1}
@test typeof(SBitVector(s4)) == SBitVector{1}
@test typeof(SBitVector(b4)) == SBitVector{1}
@test typeof(MBitVector(s4)) == MBitVector{1}
@test typeof(MBitVector(s4)) == MBitVector{1}
@test typeof(MBitVector(b4)) == MBitVector{1}

m4_ = m4;
m4_mcopy = MBitVector(m4)
m4_scopy=  SBitVector(m4);
@test m4_ === m4
@test m4_mcopy !== m4
@test m4_scopy !== m4


# cconstruction from general iterable
@test typeof(SBitVector(isodd(i) for i in 0:255)) == SBitVector{4}
@test typeof(MBitVector(isodd(i) for i in 0:255)) == MBitVector{4}

# construction from Bool
@test SBitVector(true) == SBitVector([true])
@test MBitVector(true) == MBitVector([true])

# special vectors
@test trues(SBitVector,6) == SBitVector([true, true, true, true, true, true])
@test trues(MBitVector,6) == MBitVector([true, true, true, true, true, true])
@test falses(SBitVector,5) == SBitVector([false, false, false, false, false])
@test falses(MBitVector,5) == MBitVector([false, false, false, false, false])


@info "Testing length"
@test length(SBitVector([0,1,1,0,1,0,0,1,0])) == 9


@info "Testing indexing"
# getting
@test s4[2] == false
@test s4[3] == true
@test m4[1] == false
@test m4[4] == true
@test s4[[2,3,1]] == SBitVector([false, true, false])
@test m4[4:-1:2] == MBitVector([true,false, true])

# setting
m4_[3] = true;
@test m4_ == MBitVector([false, true, true, true])
@test m4 == MBitVector([false, true, true, true])
@test m4_mcopy == MBitVector([false, true, false, true])
@test m4_scopy == SBitVector([false, true, false, true])
m4[3] = false;
@test m4 == MBitVector([false, true, false, true])
@test m4_ == MBitVector([false, true, false, true])


@info "Tesing iteration"
@test s4 == SBitVector((bit for bit in s4))


@info "Tesing vcat"
@test [s4; m4] == SBitVector([0, 0, 1, 1, 0, 1, 0, 1])


@info "Tesing widening arithmetic"
@test -s4 == [0, 0, -1, -1]
@test 1.5*m4 == [0, 1.5, 0, 1.5]
@test s4/2 == [0, 0, 0.5, 0.5]


@info "Testing hamming and parity"
v = SBitVector([1,0,0,1,1,0,1,0,0,1,0,1,1,1,0,1])
@test count(v) == 9
@test sum(v) == 9
@test parity(v) == true
@test hamming(v) == 9
@test hamming(v, MBitVector([1,0,1,1,0,0,1,0,1,1,0,1,1,0,1,0])) == 6