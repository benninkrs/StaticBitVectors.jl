using StaticBitVectors
using Test


@info "Testing Constructors"
# Construction from explicit Bool array
sc4 = SBitCol([false, false, true, true])
mc4 = MBitCol([false, true, false, true])
sr4 = SBitRow([true, false, false, true])
mr4 = MBitRow([true, false, true, false])
b4 = BitVector([false, true, true, false])
@test typeof(sc4) == SBitCol{1}
@test typeof(mc4) == MBitCol{1}
@test typeof(sr4) == SBitRow{1}
@test typeof(mr4) == MBitRow{1}


# Construction from 0,1-integer array
@test SBitCol([0,0,1,1]) == sc4
@test MBitCol([0,1,0,1]) == mc4
@test SBitRow([1,0,0,1]) == sr4
@test MBitRow([1,0,1,0]) == mr4

# Construction from other BitVectors
@test typeof(SBitCol(sc4)) == SBitCol{1}
@test typeof(SBitCol(mr4)) == SBitCol{1}
@test typeof(SBitCol(b4)) == SBitCol{1}
@test typeof(MBitCol(sc4)) == MBitCol{1}
@test typeof(MBitCol(sc4)) == MBitCol{1}
@test typeof(MBitCol(b4)) == MBitCol{1}

mc4_ = mc4;
mc4_mcopy = MBitCol(mc4)
mc4_scopy=  SBitCol(mc4);
@test mc4_ === mc4
@test mc4_mcopy !== mc4
@test mc4_scopy !== mc4


# cconstruction from general iterable
@test typeof(SBitCol(isodd(i) for i in 0:255)) == SBitCol{4}
@test typeof(MBitCol(isodd(i) for i in 0:255)) == MBitCol{4}

# construction from Bool
@test SBitCol(true) == SBitCol([true])
@test MBitCol(true) == MBitCol([true])

# special vectors
@test trues(SBitCol,6) == SBitCol([true, true, true, true, true, true])
@test trues(MBitCol,6) == MBitCol([true, true, true, true, true, true])
@test falses(SBitCol,5) == SBitCol([false, false, false, false, false])
@test falses(MBitCol,5) == MBitCol([false, false, false, false, false])
@test trues(SBitRow,6) == SBitRow([true, true, true, true, true, true])
@test trues(MBitRow,6) == MBitRow([true, true, true, true, true, true])
@test falses(SBitRow,5) == SBitRow([false, false, false, false, false])
@test falses(MBitRow,5) == MBitRow([false, false, false, false, false])


@info "Testing length"
@test length(SBitCol([0,1,1,0,1,0,0,1,0])) == 9


@info "Testing indexing"
# getting
@test sc4[2] == false
@test sc4[3] == true
@test mc4[1] == false
@test mc4[4] == true
@test sc4[[2,3,1]] == SBitCol([false, true, false])
@test mc4[4:-1:2] == MBitCol([true,false, true])

# setting
mc4_[3] = true;
@test mc4_ == MBitCol([false, true, true, true])
@test mc4 == MBitCol([false, true, true, true])
@test mc4_mcopy == MBitCol([false, true, false, true])
@test mc4_scopy == SBitCol([false, true, false, true])
mc4[3] = false;
@test mc4 == MBitCol([false, true, false, true])
@test mc4_ == MBitCol([false, true, false, true])


@info "Tesing iteration"
@test sc4 == SBitCol((bit for bit in sc4))


@info "Tesing vcat"
@test [sc4; mc4] == SBitCol([0, 0, 1, 1, 0, 1, 0, 1])
@test [sr4 mr4] == SBitRow([1, 0, 0, 1, 1, 0, 1, 0])


@info "Tesing widening arithmetic"
@test -sc4 == [0, 0, -1, -1]
@test 1.5*mc4 == [0, 1.5, 0, 1.5]
@test sc4/2 == [0, 0, 0.5, 0.5]


@info "Testing hamming and parity"
v = SBitCol([1,0,0,1,1,0,1,0,0,1,0,1,1,1,0,1])
@test count(v) == 9
@test sum(v) == 9
@test parity(v) == true
@test hamming(v) == 9
@test hamming(v, MBitCol([1,0,1,1,0,0,1,0,1,1,0,1,1,0,1,0])) == 6

@info "Testing bitmapped operations"
@test ~sc4 == SBitCol([1,1,0,0])
@test ~mr4 == MBitRow([0,1,0,1])
@test sc4 & mc4 == SBitCol([0,0,0,1])
@test sc4 | mc4 == SBitCol([0,1,1,1])
@test sc4 ⊻ mc4 == SBitCol([0,1,1,0])
@test sr4 & mr4 == SBitRow([1,0,0,0])
@test sr4 | mr4 == SBitRow([1,0,1,1])
@test sr4 ⊻ mr4 == SBitRow([0,0,1,1])