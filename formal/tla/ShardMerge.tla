----------------------------- MODULE ShardMerge -----------------------------
\* Line-aligned chunking of a byte buffer, and binary reductions of natural
\* counts. Every newline mask of a buffer of length BufLen is a state. Chunks
\* cover the buffer once, none is empty, and a chunk ends on a newline or at
\* the end of the buffer. Separately, every parenthesization of a short vector
\* of natural numbers sums to the sequential total: the sharded count merge is
\* a monoid, so no count is lost or counted twice. Floating-point addition is
\* not that monoid and is not the property checked here.

EXTENDS Integers, Sequences, FiniteSets, TLC

CONSTANT BufLen
ASSUME BufLen \in Nat /\ BufLen >= 1 /\ BufLen <= 8

RECURSIVE Pow2(_)
Pow2(n) == IF n <= 0 THEN 1 ELSE 2 * Pow2(n - 1)

IsNl(bits, p) == (bits \div Pow2(p - 1)) % 2 = 1

RECURSIVE FirstNl(_, _, _)
FirstNl(bits, p, len) ==
  IF p > len THEN 0
  ELSE IF IsNl(bits, p) THEN p
  ELSE FirstNl(bits, p + 1, len)

RECURSIVE Parts(_, _, _)
Parts(bits, start, target) ==
  IF start >= BufLen THEN << >>
  ELSE LET end0 == IF start + target > BufLen THEN BufLen ELSE start + target
           nl == FirstNl(bits, end0 + 1, BufLen)
           ended == IF end0 >= BufLen THEN BufLen ELSE IF nl = 0 THEN BufLen ELSE nl
       IN << <<start + 1, ended>> >> \o Parts(bits, ended, target)

RECURSIVE Flatten(_)
Flatten(chunks) ==
  IF chunks = << >> THEN << >>
  ELSE LET lo == Head(chunks)[1]
           hi == Head(chunks)[2]
           span == [p \in 1..(hi - lo + 1) |-> lo + p - 1]
       IN span \o Flatten(Tail(chunks))

EndsWell(bits, chunks) ==
  \A c \in 1..Len(chunks) :
    LET hi == chunks[c][2]
    IN hi = BufLen \/ IsNl(bits, hi)

CoverOk(bits, target) ==
  LET chunks == Parts(bits, 0, target)
      flat == Flatten(chunks)
  IN /\ \A c \in 1..Len(chunks) : chunks[c][1] <= chunks[c][2]
     /\ flat = [p \in 1..BufLen |-> p]
     /\ EndsWell(bits, chunks)

PartitionOk(bits) == \A target \in 1..BufLen : CoverOk(bits, target)

Left(f, i) == [j \in 1..i |-> f[j]]
Right(f, i, k) == [j \in 1..(k - i) |-> f[i + j]]

RECURSIVE SumFun(_, _)
SumFun(f, k) == IF k = 0 THEN 0 ELSE f[k] + SumFun(f, k - 1)

RECURSIVE TreeSums(_, _)
TreeSums(f, k) ==
  IF k = 0 THEN {0}
  ELSE IF k = 1 THEN {f[1]}
  ELSE UNION { { a + b : a \in TreeSums(Left(f, i), i),
                          b \in TreeSums(Right(f, i, k), k - i) }
               : i \in 1..(k - 1) }

ReductionOk ==
  \A k \in 1..3 :
    \A f \in [1..k -> 0..2] :
      TreeSums(f, k) = {SumFun(f, k)}

VARIABLE i

Init == i = 0

Next == IF i < Pow2(BufLen) - 1 THEN i' = i + 1 ELSE UNCHANGED i

Spec == Init /\ [][Next]_i

Inv == i \in 0..(Pow2(BufLen) - 1) /\ PartitionOk(i) /\ ReductionOk

=============================================================================
