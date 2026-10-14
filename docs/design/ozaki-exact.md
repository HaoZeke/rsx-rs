# Exact slice arithmetic for association tests

The only matrix product in this analysis is the Gram matrix of a marker-by-individual depth table, and that product is already an integer. Re-running an Ozaki scheme on it does not help. The same splitting idea does apply to the association tests, where a probability is a ratio of integers that the current kernels evaluate with `lgamma`. This note records what that generalisation is, what it costs on a CPU, and which pieces are machine-checked.

## The Gram is already an integer

A principal-component analysis forms \(G = X^{\mathsf T} X\) from marker depths. A typical contrast has a few dozen individuals; the shapes measured here are 40, 64, and 128. Depths are non-negative integers that fit in 16 bits, so each rank-1 update is an integer outer product. Accumulating those products in 64-bit registers is exact. On the tables below, the resulting entries also happened to be exact binary64 values, which is why a SIMD binary64 kernel agreed with the integer product bit for bit.

Ozaki scheme II rebuilds that same integer from several residue products and the Chinese Remainder Theorem (Ozaki, Uchino, and Imamura 2026; preprint arXiv:2504.08009). It is the algorithm that makes an integer matrix product run on INT8 tensor cores. On this 8-core Xeon, in a release build, it does not beat the direct product:

| markers | individuals | binary64 SIMD | scalar `u64` | scheme II | scheme II / `u64` |
| --- | ---: | ---: | ---: | ---: | ---: |
| 100000 | 40 | 17.5 ms | 48.1 ms | 138.7 ms | 2.9 |
| 20000 | 128 | 28.5 ms | 103.4 ms | 134.2 ms | 1.3 |
| 4096 | 64 | 1.54 ms | 4.92 ms | 7.90 ms | 1.6 |

Every one of those products was exact in binary64. A binary32 accumulator was exact on the two smaller shapes and was off by as much as \(4.7 \times 10^{4}\) on the \(10^{5} \times 40\) table, which is the expected overflow past \(2^{24}\). The paper's tensor-core figures, 7.4–9.8 TFLOPS on an RTX 4090 and 56.6–80.2 TFLOPS on a GH200, are the authors' measurements, not ours. This machine has no GPU. The device Gram still adds binary64 products with an atomic add, and the device Fisher and Bayes factor still call `lgamma`.

Scheme I (Ozaki, Ogita, Oishi, and Rump 2012) splits a floating-point matrix so that each slice product is exact in working precision. The later integer-unit papers (Ootomo, Ozaki, and Yokota 2024; Uchino, Ozaki, and Imamura 2025) speed that split up on INT8 units. Those are different algorithms from the CRT scheme, and neither one is a faster CPU Gram for a 40-column depth matrix.

## Where a ratio of integers is being approximated

Fisher's exact test for a \(2 \times 2\) table (Fisher 1922, 1935) sums hypergeometric probabilities of every table at least as unlikely as the observed one. The implementation evaluates each probability as a sum of `lgamma` calls and keeps a table when its log probability sits within \(10^{-9}\) of the observed value, then clamps the result to \([10^{-16}, 1]\).

The Bayes factor with a Beta prior is a ratio of beta functions, again through `lgamma`. For the uniform prior \(\mathrm{Beta}(1,1)\) the marginal does not include the binomial coefficient of the likelihood, so the factor still depends on the counts. The Yates chi-squared statistic (Yates 1934) is a short floating-point expression, and its p-value is an `erfc`. Benjamini–Hochberg adjustment (Benjamini and Hochberg 1995) is a rank scaling of whatever p-values it is given. Sharded table merges add integer counts.

The \(10^{-9}\) band and the `lgamma` evaluation are the accuracy question. The Gram is not.

## Slices and residues, without a matrix

Two facts from that literature are enough, and neither one mentions a matrix.

A 64-bit limb, times a 64-bit factor, plus a 64-bit carry, is strictly less than \(2^{128}\). Indeed
\[
(2^{64}-1)^{2} + (2^{64}-1) = (2^{64}-1)\,2^{64} = 2^{128} - 2^{64}.
\]
The low and high halves of that sum are therefore an exact split: no rounding and no lost carry. The same counting bound says a nonempty product of \(s\)-bit slices fits in \(\mathrm{length} \times s\) bits, and fits in a \(w\)-bit word whenever \(\mathrm{length} \times s \le w\). This is the integer form of the error-free product split of Dekker (1971) and of the accurate dot-product theory of Ogita, Rump, and Oishi (2005). Once the sum is an integer, the compensated floating-point summation of Demmel and Hida (2004) is unnecessary for that sum.

Scheme II's reconstruction is the Chinese Remainder Theorem. If the moduli are pairwise coprime and two integers both lie strictly below their product and agree modulo every modulus, the two integers are equal. A reconstruction that meets the congruences and the range bound is the integer it claims to be.

Both statements are proved in Lean 4 with Mathlib, from the binomial recurrence and the existing CRT lemmas, with no unfinished proof. A SymPy suite checks the same bounds on concrete integers, checks that the scaled quotient used for rounding falls in \((2^{53}, 2^{55})\), rebuilds a product from the scheme II moduli, and checks the Fisher step and the two closed forms of the uniform Bayes factor.

## Fisher's probability as one rounding

The weight of a table with fixed margins is
\[
w(a) = \binom{\mathrm{row}_1}{a}\binom{N-\mathrm{row}_1}{\mathrm{col}_1-a}.
\]
The sum of the weights is \(\binom{N}{\mathrm{col}_1}\). Stepping the free cell by one multiplies by four positive integers and divides exactly:
\[
w(a+1) = w(a)\,\frac{(\mathrm{row}_1-a)(\mathrm{col}_1-a)}{(a+1)\,(d+1)},
\]
where \(d\) is the opposite cell. The identity is an equality in the natural numbers, so the division in the implementation is exact. The two-sided probability is the sum of the weights that do not exceed the observed weight, divided by the sum of all of them. That division is performed once.

The positive rational is rounded to binary64 in round-to-nearest, ties to even. Write \(b(\cdot)\) for the bit length. With
\[
k = 54 + b(\mathrm{den}) - b(\mathrm{num}),
\]
the scaled quotient \(q = \lfloor \mathrm{num}\cdot 2^{k} / \mathrm{den}\rfloor\) lies in \([2^{53}, 2^{55})\). One extra bit is the round bit. Anything below it, including a nonzero remainder, is sticky. A subnormal result shifts relative to the least subnormal unit \(2^{-1074}\). The tie between 0 and that unit, \(2^{-1075}\), rounds to 0. A ratio whose bit lengths already place it strictly below that tie returns 0, and one strictly above \(2^{1024}\) returns infinity.

Invalid margins return 1. The value reported as a p-value is then clamped to \([10^{-16}, 1]\), which is the range the association commands already use. The unclamped probability is available separately.

On every \(2 \times 2\) table whose two group sizes run from 1 through 40 (23 820 tables), the clamped integer p-value and the `lgamma` p-value differ in their bits for 22 318 tables, and the worst relative gap is \(1.5 \times 10^{-13}\). At these sizes the \(10^{-9}\) admission band does not move a probability by a visible amount. The integer value is the correctly rounded rational; the `lgamma` value is already very close to it.

It is also slower. Both groups of size \(n\), with \(n/2\) and \(n/5\) individuals carrying the marker, release build, this host:

| \(n\) | `lgamma` Fisher | integer Fisher | ratio |
| ---: | ---: | ---: | ---: |
| 10 | 0.49 µs | 0.63 µs | 1.3 |
| 20 | 0.73 µs | 1.38 µs | 1.9 |
| 40 | 1.28 µs | 6.35 µs | 5.0 |
| 80 | 2.40 µs | 13.4 µs | 5.6 |
| 200 | 5.69 µs | 46.2 µs | 8.1 |

The integer kernel is the one to call when the p-value has to be that rounded rational. It is not a faster Fisher test, and the command path still calls `lgamma`, so the CPU and the device keep a single definition.

## Integer priors for the Bayes factor

For positive integers \(\alpha, \beta\), the beta-function ratio behind the Bayes factor is a ratio of rising factorials. The product \(\alpha(\alpha+1)\cdots(k+\alpha-1)\) is empty, and equals 1, when \(k = 0\). For the uniform prior this collapses to
\[
\frac{N+1}{(n_1+1)(n_2+1)}\cdot\frac{\binom{N}{K}}{\binom{n_1}{k_1}\binom{n_2}{k_2}},
\]
with \(K = k_1+k_2\) and \(N = n_1+n_2\), which is the same rational as the factorial form the `lgamma` kernel computes. A direct reduction of \(\mathrm{Beta}(2,3)\) on the counts “1 of 2” and “0 of 1” is \(21/20\).

Priors that are not integers in \(1,\ldots,4096\) return no value. Every such integer is an exact binary64 value. A Jeffreys prior \(\mathrm{Beta}(\tfrac12,\tfrac12)\) is outside this set and stays on `lgamma`.

The uniform integer kernel is slower at every size measured, and the gap opens once the binomial coefficient no longer fits in 128 bits:

| \(n\) | `lgamma` Bayes factor | integer Bayes factor | ratio |
| ---: | ---: | ---: | ---: |
| 10 | 0.11 µs | 0.26 µs | 2.4 |
| 20 | 0.11 µs | 0.55 µs | 5.0 |
| 40 | 0.10 µs | 1.03 µs | 10 |
| 80 | 0.10 µs | 32 µs | 320 |
| 200 | 0.11 µs | 171 µs | 1600 |

`lgamma` does not walk the factorial, so its time barely moves. The integer path walks it. Exactness is inexpensive only while the coefficient fits in one machine word; for these group sizes it does not. The log of the same rational is available for callers who want a log Bayes factor without subtracting two logarithms.

## Chi-squared, FDR, and sharded counts

The Yates statistic is a rational: \(N y^{2}\) over the four margins when \(N\) is even, and \(N(2\Delta-N)^{2}/(4\cdot\mathrm{margins})\) when \(N\) is odd, with \(\Delta = |ad-bc|\) and \(y = \Delta - N/2\). For every table with margins below 12, that rational rounds to the same binary64 value as the current expression. The unbalanced table with present counts \(9\,045\,414\) and \(2\,988\,579\,416\) and margins \(3\,273\,968\,006\) and \(3\,294\,916\,954\) differs by one unit in the last place. A much larger table whose cross-product is almost balanced, present counts \(3 \times 10^{9}\) and \(1\) out of two groups of \(4 \times 10^{9}\), rounds to the same value either way. The rational evaluation takes 1.5–1.7 µs. The floating-point expression takes about 2–5 ns. The rational form should not replace it. The p-value is still an `erfc` of the statistic.

Benjamini–Hochberg adjustment multiplies a p-value by an integer rank ratio and takes a running minimum. Its residual error is the error of the p-value.

Count shards are vectors of natural numbers. Addition of values that fit in 64 bits is a commutative monoid, and a coordinate sum that does not fit is rejected rather than wrapped. Every parenthesization therefore agrees with the left fold: a parallel reduction does not drop a count and does not count one twice. The chunks that feed that reduction are cut on newlines. A chunk ends on a newline or at the end of the buffer, chunks are non-empty, and their concatenation is the whole buffer. TLC checked every newline mask of a buffer of length 5, and every parenthesization of vectors of length at most 3 with entries in \(\{0,1,2\}\): 32 distinct states, no counterexample. Floating-point addition is not a monoid, and the specification does not claim that it is.

## A certified logarithm on the significand

The log of a positive rational is \(\ln m + e \ln 2\), with \(m \in [1, 2)\) taken from the same correctly rounded significand. Thirty-two panels of width \(1/32\) cover \([1, 2)\). Each panel is a degree-8 minimax polynomial in the centred variable \(u = m - c\), produced by Sollya's `fpminimax`. The absolute error of each polynomial is bounded by Sollya's `supnorm`, whose algorithm is that of Chevillard, Harrison, Joldes, and Lauter (2011). The worst upper endpoint is
\[
2^{-55} \times \mathtt{0x1.f48729117040bf06373d6p+0} \approx 5.43 \times 10^{-17}.
\]
The constant stored with the coefficients is \(6 \times 10^{-17}\), which sits strictly above that endpoint and is less than four times it. A Horner evaluation on 8192 points of \([1, 2)\) differed from the platform logarithm by at most \(1.11 \times 10^{-16}\). That figure includes evaluation error. The polynomial certificate does not cover it, and \(e \ln 2\) is a single binary64 product. The result is not a correctly rounded logarithm.

Correctly rounded \(\ln\) is available from CRlibm (de Dinechin, Lauter, and Muller 2007; the library is Daramy-Loirat et al. 2006) and from CORE-MATH (Sibidanov, Zimmermann, and Glondu 2022). A correctly rounded binary64 `lgamma` has since been added to CORE-MATH and to glibc 2.43. Integer arguments of the association tests do not need it: the weight stays an integer until the one division at the end.

## What was checked

Lean 4.35 and Mathlib prove the 128-bit limb bound, the split identity, the slice-product bit bound, uniqueness of a CRT reconstruction for two moduli and for a pairwise-coprime list, and the exact Fisher step. The build finishes with no unfinished proof.

SymPy checks those bounds, the rounding window, CRT reconstruction on the scheme II moduli, exact division in the Fisher step, agreement of the factorial and binomial forms of the uniform Bayes factor, and the fraction \(21/20\).

TLC checks the newline partition and the natural-number reduction on the finite model above. The recorded run reports model checking completed with no error, 32 distinct states.

Sollya 8.0 emits the 288 panel coefficients and the absolute-error endpoint. A checker refuses the build if those coefficients differ from the ones compiled into the logarithm, or if \(6 \times 10^{-17}\) does not lie just above the certified endpoint.

## What remains open

There is no tensor-core timing in this work. The device kernels are unchanged. Non-integer priors stay on `lgamma`. The correctly rounded division is tested, including the ties at \(2^{53}+1\) and at \(2^{-1075}\) and the minimum subnormal, and it is not proved in Lean. The logarithm's certificate is a polynomial bound plus a measured Horner gap. The TLC model is one finite instance. None of the integer kernels replaces the command-path implementation.

## References

Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society B*, 57(1), 289–300.

Chevillard, S., Harrison, J., Joldeş, M., and Lauter, C. (2011). Efficient and accurate computation of upper bounds of approximation errors. *Theoretical Computer Science*, 412(16), 1523–1543. https://doi.org/10.1016/j.tcs.2010.11.052

Daramy-Loirat, C., Defour, D., de Dinechin, F., Gallet, M., Gast, N., Lauter, C., and Muller, J.-M. (2006). *CR-LIBM: A library of correctly rounded elementary functions in double-precision*. Technical report, École Normale Supérieure de Lyon.

de Dinechin, F., Lauter, C., and Muller, J.-M. (2007). Fast and correctly rounded logarithms in double-precision. *RAIRO — Theoretical Informatics and Applications*, 41, 85–102. https://doi.org/10.1051/ita:2007003

Dekker, T. J. (1971). A floating-point technique for extending the available precision. *Numerische Mathematik*, 18, 224–242.

Demmel, J., and Hida, Y. (2004). Accurate and efficient floating point summation. *SIAM Journal on Scientific Computing*, 25(4), 1214–1248. https://doi.org/10.1137/S1064827502407627

Fisher, R. A. (1922). On the interpretation of \(\chi^{2}\) from contingency tables, and the calculation of P. *Journal of the Royal Statistical Society*, 85(1), 87–94.

Fisher, R. A. (1935). The logic of inductive inference. *Journal of the Royal Statistical Society*, 98(1), 39–82.

Ogita, T., Rump, S. M., and Oishi, S. (2005). Accurate sum and dot product. *SIAM Journal on Scientific Computing*, 26(6), 1955–1988. https://doi.org/10.1137/030601818

Ootomo, H., Ozaki, K., and Yokota, R. (2024). DGEMM on integer matrix multiplication unit. *The International Journal of High Performance Computing Applications*, 38(4), 297–313. https://doi.org/10.1177/10943420241239588

Ozaki, K., Ogita, T., Oishi, S., and Rump, S. M. (2012). Error-free transformations of matrix multiplication by using fast routines of matrix multiplication and its applications. *Numerical Algorithms*, 59(1), 95–118. https://doi.org/10.1007/s11075-011-9478-1

Ozaki, K., Uchino, Y., and Imamura, T. (2026). Ozaki scheme II: a GEMM-oriented emulation of floating-point matrix multiplication using an integer modular technique. *The International Journal of High Performance Computing Applications*. https://doi.org/10.1177/10943420261467787

Sibidanov, A., Zimmermann, P., and Glondu, S. (2022). The CORE-MATH project. In *29th IEEE Symposium on Computer Arithmetic (ARITH)*.

Uchino, Y., Ozaki, K., and Imamura, T. (2025). Performance enhancement of the Ozaki Scheme on integer matrix multiplication unit. *The International Journal of High Performance Computing Applications*, 39(3), 462–476. https://doi.org/10.1177/10943420241313064

Yates, F. (1934). Contingency tables involving small numbers and the \(\chi^{2}\) test. *Supplement to the Journal of the Royal Statistical Society*, 1(2), 217–235.
