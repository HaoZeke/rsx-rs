"""Exactness checks for slice bounds, CRT reconstruction, and the Fisher step."""

from sympy import Integer, binomial, factorial, gcd
from sympy.ntheory.modular import crt


def test_limb_mul_add_fits_in_128_bits():
    bound = Integer(2) ** 64
    top = bound - 1
    wide = top * top + top
    assert wide < bound**2
    assert wide == top * bound
    assert wide // bound == top
    assert wide % bound == 0


def test_low_high_split_is_an_identity():
    samples = [0, 1, 2, 17, 2**53 - 1, 2**53 + 1, 2**64 - 1, 2**80 + 12345]
    for shift in (1, 8, 32, 53, 64):
        base = Integer(2) ** shift
        for value in samples:
            value = Integer(value)
            assert (value // base) * base + (value % base) == value
            assert value % base < base


def test_nonempty_slice_product_respects_the_bit_budget():
    width = 7
    slices = [2**width - 1, 0, 3, 19, 100]
    product = Integer(1)
    for piece in slices:
        assert piece < 2**width
        product *= piece
    assert product < Integer(2) ** (len(slices) * width)
    assert len(slices) * width <= 64
    assert product < Integer(2) ** 64


def test_scaled_quotient_lands_in_the_rounding_window():
    pairs = [(1, 1), (1, 2), (3, 2), (1, 3), (2**53 + 1, 1), (2**53 + 3, 1), (10**18, 10**6), (7, 2**100)]
    for numerator, denominator in pairs:
        numerator = Integer(numerator)
        denominator = Integer(denominator)
        diff = int(numerator).bit_length() - int(denominator).bit_length()
        if diff <= -1076 or diff > 1024:
            continue
        scale = 54 + int(denominator).bit_length() - int(numerator).bit_length()
        scaled = numerator * Integer(2) ** scale / denominator
        assert Integer(2) ** 53 < scaled < Integer(2) ** 55


def test_crt_reconstruction_is_unique_below_the_product():
    moduli = [256, 255, 253, 251, 247, 239, 233, 229]
    for index, left in enumerate(moduli):
        for right in moduli[:index]:
            assert gcd(left, right) == 1
    value = Integer(10) ** 12 + 123456
    residues = [int(value % modulus) for modulus in moduli]
    reconstructed, product = crt(moduli, residues)
    assert reconstructed == value % product
    assert 0 <= reconstructed < product


def test_fisher_recurrence_divides_exactly():
    cases = [(5, 4, 9, 1), (12, 7, 20, 3), (8, 8, 16, 0), (30, 20, 40, 10)]
    for row1, col1, total, cell in cases:
        absent = total - row1
        assert 0 <= cell < row1 and cell < col1
        assert col1 - cell <= absent
        left = binomial(row1, cell)
        right = binomial(absent, col1 - cell)
        step = left * right * (row1 - cell) * (col1 - cell)
        denominator = (cell + 1) * (absent - (col1 - cell) + 1)
        assert step % denominator == 0
        nxt = binomial(row1, cell + 1) * binomial(absent, col1 - (cell + 1))
        assert step // denominator == nxt


def test_uniform_bayes_factorial_and_binomial_forms_agree():
    tables = [(1, 1, 0, 1), (3, 8, 1, 5), (10, 10, 10, 10), (4, 9, 0, 7)]
    for k1, n1, k2, n2 in tables:
        total = n1 + n2
        present = k1 + k2
        factorial_form = (
            factorial(k1)
            * factorial(n1 - k1)
            / factorial(n1 + 1)
            * factorial(k2)
            * factorial(n2 - k2)
            / factorial(n2 + 1)
            * factorial(total + 1)
            / (factorial(present) * factorial(total - present))
        )
        binomial_form = (
            Integer(total + 1)
            / ((n1 + 1) * (n2 + 1))
            * binomial(total, present)
            / (binomial(n1, k1) * binomial(n2, k2))
        )
        assert factorial_form == binomial_form


def test_integer_prior_reduces_to_twenty_one_over_twenty():
    # Beta(2, 3) on both groups and the null, counts (1 of 2) and (0 of 1).
    def rising(start, length):
        value = Integer(1)
        for offset in range(length):
            value *= start + offset
        return value

    def marginal(k, n, alpha, beta):
        return rising(alpha, k) * rising(beta, n - k) / rising(alpha + beta, n)

    ratio = marginal(1, 2, 2, 3) * marginal(0, 1, 2, 3) / marginal(1, 3, 2, 3)
    assert ratio == Integer(21) / 20
