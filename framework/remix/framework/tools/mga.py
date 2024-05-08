# https://stackoverflow.com/questions/57123194/how-to-distribute-points-evenly-on-the-surface-of-hyperspheres-in-higher-dimensi

from itertools import count
from math import cos
from math import gamma
from math import pi
from math import sin
from math import sqrt
from typing import Callable
from typing import Iterator
from typing import List


def int_sin_m(x: float, m: int) -> float:
    """Computes the integral of sin^m(t) dt from 0 to x recursively,

    Parameters
    ----------
    x : float
        Upper boundary for integration.
    m : int
        Exponent of sinus function.

    Returns
    -------
    float
        Integral value,

    Example
    -------
    >>> from math import pi
    >>> from remix.framework.tools.mga import int_sin_m
    >>> twopi = 2 * pi
    >>> twopi == int_sin_m(2 * pi, 0)
    True
    >>> int_sin_m(2 * pi, 1)
    0.0
    """
    if m == 0:
        return x
    elif m == 1:
        return 1 - cos(x)
    else:
        return (m - 1) / m * int_sin_m(x, m - 2) - cos(x) * sin(x) ** (m - 1) / m


def primes() -> Iterator[int]:
    """Returns an infinite generator of prime numbers."""
    yield from (2, 3, 5, 7)
    composites = {}
    ps = primes()
    next(ps)
    p = next(ps)
    assert p == 3
    psq = p * p
    for i in count(9, 2):
        if i in composites:  # composite
            step = composites.pop(i)
        elif i < psq:  # prime
            yield i
            continue
        else:  # composite, = p*p
            assert i == psq
            step = 2 * p
            p = next(ps)
            psq = p * p
        i += step
        while i in composites:
            i += step
        composites[i] = step


def inverse_increasing(
    func: Callable[[float], float],
    target: float,
    lower: float,
    upper: float,
    atol: float = 1e-10,
) -> float:
    """Returns function inverse value for a target value between lower and upper boundaries.

    Inverse is accurate to an absolute tolerance of atol, and must be monotonically increasing over the interval lower
    to upper,

    Parameters
    ----------
    func : Callable
        Function to be inversed,
    target : float
        Target value for the function evaluation.
    lower : float
        Lower boundary for the value.
    upper : float
        Upper boundary for the value.
    atol : float, optional
        Absolute value tolerance, by default 1e-10,

    Returns
    -------
    float
        Value so that `func(value)` is equal to the `target` (within tolerance `atol`).
    """
    mid = (lower + upper) / 2
    approx = func(mid)
    while abs(approx - target) > atol:
        if approx > target:
            upper = mid
        else:
            lower = mid
        mid = (upper + lower) / 2
        approx = func(mid)
    return mid


def uniform_hypersphere(d: int, n: int) -> List[List[float]]:
    """Generate n points over the d dimensional hypersphere,

    Parameters
    ----------
    d : int
        Dimension of the hypershpere,
    n : int
        Number of points,

    Returns
    -------
    List[List[float]]
        Generated points,
    """
    assert d > 1, f"d must be larger than 1, but is {d}"
    assert n > 0, f"n must be larger than 0, but is {n}"
    points = [[1 for _ in range(d)] for _ in range(n)]
    for i in range(n):
        t = 2 * pi * i / n
        points[i][0] *= sin(t)
        points[i][1] *= cos(t)
    for dim, prime in zip(range(2, d), primes()):
        offset = sqrt(prime)
        mult = gamma(dim / 2 + 0.5) / gamma(dim / 2) / sqrt(pi)

        def dim_func(y):
            return mult * int_sin_m(y, dim - 1)

        for i in range(n):
            deg = inverse_increasing(dim_func, i * offset % 1, 0, pi)
            for j in range(dim):
                points[i][j] *= sin(deg)
            points[i][dim] *= cos(deg)
    return points
