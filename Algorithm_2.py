"""
Copyright © Institute of Mathematical and Computational Sciences (IMACS),
Ho Chi Minh City University of Technology (HCMUT).

Contact: Prof. Phan Thanh An <thanhan@hcmut.edu.vn>

This module implements Algorithm 2 for computing the robustness index
of quasiconvex functions on a closed interval [a, b], as presented in:

    N.N. Hai, P.T. An*, N.H. Hai
    "On the Robustness of Quasiconvex Functions"
    Journal of Computational and Applied Mathematics, 2025.
    DOI: https://doi.org/10.1016/j.cam.2025.116906

*Corresponding author: P.T. An (thanhan@hcmut.edu.vn)

Source code and numerical examples are available at:
    https://github.com/hoanghaihcmut/robust-index

Title: Algorithm 2
"""

from sympy import (
    symbols, solveset, diff, Interval, oo, EmptySet,
    is_decreasing, is_increasing, is_convex
)


def is_quasiconvex(f, x, a, b):
    """
    Check whether function `f` is quasiconvex on the interval [a, b].

    Parameters:
        f : sympy expression
            The function to test.
        x : sympy Symbol
            The variable.
        a, b : float
            Interval endpoints.

    Returns:
        bool
            True if f is quasiconvex on [a, b], False otherwise.
    """
    I = Interval(a, b)
    df = diff(f, x)

    if is_decreasing(f, I) or is_increasing(f, I):
        return True

    critical_points = list(solveset(df, x, domain=I))
    if not critical_points:
        return False

    x_min = min(critical_points, key=lambda t: f.subs(x, t))
    I_left = Interval(I.start, x_min)
    I_right = Interval(x_min, I.end)

    return is_decreasing(f, I_left) and is_increasing(f, I_right)


def algorithm_2(f, x, a, b):
    """
    Compute the robustness index of a quasiconvex function on [a, b]
    using Algorithm 2 from the paper.

    Parameters:
        f : sympy expression
            The quasiconvex function.
        x : sympy Symbol
            The variable.
        a, b : float
            Interval endpoints.

    Returns:
        float
            The robustness index, or -∞ if the function is not robust.
    """
    I = Interval(a, b)
    sf = -oo

    if is_convex(f, x, domain=I):
        return oo

    df = diff(f, x)
    dff = diff(df, x)
    sol_dff_neg = solveset(dff < 0, x, domain=Interval.open(a, b))

    candidates = []
    sol_df_pos = solveset(df > 0, x, domain=sol_dff_neg)
    if sol_df_pos != EmptySet:
        candidates.append(df.subs(x, sol_df_pos.sup))

    sol_df_neg = solveset(-df > 0, x, domain=sol_dff_neg)
    if sol_df_neg != EmptySet:
        candidates.append(-df.subs(x, sol_df_neg.inf))

    if not candidates:
        return -oo

    sf_diamond = min(candidates)
    if sf_diamond == 0:
        return -oo

    sf = sf_diamond
    for t in solveset(dff, x, domain=Interval.open(a, b)):
        df_t = df.subs(x, t)
        if not is_quasiconvex(f - df_t * x, x, a, b):
            sf = min(sf, abs(df_t))

    return sf


if __name__ == "__main__":
    x = symbols('x')
    a = 0
    b = 1
    f = (1/3) * x**3 - 2 * x**2 + 4 * x

    sf = algorithm_2(f, x, a, b)

    print(f"f(x) = {f}, D = [{a}, {b}]")
    print(f"Robustness index (sf) = {sf}")
