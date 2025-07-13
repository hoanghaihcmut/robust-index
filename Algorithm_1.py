"""
Copyright © Institute of Mathematical and Computational Sciences (IMACS),
Ho Chi Minh City University of Technology (HCMUT).

Contact: Prof. Phan Thanh An <thanhan@hcmut.edu.vn>

This module implements Algorithm 1 for computing the robustness index
of quasiconvex functions on a closed interval [a, b], as presented in:

    N.N. Hai, P.T. An*, N.H. Hai
    "On the Robustness of Quasiconvex Functions"
    Journal of Computational and Applied Mathematics, 2025.
    DOI: https://doi.org/10.1016/j.cam.2025.116906

*Corresponding author: P.T. An (thanhan@hcmut.edu.vn)

Source code and numerical examples are available at:
    https://github.com/hoanghaihcmut/robust-index

Title: Algorithm 1
"""

from sympy import (
    symbols, solveset, diff, Interval, oo, is_convex
)


def algorithm_1(f, x, a, b, gamma=1e-2, zmin=1e-323, alpha=0):
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
        gamma : float, optional
            Step size increment (default: 1e-2).
        zmin : float, optional
            Minimum allowable value before declaring failure (default: 1e-323).
        alpha : float, optional
            Initial alpha value (default: 0).

    Returns:
        float
            The robustness index, or -∞ if the function is not robust.
    """
    I = Interval(a, b)

    def L_abs_df(alpha):
        return solveset(abs(diff(f, x)) <= alpha, x, domain=I)

    if is_convex(f, x, domain=I):
        return oo

    flag = False
    alpha += gamma

    while is_convex(f, x, domain=L_abs_df(alpha)):
        flag = True
        alpha += gamma

    if flag:
        return alpha - gamma

    while not is_convex(f, x, domain=L_abs_df(alpha)):
        alpha -= gamma
        if alpha <= zmin:
            return -oo
        if is_convex(f, x, domain=L_abs_df(alpha)):
            return alpha

    return -oo


if __name__ == "__main__":
    x = symbols('x')

    a = 0
    b = 1
    f = (1/3) * x**3 - 2 * x**2 + 2 * x

    sf = algorithm_1(f, x, a, b, gamma=1e-1)

    print(f"f(x) = {f}, D = [{a}, {b}]")
    print(f"Robustness index (sf) = {sf}")
