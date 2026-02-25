#!/usr/bin/env python3
"""
Routh–Hurwitz table generator with special cases:
- Zero in the first column -> epsilon substitution
- Entire row of zeros -> auxiliary polynomial method

Input: polynomial coefficients in descending powers:
    a_n s^n + a_{n-1} s^{n-1} + ... + a_0
Example:
    coeffs = [1, 2, 3, 4]  # s^3 + 2 s^2 + 3 s + 4

Output:
- Routh table
- Number of RHP roots inferred from sign changes in first column
"""

from __future__ import annotations
from dataclasses import dataclass
from typing import List, Tuple, Optional
import math


@dataclass
class RouthResult:
    table: List[List[float]]
    first_column: List[float]
    sign_changes: int
    rhp_roots: int
    notes: List[str]


def _sign(x: float, tol: float) -> int:
    if abs(x) <= tol:
        return 0
    return 1 if x > 0 else -1


def _poly_degree_from_coeffs(coeffs: List[float]) -> int:
    # Remove leading zeros safely
    i = 0
    while i < len(coeffs) and abs(coeffs[i]) == 0:
        i += 1
    if i == len(coeffs):
        raise ValueError("All coefficients are zero, polynomial is undefined.")
    return len(coeffs) - 1 - i


def _trim_leading_zeros(coeffs: List[float], tol: float) -> List[float]:
    i = 0
    while i < len(coeffs) and abs(coeffs[i]) <= tol:
        i += 1
    if i == len(coeffs):
        raise ValueError("All coefficients are (near) zero within tolerance.")
    return coeffs[i:]


def _auxiliary_polynomial(prev_row: List[float], power: int) -> List[float]:
    """
    Build auxiliary polynomial A(s) from the row above the zero row.
    If the zero row corresponds to s^power, and prev_row has terms:
        a*s^power + b*s^(power-2) + c*s^(power-4) + ...
    Return coefficients of A(s) in descending powers with gaps removed in steps of 2.
    """
    # prev_row already aligns with descending even/odd powers pattern.
    # Build as a dense list where missing odd powers are zero.
    # Example: power=4, prev_row=[a,b,c] => A(s)=a*s^4 + b*s^2 + c
    deg = power
    dense = [0.0] * (deg + 1)
    p = deg
    for val in prev_row:
        if p < 0:
            break
        dense[deg - p] = val  # place at index corresponding to s^p
        p -= 2
    return dense


def _differentiate_dense_poly(dense: List[float]) -> List[float]:
    """
    Differentiate dense polynomial coefficients (descending powers).
    dense[k] corresponds to s^(n-k).
    """
    n = len(dense) - 1
    deriv = []
    for i, a in enumerate(dense[:-1]):
        power = n - i
        deriv.append(a * power)
    return deriv  # descending powers, degree n-1


def routh_hurwitz(
    coeffs: List[float],
    tol: float = 1e-12,
    epsilon: float = 1e-6,
    use_symbolic_epsilon: bool = False,
) -> RouthResult:
    """
    Compute Routh table for polynomial with coefficients in descending powers.

    Special cases handled:
    1) First element of a row is zero (but row not all zeros):
       substitute epsilon (small positive) and proceed.
       If use_symbolic_epsilon=True, uses +epsilon with provided value but notes it.

    2) Entire row becomes (near) zero:
       use auxiliary polynomial from previous row, differentiate it,
       and replace the zero row with derivative coefficients aligned to the row.

    Returns RouthResult with table, first column, sign changes, and notes.
    """
    notes: List[str] = []
    coeffs = _trim_leading_zeros([float(c) for c in coeffs], tol)
    n = len(coeffs) - 1  # degree

    cols = (n // 2) + 1  # number of columns in Routh table
    rows = n + 1

    # Initialize table with zeros
    rt = [[0.0 for _ in range(cols)] for _ in range(rows)]

    # Fill first two rows
    # Row 0: coefficients of s^n, s^(n-2), s^(n-4), ...
    # Row 1: coefficients of s^(n-1), s^(n-3), ...
    rt[0] = [coeffs[i] if i < len(coeffs) else 0.0 for i in range(0, len(coeffs), 2)]
    rt[1] = [coeffs[i] if i < len(coeffs) else 0.0 for i in range(1, len(coeffs), 2)]
    rt[0] += [0.0] * (cols - len(rt[0]))
    rt[1] += [0.0] * (cols - len(rt[1]))

    # Helper for row "power" mapping:
    # Row index r corresponds to s^(n-r)
    def row_power(r: int) -> int:
        return n - r

    # Build subsequent rows
    for r in range(2, rows):
        # Special case: if the previous row is all zeros, apply auxiliary method.
        if all(abs(x) <= tol for x in rt[r - 1]):
            power = row_power(r - 1)
            prev = rt[r - 2]
            notes.append(
                f"Row {r} (s^{power}) is zero row, using auxiliary polynomial from row {r-1} (s^{power+1})."
            )
            aux_dense = _auxiliary_polynomial(prev, power + 1)
            deriv_dense = _differentiate_dense_poly(aux_dense)

            # Now map deriv_dense back to a Routh row aligned with powers power, power-2, ...
            # deriv_dense is degree=power, descending all powers.
            # Take coefficients at powers: power, power-2, power-4, ...
            deriv_row = []
            deg = len(deriv_dense) - 1
            # deriv_dense index i corresponds to s^(deg-i)
            for p in range(power, -1, -2):
                i = deg - p
                if 0 <= i < len(deriv_dense):
                    deriv_row.append(deriv_dense[i])
                else:
                    deriv_row.append(0.0)
            deriv_row += [0.0] * (cols - len(deriv_row))
            rt[r - 1] = deriv_row

        # If first element of the row above is zero, epsilon substitution
        if abs(rt[r - 1][0]) <= tol and not all(abs(x) <= tol for x in rt[r - 1]):
            notes.append(
                f"Row {r} (s^{row_power(r-1)}): first element is zero, substituting epsilon={epsilon}."
            )
            rt[r - 1][0] = epsilon if not use_symbolic_epsilon else epsilon

        # Compute row r
        for c in range(cols - 1):
            a = rt[r - 2][0]
            b = rt[r - 2][c + 1]
            d = rt[r - 1][0]
            e = rt[r - 1][c + 1]
            # Formula: (d*b - a*e) / d
            # where d = first element of previous row (rt[r-1][0])
            if abs(d) <= tol:
                # If still zero, force epsilon and continue
                notes.append(
                    f"Row {r} computation: pivot is zero at row {r-1}, forcing epsilon={epsilon}."
                )
                d = epsilon
                rt[r - 1][0] = d
            rt[r][c] = (d * b - a * e) / d

        # Last column remains 0 by construction

        # Clean tiny values
        for c in range(cols):
            if abs(rt[r][c]) <= tol:
                rt[r][c] = 0.0

    # First column and sign changes
    first_col = [rt[r][0] for r in range(rows)]

    # Determine sign changes ignoring zeros by treating zero as limit case:
    # Common practical approach: replace zeros with tiny epsilon keeping last nonzero sign.
    signs: List[int] = []
    last_nonzero_sign: Optional[int] = None
    for x in first_col:
        s = _sign(x, tol)
        if s == 0:
            # If leading zeros, take + sign conventionally, otherwise keep last sign
            if last_nonzero_sign is None:
                s = 1
            else:
                s = last_nonzero_sign
        else:
            last_nonzero_sign = s
        signs.append(s)

    sign_changes = 0
    for i in range(1, len(signs)):
        if signs[i] != signs[i - 1]:
            sign_changes += 1

    return RouthResult(
        table=rt,
        first_column=first_col,
        sign_changes=sign_changes,
        rhp_roots=sign_changes,
        notes=notes,
    )


def format_routh_table(table: List[List[float]], degree: int, width: int = 12, prec: int = 6) -> str:
    cols = len(table[0])
    lines = []
    for r, row in enumerate(table):
        power = degree - r
        label = f"s^{power:<2d}"
        fmt_row = []
        for v in row:
            if v == 0.0:
                fmt_row.append(f"{0:>{width}d}")
            else:
                fmt_row.append(f"{v:>{width}.{prec}g}")
        lines.append(f"{label} | " + " ".join(fmt_row[:cols]))
    return "\n".join(lines)


def main():
    # Example usage
    coeffs = [1, 2, 3, 4]  # s^3 + 2 s^2 + 3 s + 4
    res = routh_hurwitz(coeffs)

    degree = len(_trim_leading_zeros(coeffs, 1e-12)) - 1
    print("Routh–Hurwitz table:")
    print(format_routh_table(res.table, degree))
    print()
    print("First column:", res.first_column)
    print("Sign changes:", res.sign_changes)
    print("RHP roots (unstable poles):", res.rhp_roots)
    if res.notes:
        print("\nNotes:")
        for n in res.notes:
            print("-", n)


if __name__ == "__main__":
    main()
