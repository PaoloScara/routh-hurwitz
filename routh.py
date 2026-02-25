#!/usr/bin/env python3
"""
Routh–Hurwitz Stability Analyzer (symbolic + numeric)

Builds the Routh table, counts sign changes, determines stability.
When coefficients contain a parameter K and no numeric value is given,
the table is computed symbolically and stability conditions on K are derived.

Handles special cases automatically:
  - Zero in first column  → ε substitution (numeric) or noted (symbolic)
  - Entire row of zeros   → auxiliary polynomial derivative method

Usage:
    from routh import routh_hurwitz

    routh_hurwitz([1, -4, 1, 6])                        # numeric
    routh_hurwitz([1, 6, 11, 6, "K+2"], K=5)            # parametric, evaluate at K=5
    routh_hurwitz([1, 6, 11, 6, "K+2"])                  # symbolic — derive conditions on K
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Optional, Union
from fractions import Fraction

import sympy as sp

K = sp.Symbol("K", real=True)
_eps = sp.Symbol("epsilon", positive=True)


# ═══════════════════════════════════════════════════════════════════════════
#  Result
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class RouthResult:
    """Routh–Hurwitz analysis result."""
    table: list                     # List[List[sympy.Expr | Fraction]]
    degree: int
    first_column: list
    # Numeric mode
    sign_changes: Optional[int] = None
    rhp_roots: Optional[int] = None
    is_stable: Optional[bool] = None
    # Symbolic mode
    stability_conditions: Optional[list] = None
    stable_range: Optional[str] = None
    notes: List[str] = field(default_factory=list)
    symbolic: bool = False


# ═══════════════════════════════════════════════════════════════════════════
#  Numeric helpers (Fraction-based, exact)
# ═══════════════════════════════════════════════════════════════════════════

def _fmt_frac(v: Fraction) -> str:
    if v == 0:
        return "0"
    if v.denominator == 1:
        return str(v.numerator)
    return f"{v.numerator}/{v.denominator}"


def _sign_frac(x: Fraction) -> int:
    if x > 0: return 1
    if x < 0: return -1
    return 0


def _is_zero_row_frac(row: List[Fraction]) -> bool:
    return all(v == 0 for v in row)


def _aux_deriv_frac(prev_row: List[Fraction], power: int) -> List[Fraction]:
    result, p = [], power
    for val in prev_row:
        if p < 0: break
        result.append(val * p)
        p -= 2
    return result


# ═══════════════════════════════════════════════════════════════════════════
#  Symbolic helpers (sympy-based)
# ═══════════════════════════════════════════════════════════════════════════

def _is_zero_sym(expr) -> bool:
    return sp.simplify(expr) == 0


def _is_zero_row_sym(row) -> bool:
    return all(_is_zero_sym(v) for v in row)


def _aux_deriv_sym(prev_row, power: int):
    result, p = [], power
    for val in prev_row:
        if p < 0: break
        result.append(sp.simplify(val * p))
        p -= 2
    return result


# ═══════════════════════════════════════════════════════════════════════════
#  Core: routh_hurwitz
# ═══════════════════════════════════════════════════════════════════════════

def routh_hurwitz(
    coeffs: List[Union[int, float, str]],
    K_val: Optional[Union[int, float]] = None,
    show: bool = True,
) -> RouthResult:
    """
    Compute the Routh table for a polynomial.

    Parameters
    ----------
    coeffs : list
        Coefficients in descending powers of s.
        Numbers or strings with K (e.g. "K+2", "-10*K", "25+10*K").
    K_val : number or None
        If given, substitute K and compute numerically.
        If None and K appears, compute symbolically → stability conditions.
    show : bool
        Print the formatted table and verdict (default True).

    Returns
    -------
    RouthResult
    """
    # ── Parse coefficients ──────────────────────────────────────────────
    parsed = []
    has_K = False
    for c in coeffs:
        if isinstance(c, str):
            has_K = True
            expr = sp.sympify(c, locals={"K": K})
            parsed.append(expr)
        else:
            parsed.append(sp.Rational(c))

    # If K_val provided, substitute now → numeric mode
    if has_K and K_val is not None:
        parsed = [sp.Rational(v.subs(K, K_val)) if hasattr(v, "subs") else sp.Rational(v) for v in parsed]
        has_K = False

    symbolic = has_K

    if symbolic:
        return _routh_symbolic(parsed, show=show)
    else:
        return _routh_numeric(parsed, show=show)


# ═══════════════════════════════════════════════════════════════════════════
#  Numeric mode (Fraction-based, exact)
# ═══════════════════════════════════════════════════════════════════════════

def _routh_numeric(parsed, show: bool) -> RouthResult:
    resolved = [Fraction(int(v.p), int(v.q)) if hasattr(v, "p") else Fraction(v) for v in parsed]

    # Trim leading zeros
    while resolved and resolved[0] == 0:
        resolved.pop(0)
    if not resolved:
        raise ValueError("All coefficients are zero.")

    n = len(resolved) - 1
    cols = (n // 2) + 1
    rows = n + 1
    notes: List[str] = []
    epsilon = Fraction(1, 1000000)

    rt: List[List[Fraction]] = [[Fraction(0)] * cols for _ in range(rows)]

    for i in range(0, len(resolved), 2):
        if i // 2 < cols: rt[0][i // 2] = resolved[i]
    for i in range(1, len(resolved), 2):
        if i // 2 < cols: rt[1][i // 2] = resolved[i]

    row_power = lambda r: n - r

    for r in range(2, rows):
        if _is_zero_row_frac(rt[r - 1]):
            pw = row_power(r - 1)
            deriv = _aux_deriv_frac(rt[r - 2], pw + 1)
            notes.append(f"s^{pw}: zero row → auxiliary polynomial derivative")
            for j in range(min(len(deriv), cols)):
                rt[r - 1][j] = deriv[j]

        if rt[r - 1][0] == 0 and not _is_zero_row_frac(rt[r - 1]):
            notes.append(f"s^{row_power(r-1)}: zero pivot → ε substitution")
            rt[r - 1][0] = epsilon

        pivot = rt[r - 1][0]
        if pivot == 0:
            pivot = epsilon
            rt[r - 1][0] = pivot

        for c in range(cols - 1):
            rt[r][c] = (pivot * rt[r - 2][c + 1] - rt[r - 2][0] * rt[r - 1][c + 1]) / pivot

    first_col = [rt[r][0] for r in range(rows)]

    last_sign = None
    signs = []
    for x in first_col:
        s = _sign_frac(x)
        if s == 0:
            s = last_sign if last_sign is not None else 1
        else:
            last_sign = s
        signs.append(s)

    sc = sum(1 for i in range(1, len(signs)) if signs[i] != signs[i - 1])

    result = RouthResult(
        table=rt, degree=n, first_column=first_col,
        sign_changes=sc, rhp_roots=sc, is_stable=(sc == 0),
        notes=notes, symbolic=False,
    )
    if show:
        _print_numeric(result)
    return result


# ═══════════════════════════════════════════════════════════════════════════
#  Symbolic mode (sympy-based)
# ═══════════════════════════════════════════════════════════════════════════

def _routh_symbolic(parsed, show: bool) -> RouthResult:
    # Trim leading zeros
    while parsed and _is_zero_sym(parsed[0]):
        parsed.pop(0)
    if not parsed:
        raise ValueError("All coefficients are zero.")

    n = len(parsed) - 1
    cols = (n // 2) + 1
    rows = n + 1
    notes: List[str] = []

    rt = [[sp.Integer(0)] * cols for _ in range(rows)]

    for i in range(0, len(parsed), 2):
        if i // 2 < cols: rt[0][i // 2] = parsed[i]
    for i in range(1, len(parsed), 2):
        if i // 2 < cols: rt[1][i // 2] = parsed[i]

    row_power = lambda r: n - r

    for r in range(2, rows):
        if _is_zero_row_sym(rt[r - 1]):
            pw = row_power(r - 1)
            deriv = _aux_deriv_sym(rt[r - 2], pw + 1)
            notes.append(f"s^{pw}: zero row → auxiliary polynomial derivative")
            for j in range(min(len(deriv), cols)):
                rt[r - 1][j] = deriv[j]

        pivot = rt[r - 1][0]
        if _is_zero_sym(pivot) and not _is_zero_row_sym(rt[r - 1]):
            notes.append(f"s^{row_power(r-1)}: zero pivot → ε substitution")
            rt[r - 1][0] = _eps
            pivot = _eps

        if _is_zero_sym(pivot):
            rt[r - 1][0] = _eps
            pivot = _eps

        for c in range(cols - 1):
            val = (pivot * rt[r - 2][c + 1] - rt[r - 2][0] * rt[r - 1][c + 1]) / pivot
            rt[r][c] = sp.simplify(sp.cancel(val))

    first_col = [rt[r][0] for r in range(rows)]

    # ── Derive stability conditions ─────────────────────────────────────
    conditions = []
    # The first element (leading coefficient) sets the required sign
    # All elements must have the same sign as the first (assuming > 0 by convention)
    lead = first_col[0]
    for r in range(1, rows):
        expr = sp.simplify(first_col[r])
        if expr == 0:
            continue
        if expr.free_symbols:
            conditions.append(expr > 0 if lead > 0 or (hasattr(lead, "free_symbols") and not lead.free_symbols) else expr > 0)
        elif expr.is_number:
            # constant — just check it's positive
            pass

    # Assume leading coeff > 0, all first-column entries must be > 0
    all_conds = []
    for expr in first_col:
        expr_s = sp.simplify(expr)
        if expr_s.free_symbols:  # contains K
            all_conds.append(expr_s > 0)

    # Solve for K
    stable_range = None
    if all_conds:
        try:
            solution = sp.solve(all_conds, K, domain=sp.S.Reals)
            if solution is not sp.S.EmptySet:
                stable_range = str(solution)
        except Exception:
            pass

    result = RouthResult(
        table=rt, degree=n, first_column=first_col,
        stability_conditions=all_conds, stable_range=stable_range,
        notes=notes, symbolic=True,
    )
    if show:
        _print_symbolic(result)
    return result


# ═══════════════════════════════════════════════════════════════════════════
#  Pretty print — Numeric
# ═══════════════════════════════════════════════════════════════════════════

def _print_numeric(result: RouthResult) -> None:
    table = result.table
    n = result.degree
    cols = len(table[0])

    cells = [[_fmt_frac(v) for v in row] for row in table]
    col_widths = [max(len(cells[r][c]) for r in range(len(table))) for c in range(cols)]
    label_w = max(len(f"s^{n - r}") for r in range(len(table)))

    total_w = label_w + 3 + sum(col_widths) + 2 * (cols - 1)

    print()
    print("─" * (total_w + 10))
    print("  Routh–Hurwitz Table")
    print("─" * (total_w + 10))
    print()

    for r in range(len(table)):
        power = n - r
        label = f"s^{power}".rjust(label_w)
        row_str = "  ".join(cells[r][c].rjust(col_widths[c]) for c in range(cols))
        fc = result.first_column[r]
        sign_char = "+" if fc > 0 else ("−" if fc < 0 else "0")
        print(f"  {label} │ {row_str}    ({sign_char})")

    print()
    fc_str = ", ".join(_fmt_frac(v) for v in result.first_column)
    print(f"  First column : [{fc_str}]")
    print(f"  Sign changes : {result.sign_changes}")
    print()

    if result.is_stable:
        print("  ✅ STABLE — no poles in the right half-plane")
    else:
        p = "pole" if result.rhp_roots == 1 else "poles"
        print(f"  ❌ UNSTABLE — {result.rhp_roots} {p} in the right half-plane")

    if result.notes:
        print()
        for note in result.notes:
            print(f"  ⚠  {note}")
    print()


# ═══════════════════════════════════════════════════════════════════════════
#  Pretty print — Symbolic
# ═══════════════════════════════════════════════════════════════════════════

def _fmt_sym(expr) -> str:
    """Format a sympy expression compactly."""
    expr = sp.nsimplify(expr)
    if expr.is_number:
        r = sp.Rational(expr)
        if r.q == 1:
            return str(r.p)
        return f"{r.p}/{r.q}"
    # Try to express as a single fraction (a + bK)/c
    expr_f = sp.factor(expr)
    expr_c = sp.cancel(expr)
    # Pick the shortest representation
    candidates = [str(expr), str(expr_f), str(expr_c)]
    return min(candidates, key=len)


def _print_symbolic(result: RouthResult) -> None:
    table = result.table
    n = result.degree
    cols = len(table[0])

    cells = [[_fmt_sym(v) for v in row] for row in table]
    col_widths = [max(len(cells[r][c]) for r in range(len(table))) for c in range(cols)]
    label_w = max(len(f"s^{n - r}") for r in range(len(table)))

    total_w = label_w + 3 + sum(col_widths) + 2 * (cols - 1)

    print()
    print("─" * max(total_w + 10, 40))
    print("  Routh–Hurwitz Table (symbolic)")
    print("─" * max(total_w + 10, 40))
    print()

    for r in range(len(table)):
        power = n - r
        label = f"s^{power}".rjust(label_w)
        row_str = "  ".join(cells[r][c].rjust(col_widths[c]) for c in range(cols))
        print(f"  {label} │ {row_str}")

    print()

    # First column
    print("  First column:")
    for r in range(len(table)):
        power = n - r
        print(f"    s^{power}:  {_fmt_sym(result.first_column[r])}")

    # Stability conditions
    if result.stability_conditions:
        print()
        print("  Stability conditions (all must hold):")
        for cond in result.stability_conditions:
            print(f"    • {cond}")

    if result.stable_range is not None:
        print()
        range_display = result.stable_range
        try:
            # Intersect individual condition solutions as intervals
            intervals = []
            for cond in result.stability_conditions:
                sol = sp.solveset(cond, K, domain=sp.S.Reals)
                if sol is not sp.S.EmptySet and sol != sp.S.Reals:
                    intervals.append(sol)
            if intervals:
                intersection = intervals[0]
                for iv in intervals[1:]:
                    intersection = intersection.intersect(iv)
                if hasattr(intersection, 'inf') and hasattr(intersection, 'sup'):
                    lo = float(intersection.inf)
                    hi = float(intersection.sup)
                    if lo > -1e10 and hi < 1e10:
                        lo_s = f"{lo:.4g}" if lo != int(lo) else str(int(lo))
                        hi_s = f"{hi:.4g}" if hi != int(hi) else str(int(hi))
                        range_display = f"{lo_s} < K < {hi_s}"
                    else:
                        range_display = str(intersection)
                else:
                    range_display = str(intersection)
        except Exception:
            pass
        print(f"  ✅ Asymptotically stable for:  {range_display}")
    elif result.stability_conditions:
        print()
        print("  ⚠  Could not solve conditions automatically.")

    if result.notes:
        print()
        for note in result.notes:
            print(f"  ⚠  {note}")
    print()


# ═══════════════════════════════════════════════════════════════════════════
#  CLI
# ═══════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("Example: s³ − 4s² + s + 6 = 0")
    print("Roots: (s+1)(s−2)(s−3) = 0")
    routh_hurwitz([1, -4, 1, 6])

    print("\nParametric example: s⁴ + 6s³ + 11s² + 6s + (K+2) = 0")
    routh_hurwitz([1, 6, 11, 6, "K+2"])
