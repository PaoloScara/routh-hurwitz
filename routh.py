#!/usr/bin/env python3
"""
Routh–Hurwitz Stability Analyzer

Special cases:
  - Zero in first column  → symbolic ε substitution (shown as ε in output)
  - Entire row of zeros   → auxiliary polynomial derivative method

Usage:
    from routh import routh_hurwitz

    routh_hurwitz([1, -4, 1, 6])                        # numeric
    routh_hurwitz([1, 6, 11, 6, "K+2"], K_val=5)        # parametric, evaluate at K=5
    routh_hurwitz([1, 6, 11, 6, "K+2"])                  # symbolic — derive conditions on K
    routh_hurwitz([1, 0, 3, 2])                          # ε substitution (symbolic)
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Optional, Union

import sympy as sp

# Global symbols
K = sp.Symbol("K", real=True)
eps = sp.Symbol("ε", positive=True)


# ═══════════════════════════════════════════════════════════════════════════
#  Result
# ═══════════════════════════════════════════════════════════════════════════

@dataclass
class RouthResult:
    """Routh–Hurwitz analysis result."""
    table: List[List[sp.Expr]]
    degree: int
    first_column: List[sp.Expr]
    # Numeric mode (no free symbols)
    sign_changes: Optional[int] = None
    rhp_roots: Optional[int] = None
    is_stable: Optional[bool] = None
    # Symbolic mode (free symbols present)
    stability_conditions: Optional[List] = None
    stable_range: Optional[str] = None
    notes: List[str] = field(default_factory=list)
    symbolic: bool = False


# ═══════════════════════════════════════════════════════════════════════════
#  Helpers
# ═══════════════════════════════════════════════════════════════════════════

def _is_zero(expr: sp.Expr) -> bool:
    """Check if a sympy expression is zero."""
    return sp.simplify(expr) == 0


def _is_zero_row(row: List[sp.Expr]) -> bool:
    return all(_is_zero(v) for v in row)


def _has_free(expr: sp.Expr) -> bool:
    """Check if expression contains K (ignore ε)."""
    return bool(expr.free_symbols - {eps})


def _aux_derivative(prev_row: List[sp.Expr], power: int) -> List[sp.Expr]:
    """
    Given a Routh row at s^(power+1) with entries [a, b, c, ...]
    representing A(s) = a·s^power + b·s^(power-2) + ...
    return the derivative's row entries for s^(power-1).
    """
    result = []
    p = power
    for val in prev_row:
        if p < 0:
            break
        result.append(sp.expand(val * p))
        p -= 2
    return result


def _fmt(expr: sp.Expr) -> str:
    """Format a sympy expression compactly for display."""
    expr = sp.nsimplify(expr, rational=False)
    if expr.is_number:
        r = sp.Rational(expr)
        if r.q == 1:
            return str(r.p)
        return f"{r.p}/{r.q}"
    # Pick shortest among simplify/factor/cancel
    candidates = set()
    for fn in (sp.simplify, sp.factor, sp.cancel):
        try:
            candidates.add(str(fn(expr)))
        except Exception:
            pass
    candidates.add(str(expr))
    return min(candidates, key=len)


# ═══════════════════════════════════════════════════════════════════════════
#  Core
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
    parsed: List[sp.Expr] = []
    has_K = False
    for c in coeffs:
        if isinstance(c, str):
            has_K = True
            parsed.append(sp.sympify(c, locals={"K": K}))
        else:
            parsed.append(sp.Rational(c))

    # Substitute K if value provided
    if has_K and K_val is not None:
        parsed = [sp.nsimplify(v.subs(K, sp.Rational(K_val))) for v in parsed]

    # Trim leading zeros
    while parsed and _is_zero(parsed[0]):
        parsed.pop(0)
    if not parsed:
        raise ValueError("All coefficients are zero.")

    n = len(parsed) - 1  # degree
    cols = (n // 2) + 1
    rows = n + 1
    notes: List[str] = []

    # ── Build table ─────────────────────────────────────────────────────
    rt: List[List[sp.Expr]] = [[sp.Integer(0)] * cols for _ in range(rows)]

    for i in range(0, len(parsed), 2):
        if i // 2 < cols:
            rt[0][i // 2] = parsed[i]
    for i in range(1, len(parsed), 2):
        if i // 2 < cols:
            rt[1][i // 2] = parsed[i]

    def row_power(r: int) -> int:
        return n - r

    for r in range(2, rows):
        # Zero row → auxiliary polynomial derivative
        if _is_zero_row(rt[r - 1]):
            pw = row_power(r - 1)
            deriv = _aux_derivative(rt[r - 2], pw + 1)
            notes.append(f"s^{pw}: zero row → auxiliary polynomial derivative")
            for j in range(min(len(deriv), cols)):
                rt[r - 1][j] = deriv[j]
            for j in range(len(deriv), cols):
                rt[r - 1][j] = sp.Integer(0)

        # Zero pivot → symbolic ε substitution
        pivot = rt[r - 1][0]
        if _is_zero(pivot) and not _is_zero_row(rt[r - 1]):
            notes.append(f"s^{row_power(r-1)}: zero in first column → ε substitution")
            rt[r - 1][0] = eps
            pivot = eps
        elif _is_zero(pivot):
            rt[r - 1][0] = eps
            pivot = eps

        # Compute row
        for c in range(cols - 1):
            val = (pivot * rt[r - 2][c + 1] - rt[r - 2][0] * rt[r - 1][c + 1]) / pivot
            rt[r][c] = sp.cancel(sp.expand(val))

    # ── First column ────────────────────────────────────────────────────
    first_col = [rt[r][0] for r in range(rows)]

    # Determine mode: symbolic if any first-column entry has K
    is_symbolic = any(_has_free(e) for e in first_col)

    # ── Numeric analysis ────────────────────────────────────────────────
    if not is_symbolic:
        # Evaluate ε → 0⁺ for sign analysis
        signs = []
        last_sign = None
        for expr in first_col:
            # Take limit ε → 0⁺
            val = sp.limit(expr, eps, 0, "+") if eps in expr.free_symbols else expr
            if val == sp.oo or val > 0:
                s = 1
            elif val == -sp.oo or val < 0:
                s = -1
            elif val == 0:
                s = last_sign if last_sign is not None else 1
            else:
                s = last_sign if last_sign is not None else 1
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

    # ── Symbolic analysis ───────────────────────────────────────────────
    # All first-column entries with K must be > 0 (assuming leading coeff > 0)
    all_conds = []
    for expr in first_col:
        expr_s = sp.simplify(expr)
        if _has_free(expr_s):
            all_conds.append(expr_s > 0)

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

    cells = [[_fmt(v) for v in row] for row in table]
    col_widths = [max(len(cells[r][c]) for r in range(len(table))) for c in range(cols)]
    label_w = max(len(f"s^{n - r}") for r in range(len(table)))
    total_w = label_w + 3 + sum(col_widths) + 2 * (cols - 1)

    # Compute signs via limit ε → 0⁺
    sign_chars = []
    last_sign = None
    for expr in result.first_column:
        val = sp.limit(expr, eps, 0, "+") if eps in expr.free_symbols else expr
        if val == sp.oo or val > 0:
            ch, last_sign = "+", 1
        elif val == -sp.oo or val < 0:
            ch, last_sign = "−", -1
        elif val == 0:
            ch = "0"
        else:
            ch = "?"
        sign_chars.append(ch)

    print()
    print("─" * (total_w + 10))
    print("  Routh–Hurwitz Table")
    print("─" * (total_w + 10))
    print()

    for r in range(len(table)):
        power = n - r
        label = f"s^{power}".rjust(label_w)
        row_str = "  ".join(cells[r][c].rjust(col_widths[c]) for c in range(cols))
        print(f"  {label} │ {row_str}    ({sign_chars[r]})")

    print()
    fc_str = ", ".join(_fmt(v) for v in result.first_column)
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

def _print_symbolic(result: RouthResult) -> None:
    table = result.table
    n = result.degree
    cols = len(table[0])

    cells = [[_fmt(v) for v in row] for row in table]
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
    print("  First column:")
    for r in range(len(table)):
        power = n - r
        print(f"    s^{power}:  {_fmt(result.first_column[r])}")

    if result.stability_conditions:
        print()
        print("  Stability conditions (all must hold):")
        for cond in result.stability_conditions:
            print(f"    • {cond}")

    if result.stable_range is not None:
        print()
        range_display = result.stable_range
        try:
            intervals = []
            for cond in result.stability_conditions:
                sol = sp.solveset(cond, K, domain=sp.S.Reals)
                if sol is not sp.S.EmptySet and sol != sp.S.Reals:
                    intervals.append(sol)
            if intervals:
                intersection = intervals[0]
                for iv in intervals[1:]:
                    intersection = intersection.intersect(iv)
                if hasattr(intersection, "inf") and hasattr(intersection, "sup"):
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

    print("\nParametric: s⁴ + 6s³ + 11s² + 6s + (K+2) = 0")
    routh_hurwitz([1, 6, 11, 6, "K+2"])
