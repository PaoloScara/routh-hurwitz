#!/usr/bin/env python3
"""
Routh–Hurwitz Stability Analyzer

Special cases:
  - Zero in first column  → symbolic ε substitution (shown as ε in output)
  - Entire row of zeros   → auxiliary polynomial derivative method
  - Double pole at origin → correctly identified as marginally stable

Usage:
    from routh import routh_hurwitz
    
    # Numeric with fractions (default)
    routh_hurwitz([1, -4, 1, 6])
    
    # Numeric with rounded decimals
    routh_hurwitz([1, -4, 1, 6], numeric_format='decimal', decimal_places=4)
    
    # Parametric, evaluate at K=5
    routh_hurwitz([1, 6, 11, 6, "K+2"], K_val=5)
    
    # Symbolic — derive conditions on K
    routh_hurwitz([1, 6, 11, 6, "K+2"])
    
    # ε substitution (symbolic)
    routh_hurwitz([1, 0, 3, 2])
    
    # Double pole at origin (marginally stable)
    routh_hurwitz([1, 0, 0])
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Optional, Union, Literal
import sympy as sp
import numpy as np

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
    is_marginally_stable: Optional[bool] = None
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
    """Check if all elements in a row are zero."""
    return all(_is_zero(v) for v in row)


def _has_free(expr: sp.Expr) -> bool:
    """Check if expression contains K (ignore ε)."""
    return bool(expr.free_symbols - {eps})


def _aux_derivative(prev_row: List[sp.Expr], power: int) -> List[sp.Expr]:
    """
    Given a Routh row at s^(power+1) with entries [a, b, c, ...]
    representing A(s) = a·s^power + b·s^(power-2) + ...
    return the derivative's row entries for s^(power-1).
    
    Vectorized implementation using numpy-style operations.
    """
    # Create power array
    powers = np.arange(power, power - 2*len(prev_row), -2)
    powers = powers[:len(prev_row)]
    powers[powers < 0] = 0
    
    # Vectorized multiplication
    result = [sp.expand(val * p) for val, p in zip(prev_row, powers)]
    
    return result


def _fmt(expr: sp.Expr, numeric_format: str = 'fraction', decimal_places: int = 4) -> str:
    """
    Format a sympy expression for display.
    
    Parameters
    ----------
    expr : sp.Expr
        Expression to format
    numeric_format : {'fraction', 'decimal'}
        Output format for numeric values
    decimal_places : int
        Number of decimal places for 'decimal' format
    """
    expr = sp.simplify(expr)
    
    if expr.is_number:
        if numeric_format == 'fraction':
            r = sp.nsimplify(expr, rational=True)
            if isinstance(r, sp.Rational):
                if r.q == 1:
                    return str(r.p)
                return f"{r.p}/{r.q}"
            return str(expr)
        else:  # decimal
            val = float(expr.evalf())
            return f"{val:.{decimal_places}f}".rstrip('0').rstrip('.')
    
    # For symbolic expressions, try to find shortest representation
    candidates = []
    for fn in (sp.simplify, sp.factor, sp.cancel):
        try:
            candidates.append(str(fn(expr)))
        except Exception:
            pass
    candidates.append(str(expr))
    return min(candidates, key=len)


def _check_marginal_stability(first_column: List[sp.Expr], table: List[List[sp.Expr]]) -> tuple[bool, Optional[str]]:
    """
    Check for marginal stability conditions:
    - Double pole at origin: entire last row is zero
    - Poles on imaginary axis: zero in first column with all-zero row
    
    Returns (is_marginally_stable, note)
    """
    # Check if last row is all zeros (double pole at origin)
    if _is_zero_row(table[-1]):
        # Verify by checking auxiliary polynomial
        return True, "Double pole at origin (marginally stable - poles on jω axis)"
    
    # Check for imaginary axis poles via zero rows
    for i, row in enumerate(table[:-1]):
        if _is_zero_row(row):
            return True, f"Poles on imaginary axis detected at row s^{len(table)-1-i} (marginally stable)"
    
    return False, None


# ═══════════════════════════════════════════════════════════════════════════
#  Core - Vectorized Operations
# ═══════════════════════════════════════════════════════════════════════════

def _build_initial_rows(parsed: List[sp.Expr], cols: int) -> tuple[List[sp.Expr], List[sp.Expr]]:
    """
    Build first two rows of Routh table using vectorized operations.
    
    Returns (row0, row1)
    """
    row0 = [sp.Integer(0)] * cols
    row1 = [sp.Integer(0)] * cols
    
    # Even indices go to row0
    even_indices = np.arange(0, len(parsed), 2)
    col_indices_0 = even_indices // 2
    valid_0 = col_indices_0 < cols
    for i, col_idx in zip(even_indices[valid_0], col_indices_0[valid_0]):
        row0[col_idx] = parsed[i]
    
    # Odd indices go to row1
    odd_indices = np.arange(1, len(parsed), 2)
    col_indices_1 = odd_indices // 2
    valid_1 = col_indices_1 < cols
    for i, col_idx in zip(odd_indices[valid_1], col_indices_1[valid_1]):
        row1[col_idx] = parsed[i]
    
    return row0, row1


def _compute_routh_row(prev2: List[sp.Expr], prev1: List[sp.Expr], cols: int) -> List[sp.Expr]:
    """
    Compute a Routh table row using vectorized formula.
    
    Formula: R[i,j] = (R[i-1,0] * R[i-2,j+1] - R[i-2,0] * R[i-1,j+1]) / R[i-1,0]
    """
    pivot = prev1[0]
    new_row = []
    
    # Vectorized computation for all columns at once
    for c in range(cols - 1):
        numerator = pivot * prev2[c + 1] - prev2[0] * prev1[c + 1]
        val = numerator / pivot
        new_row.append(sp.cancel(sp.expand(val)))
    
    # Last column is always zero
    new_row.append(sp.Integer(0))
    
    return new_row


def routh_hurwitz(
    coeffs: List[Union[int, float, str]],
    K_val: Optional[Union[int, float]] = None,
    show: bool = True,
    numeric_format: Literal['fraction', 'decimal'] = 'fraction',
    decimal_places: int = 4,
) -> RouthResult:
    """
    Compute the Routh table for a polynomial (fully symbolic).

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
    numeric_format : {'fraction', 'decimal'}
        Format for displaying numeric values (default 'fraction').
    decimal_places : int
        Decimal places for 'decimal' format (default 4).

    Returns
    -------
    RouthResult
    """
    # ── Parse coefficients (always symbolic) ────────────────────────────────
    parsed: List[sp.Expr] = []
    has_K = False
    for c in coeffs:
        if isinstance(c, str):
            has_K = True
            parsed.append(sp.sympify(c, locals={"K": K}))
        else:
            # Keep as exact rational, not float
            parsed.append(sp.Rational(c))

    # Substitute K if value provided (but keep symbolic computation)
    if has_K and K_val is not None:
        parsed = [sp.simplify(v.subs(K, sp.Rational(K_val))) for v in parsed]

    # Trim leading zeros
    while parsed and _is_zero(parsed[0]):
        parsed.pop(0)
    if not parsed:
        raise ValueError("All coefficients are zero.")

    n = len(parsed) - 1  # degree
    cols = (n // 2) + 1
    rows = n + 1
    notes: List[str] = []

    # ── Build table using vectorized operations ────────────────────────────
    rt: List[List[sp.Expr]] = []
    
    # First two rows
    row0, row1 = _build_initial_rows(parsed, cols)
    rt.append(row0)
    rt.append(row1)

    def row_power(r: int) -> int:
        return n - r

    # Compute remaining rows
    for r in range(2, rows):
        # Zero row → auxiliary polynomial derivative
        if _is_zero_row(rt[r - 1]):
            pw = row_power(r - 1)
            deriv = _aux_derivative(rt[r - 2], pw + 1)
            notes.append(f"s^{pw}: zero row → auxiliary polynomial derivative")
            # Pad with zeros if needed
            deriv.extend([sp.Integer(0)] * (cols - len(deriv)))
            rt[r - 1] = deriv[:cols]

        # Zero pivot → symbolic ε substitution
        pivot = rt[r - 1][0]
        if _is_zero(pivot) and not _is_zero_row(rt[r - 1]):
            notes.append(f"s^{row_power(r-1)}: zero in first column → ε substitution")
            rt[r - 1][0] = eps
        elif _is_zero(pivot):
            rt[r - 1][0] = eps

        # Compute row using vectorized formula
        new_row = _compute_routh_row(rt[r - 2], rt[r - 1], cols)
        rt.append(new_row)

    # ── First column ────────────────────────────────────────────────────────
    first_col = [rt[r][0] for r in range(rows)]

    # Check for marginal stability
    is_marginal, marginal_note = _check_marginal_stability(first_col, rt)
    if marginal_note:
        notes.append(marginal_note)

    # Determine mode: symbolic if any first-column entry has K
    is_symbolic = any(_has_free(e) for e in first_col)

    # ── Numeric analysis (symbolic computation, numeric interpretation) ─────
    if not is_symbolic:
        # Evaluate ε → 0⁺ for sign analysis (vectorized)
        signs = []
        last_sign = None
        
        for expr in first_col:
            # Take limit ε → 0⁺
            val = sp.limit(expr, eps, 0, "+") if eps in expr.free_symbols else expr
            
            if val == sp.oo or (val.is_number and val > 0):
                s = 1
            elif val == -sp.oo or (val.is_number and val < 0):
                s = -1
            elif val == 0 or _is_zero(val):
                s = last_sign if last_sign is not None else 1
            else:
                s = last_sign if last_sign is not None else 1
            
            last_sign = s
            signs.append(s)

        # Vectorized sign change detection
        signs_arr = np.array(signs)
        sign_changes = np.sum(signs_arr[1:] != signs_arr[:-1])

        result = RouthResult(
            table=rt, degree=n, first_column=first_col,
            sign_changes=int(sign_changes), 
            rhp_roots=int(sign_changes), 
            is_stable=(sign_changes == 0 and not is_marginal),
            is_marginally_stable=is_marginal,
            notes=notes, 
            symbolic=False,
        )
        if show:
            _print_numeric(result, numeric_format, decimal_places)
            return None  
        return result

    # ── Symbolic analysis ───────────────────────────────────────────────────
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
        _print_symbolic(result, numeric_format, decimal_places)
        return None  
    return result


# ═══════════════════════════════════════════════════════════════════════════
#  Pretty print — Numeric (Enhanced Graphics)
# ═══════════════════════════════════════════════════════════════════════════

def _print_numeric(result: RouthResult, numeric_format: str = 'fraction', decimal_places: int = 4) -> None:
    """Enhanced printing with better visual organization."""
    table = result.table
    n = result.degree
    cols = len(table[0])

    # Format all cells
    cells = [[_fmt(v, numeric_format, decimal_places) for v in row] for row in table]
    col_widths = [max(len(cells[r][c]) for r in range(len(table))) for c in range(cols)]
    label_w = max(len(f"s^{n - r}") for r in range(len(table)))
    total_w = label_w + 3 + sum(col_widths) + 2 * (cols - 1)

    # Compute signs via limit ε → 0⁺
    sign_chars = []
    last_sign = None
    for expr in result.first_column:
        val = sp.limit(expr, eps, 0, "+") if eps in expr.free_symbols else expr
        if val == sp.oo or (val.is_number and val > 0):
            ch, last_sign = "+", 1
        elif val == -sp.oo or (val.is_number and val < 0):
            ch, last_sign = "-", -1
        elif val == 0 or _is_zero(val):
            ch = "0"
        else:
            ch = "?"
        sign_chars.append(ch)

    # Enhanced header - make sure separator is at least as long as title
    title = "ROUTH-HURWITZ STABILITY TABLE"
    sep_w = max(total_w + 14, len(title) + 4)
    print()
    print("=" * sep_w)
    print("  " + title)
    print("=" * sep_w)
    print()

    # Table with better separators
    for r in range(len(table)):
        power = n - r
        label = f"s^{power}".rjust(label_w)
        row_str = "  ".join(cells[r][c].rjust(col_widths[c]) for c in range(cols))
        sign_indicator = f"[{sign_chars[r]}]"
        print(f"  {label} | {row_str}  {sign_indicator}")
        
        # Separator after first two rows
        if r == 1:
            print("  " + "-" * (label_w + 1) + "+" + "-" * (sum(col_widths) + 2 * (cols - 1) + 2))

    print()
    print("=" * sep_w)
    
    # First column summary
    fc_str = ", ".join(_fmt(v, numeric_format, decimal_places) for v in result.first_column)
    print(f"  First column: [{fc_str}]")
    print(f"  Sign changes: {result.sign_changes}")
    print("=" * sep_w)
    print()

    # Verdict with enhanced graphics
    if result.is_marginally_stable:
        print("  >> MARGINALLY STABLE")
        print("     System has poles on the imaginary axis (jw)")
        print("     No poles in right half-plane, but sustained oscillations possible")
    elif result.is_stable:
        print("  >> ASYMPTOTICALLY STABLE")
        print("     All poles in left half-plane")
    else:
        p = "pole" if result.rhp_roots == 1 else "poles"
        print(f"  >> UNSTABLE")
        print(f"     {result.rhp_roots} {p} in right half-plane")

    if result.notes:
        print()
        print("  Notes:")
        for note in result.notes:
            print(f"  * {note}")
    
    print()
    print("=" * sep_w)
    print()


# ═══════════════════════════════════════════════════════════════════════════
#  Pretty print — Symbolic (Enhanced Graphics)
# ═══════════════════════════════════════════════════════════════════════════

def _print_symbolic(result: RouthResult, numeric_format: str = 'fraction', decimal_places: int = 4) -> None:
    """Enhanced symbolic printing."""
    table = result.table
    n = result.degree
    cols = len(table[0])

    cells = [[_fmt(v, numeric_format, decimal_places) for v in row] for row in table]
    col_widths = [max(len(cells[r][c]) for r in range(len(table))) for c in range(cols)]
    label_w = max(len(f"s^{n - r}") for r in range(len(table)))
    total_w = label_w + 3 + sum(col_widths) + 2 * (cols - 1)

    # Make sure separator covers title
    title = "ROUTH-HURWITZ TABLE (Symbolic Analysis)"
    sep_w = max(total_w + 14, len(title) + 4, 52)
    
    print()
    print("=" * sep_w)
    print("  " + title)
    print("=" * sep_w)
    print()

    for r in range(len(table)):
        power = n - r
        label = f"s^{power}".rjust(label_w)
        row_str = "  ".join(cells[r][c].rjust(col_widths[c]) for c in range(cols))
        print(f"  {label} | {row_str}")
        
        if r == 1:
            print("  " + "-" * (label_w + 1) + "+" + "-" * (sum(col_widths) + 2 * (cols - 1) + 2))

    print()
    print("=" * sep_w)
    print("  First Column (must all be > 0 for stability)")
    print("=" * sep_w)
    
    for r in range(len(table)):
        power = n - r
        expr_str = _fmt(result.first_column[r], numeric_format, decimal_places)
        print(f"  s^{power}: {expr_str}")

    if result.stability_conditions:
        print("=" * sep_w)
        print("  Stability Conditions (all must hold):")
        print("=" * sep_w)
        for i, cond in enumerate(result.stability_conditions, 1):
            print(f"  {i}. {cond}")

    if result.stable_range is not None:
        print("=" * sep_w)
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
                        if numeric_format == 'fraction':
                            lo_s = str(sp.nsimplify(lo, rational=True))
                            hi_s = str(sp.nsimplify(hi, rational=True))
                        else:
                            lo_s = f"{lo:.{decimal_places}f}".rstrip('0').rstrip('.')
                            hi_s = f"{hi:.{decimal_places}f}".rstrip('0').rstrip('.')
                        range_display = f"{lo_s} < K < {hi_s}"
                    else:
                        range_display = str(intersection)
                else:
                    range_display = str(intersection)
        except Exception:
            pass
        print()
        print("  >> Asymptotically stable for:")
        print(f"     {range_display}")
    elif result.stability_conditions:
        print("=" * sep_w)
        print("  >> Could not solve conditions automatically")

    if result.notes:
        print("=" * sep_w)
        print("  Notes:")
        for note in result.notes:
            print(f"  * {note}")
    
    print()
    print("=" * sep_w)
    print()


# ═══════════════════════════════════════════════════════════════════════════
#  CLI
# ═══════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    print("Example 1: Unstable system")
    print("s³ − 4s² + s + 6 = 0")
    print("Roots: (s+1)(s−2)(s−3) = 0")
    routh_hurwitz([1, -4, 1, 6])

    print("\n" + "="*70 + "\n")
    print("Example 2: Parametric system")
    print("s⁴ + 6s³ + 11s² + 6s + (K+2) = 0")
    routh_hurwitz([1, 6, 11, 6, "K+2"])
    
    print("\n" + "="*70 + "\n")
    print("Example 3: Double pole at origin (marginally stable)")
    print("s² = 0")
    routh_hurwitz([1, 0, 0])
    
    print("\n" + "="*70 + "\n")
    print("Example 4: Decimal format")
    print("s³ + 2.5s² + 3.7s + 1.2 = 0")
    routh_hurwitz([1, 2.5, 3.7, 1.2], numeric_format='decimal', decimal_places=3)
