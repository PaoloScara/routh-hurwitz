#!/usr/bin/env python3
from routh import routh_hurwitz

def header(title: str):
    print("\n" + "═" * 64)
    print(f"  {title}")
    print("═" * 64)

# ═══════════════════════════════════════════════════════════════════════
#  NUMERIC EXAMPLES
# ═══════════════════════════════════════════════════════════════════════

header("Example 1:   s³ − 4s² + s + 6 = 0")
routh_hurwitz([1, -4, 1, 6])


header("Example 2:   2s⁴ + s³ + 3s² + 5s + 10 = 0")
routh_hurwitz([2, 1, 3, 5, 10])


header("Example 4:   s³ + 3s + 2 = 0  (missing s² → ε substitution)")
routh_hurwitz([1, 0, 3, 2])


header("Example 5:   s³ + s² + s = 0  (zero last row, pole at origin)")
print("Factored: s(s² + s + 1) = 0  →  marginally stable")
routh_hurwitz([1, 1, 1, 0])


header("Example 6:   s³ + s² = 0  (double pole at origin)")
print("Factored: s²(s + 1) = 0")
routh_hurwitz([1, 1, 0, 0])


header("Example 7:   s⁴ + s³ − 3s² − s + 2 = 0  (zero row)")
print("Auxiliary eq: −2s² + 2 = 0  →  s = ±1")
routh_hurwitz([1, 1, -3, -1, 2])


# ═══════════════════════════════════════════════════════════════════════
#  SYMBOLIC (PARAMETRIC) EXAMPLES — K left free
# ═══════════════════════════════════════════════════════════════════════

header("Example 11 (symbolic):   s⁴ + 6s³ + 11s² + 6s + (K+2) = 0")
print("Find the range of K for asymptotic stability.")
routh_hurwitz([1, 6, 11, 6, "K+2"])


header("Example 12 (symbolic):   s⁴ + 9s³ + 33s² + (25+10*K)s − 10*K = 0")
print("Find the range of K for asymptotic stability.")
routh_hurwitz([1, 9, 33, "25+10*K", "-10*K"])


# ═══════════════════════════════════════════════════════════════════════
#  PARAMETRIC with K substituted — verify specific values
# ═══════════════════════════════════════════════════════════════════════

header("Example 11:   verify specific K values")

print("  K = 5  (stable, −2 < K < 8):")
routh_hurwitz([1, 6, 11, 6, "K+2"], K_val=5)

print("  K = −5  (unstable, K < −2):")
routh_hurwitz([1, 6, 11, 6, "K+2"], K_val=-5)

print("  K = 10  (unstable, K > 8):")
routh_hurwitz([1, 6, 11, 6, "K+2"], K_val=10)

print("  K = −2  (marginally stable, pole at origin):")
routh_hurwitz([1, 6, 11, 6, "K+2"], K_val=-2)

print("  K = 8  (marginally stable, aux eq → s = ±j):")
routh_hurwitz([1, 6, 11, 6, "K+2"], K_val=8)


# ═══════════════════════════════════════════════════════════════════════

print("\n" + "═" * 64)
print("  All examples completed!")
print("═" * 64)
