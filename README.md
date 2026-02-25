# Routh–Hurwitz Stability Analyzer

A Python tool for analyzing polynomial stability using the **Routh–Hurwitz criterion**, fully symbolic via `sympy`.

## Quick Start

```bash
git clone https://github.com/PaoloScara/routh-hurwitz.git
cd routh-hurwitz
pip install sympy

# Run examples
python examples.py
```

## Usage

```python
from routh import routh_hurwitz

# Numeric — just pass coefficients (descending powers of s)
routh_hurwitz([1, -4, 1, 6])              # s³ − 4s² + s + 6

# Parametric with K evaluated
routh_hurwitz([1, 6, 11, 6, "K+2"], K_val=5)

# Symbolic — leave K_val=None to derive stability conditions
routh_hurwitz([1, 6, 11, 6, "K+2"])

```

### Numeric output

```
  s^3 │   1  1    (+)
  s^2 │  -4  6    (−)
  s^1 │ 5/2  0    (+)
  s^0 │   6  0    (+)

  First column : [1, -4, 5/2, 6]
  Sign changes : 2

  ❌ UNSTABLE — 2 poles in the right half-plane
```

### Symbolic ε (zero pivot)

```
  s^3 │       1  3    (+)
  s^2 │       ε  2    (0)
  s^1 │ 3 - 2/ε  0    (−)
  s^0 │       2  0    (+)

  ⚠  s^2: zero in first column → ε substitution
```

### Parametric K output

```
  s^4 │            1     11  K + 2
  s^3 │            6      6      0
  s^2 │           10  K + 2      0
  s^1 │ 24/5 - 3*K/5      0      0
  s^0 │        K + 2      0      0

  Stability conditions (all must hold):
    • 24/5 - 3*K/5 > 0
    • K + 2 > 0

  ✅ Asymptotically stable for:  -2 < K < 8
```

## API

```python
result = routh_hurwitz(coeffs, K_val=None, show=True)
```

## Special Cases

Handled automatically:

1. **Zero in first column** → symbolic ε substitution (`3 - 2/ε`), sign determined via limit ε → 0⁺
2. **Entire row of zeros** → replaced with derivative of auxiliary polynomial

## Examples

| Polynomial | Result |
|-----------|--------|
| s³ − 4s² + s + 6 | Unstable (2 RHP) |
| 2s⁴ + s³ + 3s² + 5s + 10 | Unstable (2 RHP) |
| s³ + 3s + 2 | Unstable (ε case) |
| s³ + s² + s | Marginally stable |
| s⁴ + s³ − 3s² − s + 2 | Unstable (zero row) |
| s⁴ + 6s³ + 11s² + 6s + (K+2) | **−2 < K < 8** |
| s⁴ + 9s³ + 33s² + (25+10K)s − 10K | **−1.956 < K < 0** |

## Dependencies

- **Python 3.6+**
- **sympy**

## References

- Routh, E.J. (1877). *A Treatise on the Stability of a Given State of Motion*
- Hurwitz, A. (1895). *Über die Bedingungen, unter welchen eine Gleichung nur Wurzeln mit negativen reellen Teilen besitzt*

## License

MIT — use freely for research and education.
