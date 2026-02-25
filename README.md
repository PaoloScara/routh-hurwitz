# Routh–Hurwitz Stability Analyzer

A Python tool for analyzing the stability of linear control systems using the classical **Routh–Hurwitz criterion**.

## What It Does

Given a polynomial (characteristic equation of a control system), it:
- **Builds the Routh table** automatically
- **Counts sign changes** in the first column
- **Determines the number of unstable poles** (right half-plane roots)
- **Handles special cases**: zero rows, zero first columns with epsilon substitution
- **Provides detailed notes** on computation steps

## Key Features

✅ **No dependencies** – Pure Python 3, no external packages  
✅ **Automatic special case handling** – Auxiliary polynomial method for zero rows  
✅ **Clear output formatting** – Easy-to-read Routh table display  
✅ **Numerical stability** – Configurable tolerance for numerical errors  
✅ **Well-documented** – Includes examples and detailed notes on computation  

## Installation

```bash
git clone https://github.com/YOUR_USERNAME/routh-hurwitz.git
cd routh-hurwitz
```

**No dependencies to install!** Just clone and run.

## Quick Start

### Run Examples

```bash
python examples.py
```

### Use in Your Code

```python
from routh import routh_hurwitz, format_routh_table

# Polynomial: s³ + 2s² + 3s + 4
coeffs = [1, 2, 3, 4]

result = routh_hurwitz(coeffs)

# Display the table
print(format_routh_table(result.table, degree=3))

# Check stability
if result.rhp_roots == 0:
    print("✓ System is STABLE")
else:
    print(f"✗ System is UNSTABLE ({result.rhp_roots} unstable poles)")
```

## Usage

### Basic Function

```python
result = routh_hurwitz(coeffs)
```

**Parameters:**
- `coeffs` (list): Polynomial coefficients in **descending powers**
  - Example: `[1, 2, 3, 4]` represents `s³ + 2s² + 3s + 4`
- `tol` (float, default=1e-12): Tolerance for treating values as zero
- `epsilon` (float, default=1e-6): Small positive value for zero substitution
- `use_symbolic_epsilon` (bool): Whether to mark epsilon symbolically

**Returns:** `RouthResult` with:
- `table` – Complete Routh array (list of lists)
- `first_column` – Leftmost column (used for stability analysis)
- `sign_changes` – Number of sign changes in first column
- `rhp_roots` – **Number of poles in right half-plane** ← Use this!
- `notes` – Messages about special cases encountered

### Stability Criterion

- **Stable system**: All roots in left half-plane → `rhp_roots == 0`
- **Unstable system**: Any roots in right half-plane → `rhp_roots > 0`
- Number of unstable poles = number of sign changes in first column

### Example: Check Stability

```python
from routh import routh_hurwitz

coeffs = [1, 3, 3, 1]  # s³ + 3s² + 3s + 1
result = routh_hurwitz(coeffs)

print(f"Unstable poles: {result.rhp_roots}")
print(f"Status: {'STABLE ✓' if result.rhp_roots == 0 else 'UNSTABLE ✗'}")
```

## How It Works

### The Routh–Hurwitz Criterion

For a polynomial `P(s) = aₙsⁿ + aₙ₋₁sⁿ⁻¹ + ... + a₀`:

1. **Build the Routh table** using recursive formula
2. **Count sign changes** in the first column
3. **Result**: Number of sign changes = Number of poles in right half-plane

### Special Cases Handled

#### Case 1: Zero in First Column
When the first element of a row is zero (but row is not all zeros):
- **Solution**: Substitute small epsilon value
- **Result**: Allows computation to continue

#### Case 2: Entire Row of Zeros
When all elements of a row become zero:
- **Solution**: Use auxiliary polynomial from previous row and differentiate
- **Result**: Extracts information about marginally stable modes

## Examples

### Example 1: Stable System

```python
from routh import routh_hurwitz, format_routh_table

coeffs = [1, 3, 3, 1]  # s³ + 3s² + 3s + 1
result = routh_hurwitz(coeffs)

print(format_routh_table(result.table, degree=3))
# Output:
# s^3  |            1            3
# s^2  |            3            1
# s^1  |      2.66667            0
# s^0  |            1            0
#
# Sign changes: 0
# Unstable poles: 0
# ✓ STABLE
```

### Example 2: Unstable System

```python
coeffs = [1, -1, 2, -2]  # s³ - s² + 2s - 2
result = routh_hurwitz(coeffs)

print(f"RHP roots: {result.rhp_roots}")
# Output: RHP roots: 1
# ✗ UNSTABLE (1 unstable pole)
```

### Example 3: With Special Cases

```python
coeffs = [1, 1, 2, 2, 1]  # s⁴ + s³ + 2s² + 2s + 1
result = routh_hurwitz(coeffs)

print(result.notes)
# Output: Notes about zero handling, if any
```

## Command Line

Run the script directly:

```bash
python routh.py
```

This runs a default example: `s³ + 2s² + 3s + 4`

## Common Use Cases

### Control System Design
Quickly verify if a closed-loop transfer function is stable:

```python
# Closed-loop denominator
coeffs = [1, 10, 5, 2]  # D(s) = s³ + 10s² + 5s + 2
result = routh_hurwitz(coeffs)
assert result.rhp_roots == 0, "System must be stable!"
```

### Academic Research
Analyze polynomial stability in networked systems, power systems, manufacturing:

```python
for test_case in test_cases:
    result = routh_hurwitz(test_case.polynomial)
    print(f"Test {test_case.name}: {result.rhp_roots} unstable poles")
```

### Educational
Learn the Routh–Hurwitz method step-by-step:

```bash
python examples.py  # See how the table is built
```

## Files

```
routh-hurwitz/
├── routh.py              # Main module
├── examples.py           # 5 example cases
├── README.md             # This file
├── requirements.txt      # (empty - no dependencies)
└── .gitignore           # Python cache files
```

## Testing

Run all examples:

```bash
python examples.py
```

Expected output: 5 examples showing stable, unstable, and special cases.

## Mathematical Background

The Routh–Hurwitz criterion states:

> **The number of roots of P(s) in the right half-plane equals the number of sign changes in the first column of the Routh table.**

This avoids computing roots explicitly and is useful for polynomials where numerical root-finding is unstable.

## Numerical Considerations

- **Tolerance** (`tol`): Set higher for ill-conditioned polynomials
- **Epsilon** (`epsilon`): Used when first column element is zero; affects precision
- All computation uses floating-point arithmetic

## License

Open source – use freely for research and education.

## References

- Routh, E. J. (1877). "A Treatise on the Stability of a Given State of Motion"
- Hurwitz, A. (1895). "Über die Bedingungen, unter welchen eine Gleichung nur Wurzeln mit negativen reellen Teilen besitzt"

## Support

For issues or questions:
1. Check `examples.py` for usage patterns
2. Review function docstrings in `routh.py`
3. Run with your own polynomial and inspect `result.notes`

---

**Made for control engineers and researchers who need quick stability analysis.**
