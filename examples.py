#!/usr/bin/env python3
"""
Example usage of the Routh-Hurwitz stability analyzer.
Run: python examples.py
"""

from routh import routh_hurwitz, format_routh_table


def example_1_stable():
    """Example 1: Stable system (all poles in left half-plane)"""
    print("=" * 60)
    print("EXAMPLE 1: Stable System")
    print("=" * 60)
    print("Polynomial: s³ + 3s² + 3s + 1")
    print()

    coeffs = [1, 3, 3, 1]
    result = routh_hurwitz(coeffs)

    print(format_routh_table(result.table, degree=3))
    print()
    print(f"First column: {result.first_column}")
    print(f"Sign changes: {result.sign_changes}")
    print(f"Unstable poles (RHP roots): {result.rhp_roots}")
    print(f"Verdict: {'✓ STABLE' if result.rhp_roots == 0 else '✗ UNSTABLE'}")
    if result.notes:
        print("\nNotes:")
        for note in result.notes:
            print(f"  - {note}")
    print()


def example_2_unstable():
    """Example 2: Unstable system (poles in right half-plane)"""
    print("=" * 60)
    print("EXAMPLE 2: Unstable System")
    print("=" * 60)
    print("Polynomial: s³ - s² + 2s - 2")
    print()

    coeffs = [1, -1, 2, -2]
    result = routh_hurwitz(coeffs)

    print(format_routh_table(result.table, degree=3))
    print()
    print(f"First column: {result.first_column}")
    print(f"Sign changes: {result.sign_changes}")
    print(f"Unstable poles (RHP roots): {result.rhp_roots}")
    print(f"Verdict: {'✓ STABLE' if result.rhp_roots == 0 else '✗ UNSTABLE'}")
    if result.notes:
        print("\nNotes:")
        for note in result.notes:
            print(f"  - {note}")
    print()


def example_3_marginal():
    """Example 3: Marginally stable (oscillatory, boundary case)"""
    print("=" * 60)
    print("EXAMPLE 3: Marginally Stable System")
    print("=" * 60)
    print("Polynomial: s³ + 2s² + 2s + 1")
    print()

    coeffs = [1, 2, 2, 1]
    result = routh_hurwitz(coeffs)

    print(format_routh_table(result.table, degree=3))
    print()
    print(f"First column: {result.first_column}")
    print(f"Sign changes: {result.sign_changes}")
    print(f"Unstable poles (RHP roots): {result.rhp_roots}")
    print(f"Verdict: {'✓ STABLE' if result.rhp_roots == 0 else '✗ UNSTABLE'}")
    if result.notes:
        print("\nNotes:")
        for note in result.notes:
            print(f"  - {note}")
    print()


def example_4_higher_order():
    """Example 4: Higher order stable system"""
    print("=" * 60)
    print("EXAMPLE 4: 4th Order Stable System")
    print("=" * 60)
    print("Polynomial: s⁴ + 4s³ + 7s² + 6s + 2")
    print()

    coeffs = [1, 4, 7, 6, 2]
    result = routh_hurwitz(coeffs)

    print(format_routh_table(result.table, degree=4))
    print()
    print(f"First column: {result.first_column}")
    print(f"Sign changes: {result.sign_changes}")
    print(f"Unstable poles (RHP roots): {result.rhp_roots}")
    print(f"Verdict: {'✓ STABLE' if result.rhp_roots == 0 else '✗ UNSTABLE'}")
    if result.notes:
        print("\nNotes:")
        for note in result.notes:
            print(f"  - {note}")
    print()


def example_5_with_special_case():
    """Example 5: System with special case (zero in first column)"""
    print("=" * 60)
    print("EXAMPLE 5: Special Case - Zero in First Column")
    print("=" * 60)
    print("Polynomial: s⁴ + s³ + 2s² + 2s + 1")
    print()

    coeffs = [1, 1, 2, 2, 1]
    result = routh_hurwitz(coeffs)

    print(format_routh_table(result.table, degree=4))
    print()
    print(f"First column: {result.first_column}")
    print(f"Sign changes: {result.sign_changes}")
    print(f"Unstable poles (RHP roots): {result.rhp_roots}")
    print(f"Verdict: {'✓ STABLE' if result.rhp_roots == 0 else '✗ UNSTABLE'}")
    if result.notes:
        print("\nNotes:")
        for note in result.notes:
            print(f"  - {note}")
    print()


def main():
    """Run all examples"""
    example_1_stable()
    example_2_unstable()
    example_3_marginal()
    example_4_higher_order()
    example_5_with_special_case()

    print("=" * 60)
    print("All examples completed!")
    print("=" * 60)


if __name__ == "__main__":
    main()
