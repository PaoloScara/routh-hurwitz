#!/usr/bin/env python3
"""
Examples for Routh-Hurwitz Stability Analyzer
"""

from routh import routh_hurwitz

def main():
    print("="*80)
    print(" ROUTH-HURWITZ STABILITY ANALYZER - EXAMPLES")
    print("="*80)
    
    # Example 1: Unstable system
    print("\n" + "─"*80)
    print("EXAMPLE 1: Unstable System")
    print("─"*80)
    print("Polynomial: s³ − 4s² + s + 6 = 0")
    print("Factored: (s+1)(s−2)(s−3) = 0")
    print("Expected: UNSTABLE (2 roots in RHP)")
    routh_hurwitz([1, -4, 1, 6])
    
    # Example 2: Stable system
    print("\n" + "─"*80)
    print("EXAMPLE 2: Stable System")
    print("─"*80)
    print("Polynomial: s³ + 3s² + 3s + 1 = 0")
    print("Factored: (s+1)³ = 0")
    print("Expected: STABLE (all roots in LHP)")
    routh_hurwitz([1, 3, 3, 1])
    
    # Example 3: Parametric system
    print("\n" + "─"*80)
    print("EXAMPLE 3: Parametric System")
    print("─"*80)
    print("Polynomial: s⁴ + 6s³ + 11s² + 6s + (K+2) = 0")
    print("Expected: Stability conditions on K")
    routh_hurwitz([1, 6, 11, 6, "K+2"])
    
    # Example 4: Zero in first column (ε substitution)
    print("\n" + "─"*80)
    print("EXAMPLE 4: Zero in First Column")
    print("─"*80)
    print("Polynomial: s⁴ + s³ + 3s² + 2s + 1 = 0")
    print("Expected: ε substitution")
    routh_hurwitz([1, 1, 3, 2, 1])
    
    # Example 5: Double pole at origin (marginally stable)
    print("\n" + "─"*80)
    print("EXAMPLE 5: Double Pole at Origin")
    print("─"*80)
    print("Polynomial: s² = 0")
    print("Expected: MARGINALLY STABLE")
    routh_hurwitz([1, 0, 0])
    
    # Example 6: Another marginally stable case
    print("\n" + "─"*80)
    print("EXAMPLE 6: Poles on Imaginary Axis")
    print("─"*80)
    print("Polynomial: s⁴ + s² + 1 = 0")
    print("Expected: MARGINALLY STABLE")
    routh_hurwitz([1, 0, 1, 0, 1])
    
    # Example 7: Zero row (auxiliary polynomial)
    print("\n" + "─"*80)
    print("EXAMPLE 7: Zero Row - Auxiliary Polynomial")
    print("─"*80)
    print("Polynomial: s⁵ + 2s⁴ + 3s³ + 6s² + 2s + 4 = 0")
    routh_hurwitz([1, 2, 3, 6, 2, 4])
    
    # Example 8: Parametric evaluation
    print("\n" + "─"*80)
    print("EXAMPLE 8: Parametric with K=5")
    print("─"*80)
    print("Polynomial: s⁴ + 6s³ + 11s² + 6s + (K+2) = 0, K=5")
    routh_hurwitz([1, 6, 11, 6, "K+2"], K_val=5)
    
    # Example 9: Decimal format
    print("\n" + "─"*80)
    print("EXAMPLE 9: Decimal Format")
    print("─"*80)
    print("Polynomial: s³ + 2.5s² + 3.7s + 1.2 = 0")
    routh_hurwitz([1, 2.5, 3.7, 1.2], numeric_format='decimal', decimal_places=3)
    
    # Example 10: Fraction format (default)
    print("\n" + "─"*80)
    print("EXAMPLE 10: Fraction Format")
    print("─"*80)
    print("Polynomial: s³ + 0.5s² + 0.25s + 0.125 = 0")
    routh_hurwitz([1, 0.5, 0.25, 0.125], numeric_format='fraction')
    
    # Example 11: Complex parametric
    print("\n" + "─"*80)
    print("EXAMPLE 11: Complex Parametric")
    print("─"*80)
    print("Polynomial: s³ + 2s² + (K+1)s + K = 0")
    routh_hurwitz([1, 2, "K+1", "K"])
    
    # Example 12: High-order system
    print("\n" + "─"*80)
    print("EXAMPLE 12: Fifth-order System")
    print("─"*80)
    print("Polynomial: s⁵ + 5s⁴ + 11s³ + 15s² + 10s + 4 = 0")
    routh_hurwitz([1, 5, 11, 15, 10, 4])
    
    print("\n" + "="*80)
    print(" ALL EXAMPLES COMPLETED")
    print("="*80)


if __name__ == "__main__":
    main()
