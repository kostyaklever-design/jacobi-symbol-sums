# Jacobi-symbol sums S(n) = Σ ((x² − 1)/n): mean squares and oscillations

Research code and paper sources by **Kostiantyn Koshelenko** (independent researcher, the Netherlands).
ORCID [0009-0009-5184-9737](https://orcid.org/0009-0009-5184-9737) · kostiantynkoshelenko@gmail.com

## The object

For odd n let

    S(n) = Σ_{x mod n} ( (x² − 1) / n )        (Jacobi symbol)

S is multiplicative, and its sign is the Liouville function λ(n). First values for n = 1, 3, 5, 7, …:
`1, −1, −1, −1, 3, −1, −1, 1, −1, −1, 1, −1, 15, −9, …`.
Its Dirichlet series contains the factor ζ(2s−2)/ζ(2s−1)², so every zero ρ of ζ gives a pole of order ≥ 2 at (1+ρ)/2.
The family S_m(n) = Σ ((x² − m)/n) replaces one ζ by a Dirichlet L-function.

## Results

Writing Σ_{n≤x, n odd} S(n) = C·x^{3/2} + E(x):

| Result | Statement | Conditions |
|---|---|---|
| Theorem A | logarithmic mean square of E(x)/x equals **1/π²**; for squarefree m an explicit constant, independent of χ_m(p) | unconditional |
| Theorem B | mean square of a Riesz-smoothed error term normalised by X^{3/4} log X equals Λ = 7.507…·10⁻⁶ | RH, simple zeros, a moment bound for 1/ζ′(ρ) |
| Theorem C | the smoothed error term has ≫ log T sign changes of size X^{3/4}(log X)^ν in [2, T] | unconditional |
| Theorem D | **E(x) = Ω±(x log log x)** (Montgomery's method; for Euler's φ the analogue is only √(log log x)) | unconditional |
| Theorem D_m | E_m(x) = Ω±(x √(log log x)) for squarefree m ≠ 1 | unconditional |
| Theorem C_q (𝔽_q[T]) | function-field analogue: oscillations of size q^{3n/4} from the zeros of L(u, χ_m), proved via Weil's theorem | unconditional (q odd, m squarefree, deg m ≥ 3; for Ω± also L(q^{−1/2}, χ_m) ≠ 0) |

None of these results proves or disproves the Riemann hypothesis. Theorem B assumes it; the papers also note the known-type equivalence "smoothed error term = O(X^{3/4+ε}) ⇔ RH".

## Status — please read

These are **preprints, not peer-reviewed**. The work was carried out in close collaboration with an AI system (Claude, Anthropic): the author introduced S(n) and S_m(n), chose the questions and directed the work; the AI drafted proofs and text and wrote the computations. The proofs of the mean-value lemma and of Theorems A, D and D_m were gone through step by step with an independent numerical test of each step (scripts below); for Theorem C the constants and the leading coefficient were re-computed, while the lemmas it takes from Kaczorowski–Wiertelak and Kaczorowski–Pintz are applied, not re-proved. **No human specialist has checked the proofs yet.** For Theorem B (which assumes RH, simple zeros and a moment hypothesis (H)) the explicit formula was confirmed numerically to 0.1%. The function-field paper was gone through as well: the factorisation (Lemma 2) and the golden-ratio non-vanishing step were checked symbolically, and the theorem was tested on a case not in the paper (q = 7, m = T³+T+1: residues agree to 10⁻²²). All discrepancies found during these passes were errors in the checking scripts, not in the papers. Corrections and comments are welcome.

## Papers

| | Source | Zenodo |
|---|---|---|
| Integers (v4): Theorems A–D, D_m | `papers/integers/` | [10.5281/zenodo.24369317](https://doi.org/10.5281/zenodo.24369317) |
| Function fields 𝔽_q[T] (v2) | `papers/function-field/` | [10.5281/zenodo.24378507](https://doi.org/10.5281/zenodo.24378507) |

## Reproducing the checks

```
pip install -r requirements.txt
cd scripts/integers && python3 thmA_mean_square.py
```

Each script prints what it checks; run times are seconds to a few minutes (sieves go up to 10⁷).

| Script | Checks |
|---|---|
| `integers/thmA_mean_square.py` | log mean square of E(x)/x → 1/π² (sieve to 10⁷; matches the table in the paper) |
| `integers/lemma_mv_mollifier.py` | step (ii): mollifier error for 1/ζ decays like 1/y |
| `integers/lemma_mv_product.py` | step (iii): replacing the product by short sums |
| `integers/lemma_mv_diagonal.py` | step (iv): diagonal vs full mean square vs formal series |
| `integers/lemma_mv_offdiagonal.py` | step (iv): solutions of ℓk² − ℓ′k′² = Δ ≤ τ(Δ) log U |
| `integers/lemma_mv_local_factors.py` | step (v): local factors and continuity in ε |
| `integers/thmB_explicit_formula.py` | Theorem B: Λ = 7.507·10⁻⁶ from the zeros; the explicit formula (zeros + line integral on Re s = 0.7) matches R₂ from a sieve to 0.1% for X ≤ 10⁷ (~15 min) |
| `integers/thmC_first_zero.py` | Theorem C: leading coefficient B_ρ₀ = 0.001775 for the first zero, loss exponent 1.147, gap Δ* = 0.394, non-cancellation of poles |
| `integers/thmD_lemma_D0.py` | exact identity of Lemma D0 and abs g(ab²) = 2^{ω(b)} b |
| `integers/thmD_lemma_D1.py` | cancellation in Σ g(d)/d |
| `integers/thmD_lemma_D1b.py` | E(x)/x vs truncated sawtooth sums |
| `integers/thmD_lemma_D2.py` | Liouville-type bound ‖hq/(2√d)‖ ≥ 1/(6(hq)²) |
| `integers/thmD_lemma_D3.py` | average over progressions = A(q, α); growth ≈ 0.88 log z |
| `integers/thmDm_ingredients.py` | Theorem D_m: local weights, density, exact sign rule |
| `function-field/ff_jacobi.py` | 𝔽_q[T]: local formula by brute force, exact A(n), main term |
| `function-field/ff_twist2.py` | Theorem C_q for a given q and m: residues, Weil oscillations, Λ_q |
| `function-field/ff_twist.py`, `ff_twist_res.py`, `ff_double.py` | the case q = 3, m = T³−T+1, and a double-zero case (q = 5) |

## Citing

Please cite the Zenodo record (see `CITATION.cff`).

## License

Code: MIT (see `LICENSE`). Papers: © the author, under the licence stated on the Zenodo records.
