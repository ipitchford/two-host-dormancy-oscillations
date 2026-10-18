# Agent navigation

## Statement and status
Read [README.md](README.md), [ASSURANCE.md](ASSURANCE.md), and `manuscript.md` before reusing a claim. This is an anonymous written-proof candidate with computational supporting certificates, not a formally verified theorem or human peer-reviewed article.

## Claim-to-evidence map

| Claim | Inspect | Trust boundary |
|---|---|---|
| Equilibrium equations and Jacobian | [model.py](model.py), [verify_symbolic.py](verify_symbolic.py) | Model transcription must match manuscript/source equations |
| Exact counterexample and resident stability | [exact_counterexample.js](exact_counterexample.js), [family_symbolic_certificate.json](family_symbolic_certificate.json) | Finite rational identities plus written Routh-Hurwitz interpretation |
| All C >= 1/2 instability | [all_c_hurwitz.js](all_c_hurwitz.js), [verify_symbolic.py](verify_symbolic.py) | Coefficient sign certificate after reparameterization; not a numerical grid |
| Degree-17 Hopf polynomial | [hopf_polynomial.js](hopf_polynomial.js), [hopf_polynomial_coefficients.json](hopf_polynomial_coefficients.json) | Exact polynomial identities |
| Local supercritical Hopf | [certify_hopf.py](certify_hopf.py), [interval_arithmetic.py](interval_arithmetic.py), [hopf_interval_certificate.json](hopf_interval_certificate.json) | Outward-rounded integer intervals plus written bifurcation theorem |
| Almost-everywhere nonconvergence | `manuscript.md` | Stable-manifold argument, not orbit sampling |

## Reproduction and falsification
Use [run_checks.py](run_checks.py) with pinned [requirements.txt](requirements.txt), both normally and under `-O`. [negative_controls.py](negative_controls.py) rejects a corrupted Hopf polynomial, altered vector field and inverted interval. Logs are in `reproduction/`; publication-time retained checks are summarized in `REPLAY_RECEIPT.md`.

## Safe reuse
Invasion Schur complements determine signs but are not eigenvalue magnitudes. Local periodic solutions do not establish a global attractor or convergence of almost every trajectory to a cycle. Consult [REVIEW_RESPONSE.md](REVIEW_RESPONSE.md), [PROVENANCE.md](PROVENANCE.md), and [LICENSES.md](LICENSES.md). No historic search completeness, independent empirical reproduction or universal formal verification is claimed.
