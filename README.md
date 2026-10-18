# Two-host dormancy: exact instability and local stable oscillations

Scholarly creator: Anonymous. Written-proof candidate with exact symbolic and outward-rounded interval certificates; not externally peer reviewed or formally verified.

The manuscript constructs an admissible two-host dormancy family for every competition parameter C >= 1/2. Each resident equilibrium is stable within its resident subsystem and invasible by the missing host, yet the unique positive coexistence equilibrium has two eigenvalues in the open right half-plane. Every nonnegative equilibrium is hyperbolic and unstable; almost every strictly positive initial condition fails to converge to any equilibrium. This refutes the specified deterministic convergence conjecture, not the source's stochastic invasion results.

Varying only the dormancy probability near q = 0.97691614864474158 yields a certified supercritical Hopf bifurcation and locally asymptotically stable positive periodic solutions on its specified side. No explicit parameter radius or global periodic attractor is certified. In particular, the periodic-orbit conclusion is not asserted at q = 0.98.

## Read and replay

- `manuscript.md` and `paper.pdf`: complete statement, proofs and references (added by the publication build).
- [AI_INDEX.md](AI_INDEX.md): claim-to-evidence navigation.
- [ASSURANCE.md](ASSURANCE.md): exact scope and limitations.
- [REVIEW_RESPONSE.md](REVIEW_RESPONSE.md): substantive review dispositions.
- [PROVENANCE.md](PROVENANCE.md) and [LICENSES.md](LICENSES.md): origin and reuse.

Install Python 3.10+, Node.js 18+, and the versions in `requirements.txt`, then run:

```sh
python3 -m pip install -r requirements.txt
python3 -B run_checks.py
python3 -B -O run_checks.py
python3 -B negative_controls.py
python3 -B -O negative_controls.py
```

The main runner executes seven programs and regenerates the three certificate JSON files plus `reproduction/` logs. Run in a disposable copy when preserving deposit bytes. Python subprocesses inherit optimization and disable bytecode writes. All release predicates use explicit exceptions rather than removable `assert` statements. Semantic controls mutate isolated temporary copies only.

Exact arithmetic and interval replay support the stated finite algebraic and enclosure claims. The all-parameter inference, dynamical theorem application and measure-zero argument also require the written proof; successful execution is not a formal proof of the complete paper.
