# Internal replay receipt

Executed 18 October 2026 in the selected public derivative, using Python 3.14.7, SymPy 1.14.0, mpmath 1.3.0 and Node v25.2.1, as retained in the machine run records.

- `python3 -B run_checks.py`: PASS, seven programs; `reproduction-normal/run_record.json` retains commands, source hashes, runtimes and certificate hashes.
- `python3 -B -O run_checks.py`: PASS, seven programs; `reproduction-optimized/run_record.json` confirms optimized Python children and bytecode suppression.
- `python3 -B negative_controls.py`: PASS; `negative-normal.json`.
- `python3 -B -O negative_controls.py`: PASS; `negative-optimized.json`.

Checks include seven exact C interpolation points, the all-C Hurwitz sign pattern ++++--, thirteen q samples, the degree-17 polynomial, 99 rational arithmetic tests and 60-decimal-grid interval enclosures. Samples do not replace the all-C symbolic certificate. Every corruption control exited nonzero at the expected explicit mathematical check.

A development attempt incorrectly called the interval constructor with two arguments; the control harness rejected that unrelated TypeError, as intended, rather than counting it as a semantic pass. The corrected control calls `I.raw(2,1)` and fails on the endpoint-order predicate. Final retained receipts describe corrected controls only.

These records demonstrate local execution, not external independence or formal proof. A fresh extraction replay and final manifest/PDF checks are additional publication gates.
