# Assurance and limitations

The scientific output combines a written proof, exact rational/symbolic computations and fixed-grid outward-rounded interval arithmetic. Re-execution is internal reproducibility checking, not independent reproduction or external peer review. The interval implementation is itself part of the trusted code base; arithmetic self-tests and semantic controls reduce, but do not eliminate, implementation risk.

The all-C family is proved using coefficient signs, not inferred from samples. The Hopf proof includes positivity, root isolation, noncritical spectral stability, transversality and a negative first Lyapunov coefficient in the declared normalization. It gives only a sufficiently small one-sided neighborhood. No numeric size of that neighborhood is provided. Nonconvergence to equilibria does not characterize all limit sets.

The local source review supports the precise conjecture comparison; novelty is a literature assessment, not a machine-certified property. Source-paper statements and hypotheses should be consulted when transferring the result to another model.
