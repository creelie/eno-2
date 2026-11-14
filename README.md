# Minimal Sand in Manna's Sandpile and Exact Relaxation of the Oslo Rice Pile

Deep Bhattacharjee

Programs that check the finite statements of the paper. No proof in the paper depends on them.

- `lean/Sandpile.lean`: Lean 4 certificate (core library only): the recursion of G, the same-parity families, the legal histories of the chain and of the Gray cube, the matrix identities of the Oslo model, the count of recurrent configurations and the relaxation times on small piles. Run `lean lean/Sandpile.lean` (Lean 4.15).
- `julia/manna.jl`: exhaustive recurrent sets on small boxes and cylinders, the constructions executed with every toppling checked, the families and the recursion of G. Run `julia julia/manna.jl`.
- `julia/oslo.jl`: exact rational laws of the Oslo model for L ≤ 6, the steady state, the recurrent set and every relaxation time. Run `julia julia/oslo.jl`.
- `python/construct.py`: the layer construction, executed with legality checks.
- `python/ilp_lower.py`: the integer programme of the appendix (numpy, scipy).
- `python/oslo_support.py`: the support of the Oslo dynamics for L ≤ 8.
- `python/make_figures.py`: draws the figures (matplotlib).

License: MIT.
