# AIdmix methodology

## Model

Let `a_k` be the ancestry proportion for population `k`, constrained by

\[
a_k \ge 0, \qquad \sum_k a_k=1.
\]

At marker `j`, let `f_jk` be the population ALT frequency. The ancestry mixture
implies sample ALT frequency

\[
q_j=\sum_k a_k f_{jk}=f_j^T a.
\]

Assuming diploid Hardy-Weinberg equilibrium,

\[
P(G_j=0)=(1-q_j)^2,\quad
P(G_j=1)=2q_j(1-q_j),\quad
P(G_j=2)=q_j^2.
\]

## Read and genotype likelihoods

For an observed read base with Phred quality `Q`, define `e=10^{-Q/10}`. The
implementation bounds `e` to `[0.005, 0.25]`. For genotype dosage `g`, let
`p_g=g/2`. Then

\[
P(b=REF\mid g,Q)=(1-p_g)(1-e)+p_g e/3,
\]

\[
P(b=ALT\mid g,Q)=p_g(1-e)+(1-p_g)e/3.
\]

Log probabilities are summed over retained reads at the marker:

\[
\log L_j(g)=\sum_i \log P(b_{ij}\mid g,Q_{ij}).
\]

These are genotype likelihoods, not hard genotype calls or genotype posterior
probabilities.

## Ancestry likelihood

Marginalize the latent genotype:

\[
M_j(q_j)=\sum_{g=0}^2 L_j(g)P(G_j=g\mid q_j).
\]

The weighted log likelihood is

\[
\ell(a)=\sum_j w_j\log M_j(q_j).
\]

Background markers normally have weight 1. Target markers can be downweighted
in a joint adaptive-sampling sensitivity estimate. The estimator maximizes
`ell(a)` jointly over ancestry proportions; it does not assign each SNP or read
to a population independently.

## Closed-form gradient

For the three HWE priors,

\[
p'_0=-2(1-q),\quad p'_1=2-4q,\quad p'_2=2q.
\]

Thus

\[
M'_j(q)=\sum_g L_j(g)p'_g(q),
\qquad
\frac{\partial\ell_j}{\partial q_j}=\frac{M'_j(q_j)}{M_j(q_j)}.
\]

Because `q_j=f_j^T a`,

\[
\frac{\partial\ell}{\partial a_k}
=\sum_j w_j f_{jk}\frac{\partial\ell_j}{\partial q_j}.
\]

For `N` markers and `K` populations, the implementation computes this as

```python
q = frequencies @ ancestry
dlog_dq = ...
gradient = frequencies.T @ (weights * dlog_dq)
```

The optimizer sees only `K` parameters; all marker contributions are summed in
vectorized matrix operations.

## Closed-form Hessian

The HWE-prior second derivatives are

\[
p''_0=2,\quad p''_1=-4,\quad p''_2=2.
\]

Therefore

\[
M''_j(q)=\sum_g L_j(g)p''_g(q),
\]

\[
\frac{\partial^2\ell_j}{\partial q_j^2}
=\frac{M''_j(q_j)}{M_j(q_j)}
-\left(\frac{M'_j(q_j)}{M_j(q_j)}\right)^2,
\]

and

\[
H(a)=\sum_j w_j
\frac{\partial^2\ell_j}{\partial q_j^2}f_jf_j^T.
\]

The exact Hessian is only `K x K`. The current implementation supplies the
closed-form gradient and lets quasi-Newton optimizers approximate curvature;
the formula remains useful for diagnostics or a future Newton/trust-region
implementation.

## Optimization and validation

The robust fitter uses bounded logit-space L-BFGS-B as its primary exploration
and retains constrained SLSQP candidates as a fallback or diagnostic. Starts
include a uniform mixture and near-population vertices. Candidate solutions are
re-evaluated in ancestry space and checked with scale-normalized simplex KKT
conditions before the highest valid likelihood is selected.

Near-zero logit-space components approximate exact simplex boundaries. KKT
validation is more informative than the optimizer's status flag alone because
constrained optimizers can report false failures at valid boundary optima.

## Panel orientation

All frequencies must describe the panel ALT allele. When target-reference REF
equals panel ALT, swap REF/ALT and replace every population frequency `f` with
`1-f`. If neither panel allele equals the target-reference base, exclude the
marker. Coordinate liftover alone is insufficient.

## Interpretation limits

The estimates are coordinates relative to the selected reference panel, not
race, ethnicity, or universally defined biological populations. Results from
different marker sets or population definitions are not automatically
comparable. Sparse components can have unstable correlation while retaining
small absolute error. Convergence does not compensate for inadequate coverage,
linkage, reference-panel misspecification, or platform-specific quality
miscalibration.
