# OTT-Fixed

Two defects repaired in MATLAB OTT. Functional diff: **two files**.

---

## 1. Missing Condon–Shortley phase — `+ott/+utils/spharm.m`

`spharm` returned `(−1)^m Y_l^m` instead of `Y_l^m`. The Legendre routine beneath
it, `legendrerow`, implements the Holmes & Featherstone (2002) recursion from
geodesy, where the phase is conventionally omitted; nothing downstream restored
it. Meanwhile `ott.forcetorque` implements Farsund & Felderhof's expressions,
which assume it — the ladder coefficient `sqrt((n−m)(n+m+1))` appears throughout.

Three lines restore it (the derivative recursions pick up one extra sign):

```matlab
Y      = Y .* (-1).^mv;
Ytheta = -Ytheta;
Yphi   = -Yphi;
```

**Evidence.** A sphere in a plane wave must feel a force exactly along `k̂` —
symmetry, no formula required:

| θ | unpatched `F·k̂/|F|` | patched |
|---|---|---|
| 30° | 0.50000 | 1.00000000 |
| 45° | 0.00000 | 1.00000000 |
| 60° | −0.50000 | 1.00000000 |

Unpatched gives `cos 2θ` at every angle — the signature of `F → diag(−1,−1,1)·F`.
At 45° the force is perpendicular to the beam. `F_z` is bit-identical between
builds; only `F_x, F_y` flip.

## 2. SMARTIES T-matrix import — `+ott/TmatrixSmarties.m`

`getTmatrixData` unfolded SMARTIES' `|m|` blocks onto `±m` by hand and got two
things wrong: it dropped the `(−1)^(s+s′)` that achirality of a body of
revolution requires for `m < 0`, and it indexed through `meshgrid(rows, cols)`,
whose first output varies along columns — storing every block transposed.

Replaced by SMARTIES' own `sparseTmatrix`, which applies the sign and returns
Nieminen ordering directly. The `internal` option routes the `st4MR` blocks
through the same call.

**Evidence.** Against `ott.TmatrixEbcm`, an independent null-field solver, on an
oblate spheroid: unpatched `9.98e−02`, patched `3.82e−03`. Unitarity for a
lossless particle: `4.03e−03` → `6.95e−04`.

---

## What changes for users

The phase is a basis gauge change, so most of the toolbox is invariant.

| unchanged | changed |
|---|---|
| all field evaluation, far fields, `power` | `F_x, F_y, T_x, T_y` for beams spanning several `m` |
| `F_z`, `T_z` everywhere | `translateXyz` off the z axis |
| `TmatrixMie`, `TmatrixEbcm` | every spheroid result |
| `translateZ`, `axial_equilibrium` | `rotateX`, `rotateY`, `scatter(…,'rotation')` |

`rotateZ` is unaffected. OTT's own test suite gives identical results on both
builds; its three failures are pre-existing upstream.

**`translateXyz` was internally inconsistent** and is fixed as a consequence, not
by design: `wigner_rotation_matrix` never reads `spharm`, so the basis moved
under it. Requesting `+0.3λ` unpatched moved the beam to `+0.3` along x but
`−0.3` along z — two branches of the same function, only one of which uses the
rotation path. Patched, all axes obey one rule to `5e−15`.

**The two defects cancel for the commonest calculation.** A laterally displaced
on-axis beam has the beam put on the wrong side *and* `F_x` negated, and the
product is right: the lateral force curve is identical in both builds to
`0.000e+00`. Axial and lateral trap stiffness with an on-axis beam were never
affected — which is how this survived. What was wrong is narrower: beams
**constructed** at an angle, rotated non-axisymmetric particles, and the
**position** of a laterally translated beam if you looked at its field.

Scripts that compensated for the old convention must drop the workaround.

---

## Root cause

`spharm` never adds the phase, and neither does `legendrerow`. In the usual
framing — *the phase may sit in the associated Legendre polynomials or in the
spherical harmonics, but must not be counted twice* — **OTT counts it zero
times**.

`legendrerow` cites Holmes & Featherstone (2002) and Jekeli (2007), both geodesy,
where Condon–Shortley is conventionally absent. Every term is built from positive
square roots. Measured against MATLAB's unnormalised `legendre`, the ratio is
exactly `(−1)^m` at every `(n,m,θ)`.

The `(−1)^m` that *is* present in `spharm` applies only to `m < 0` rows: that is
the `Y_{n,−m} = (−1)^m Y_{n,m}*` relation, not the phase.

**Trap when re-deriving this.** MATLAB's `legendre(…,'norm')` carries its own
`(−1)^m`, which cancels the one inside `P_n^m`. Comparing `legendrerow` against
`'norm'` makes them look identical and hides the defect. Compare against the
unnormalised form or the closed-form harmonics.

So: a geodesy Legendre routine feeding a quantum-mechanics force formula — two
conventions meeting at a boundary nobody tested.

## Why only the transverse components broke

Farsund's sums are angular-momentum ladder algebra. The coefficient throughout
them, `Λ_lm = sqrt((l−m)(l+m+1))`, *is* the `L₊` matrix element, correct only in
a basis carrying the Condon–Shortley phases.

Feed them a phase-free basis and every term that moves `m` by `±1` picks up

```
(−1)^m · (−1)^(m±1) = −1
```

while every term keeping `m` fixed is immune. By the selection rules `F_z` and
`T_z` pair `m` with `m`, and `F_x, F_y, T_x, T_y` pair `m` with `m±1`. Hence
axial components correct, transverse components negated — and only when the
coefficients span several `m`, which is why on-axis beams on axisymmetric
particles never showed it.

## Where the phase comes from

Nothing forces it. Not the differential equation (solutions are fixed only up to
a constant), not normalisation (a unit-modulus phase preserves `∫|Y|²dΩ = 1`),
not orthogonality. It is a choice, made so the ladder operators have real,
positive matrix elements:

```
L± |l,m⟩ = sqrt( (l∓m)(l±m+1) ) |l,m±1⟩
```

Everything built on angular-momentum algebra inherits it: 3j symbols, the
Wigner–Eckart theorem, the Wigner d-matrices — and Farsund's force sums.

---

## Validation

| test | unpatched | patched |
|---|---|---|
| `spharm` vs closed-form `Y_1^1` | ratio −1 | ratio +1 |
| sphere, `F·k̂/\|F\|` at 45° | 0.000 | 1.000 |
| `translateXyz`, x vs z | inconsistent | consistent, `5e−15` |
| spheroid T-matrix vs EBCM | `9.98e−02` | `3.82e−03` |
| internal T-matrix vs analytic Mie (sphere limit) | — | `6.96e−16` |
| COMSOL, force sign, 11 wavelengths | 0/11 | **11/11** |
| COMSOL, force direction (cosine) | 0.033 | **0.994** |
| COMSOL, absorbed fraction | 11.7% | **2.2%** |

`swforce_compare/matlab/CHECK_EVERYTHING.m` runs the first five in about a
minute. None of them compares one of my implementations against another.

**Not established:** the internal T-matrix for a genuinely non-spherical
particle. The sphere limit is exact but blind to both defects (a sphere's
T-matrix is diagonal). A boundary-condition check is unavailable because neither
VSWF series converges on a spheroid's surface — the outgoing series needs
`r > a`, the regular one `r < c`, and the surface spans both. `TmatrixEbcm`
disagrees by 16%, but is the weaker method here.
