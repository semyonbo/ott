# OTT-Fixed

MATLAB OTT with two convention/coding defects repaired. Diff against the upstream
checkout at `~/git/StandingWave/ott` is **two files**.

## FIX 1 — Condon–Shortley phase (`+ott/+utils/spharm.m`)

As shipped, `spharm` returns `(-1)^m` times the standard harmonic. Verified
against the closed forms:

| harmonic | textbook (CS) | shipped `spharm` | ratio |
|---|---|---|---|
| Y₁⁰ | +0.373704 | +0.373704 | +1 |
| Y₁¹ | −0.205004 − 0.086674i | +0.205004 + 0.086674i | −1 |
| Y₂¹ | −0.350605 − 0.148234i | +0.350605 + 0.148234i | −1 |
| Y₂² | +0.111689 + 0.114999i | +0.111689 + 0.114999i | +1 |

The Farsund recursions inside `ott.forcetorque` were derived **with** the phase.
The mismatch is bilinear, so it cancels in `F_z`, `T_z` (which pair `m` with `m`)
and **survives as a sign in `F_x, F_y, T_x, T_y`** (which pair `m` with `m±1`)
whenever the coefficients span several `m`.

The derivatives `Ytheta`/`Yphi` are built by recursions linking `Y(n,m)` to
`Y(n,m±1)`; those neighbours now carry `−(−1)^m`, so each acquires one extra
sign, undone explicitly. Checked against analytic `dY/dθ` and `(1/sinθ)dY/dφ`.

## FIX 2 — spheroid T-matrix (`+ott/TmatrixSmarties.m`)

`getTmatrixData` unfolded SMARTIES' `|m|` blocks onto `±m` itself, with two
defects: it dropped the `(−1)^(s+s')` required by the achirality of a body of
revolution, and it indexed through `meshgrid(rows, cols)` — whose first argument
varies along the **columns** — storing every block transposed.

Measured: achirality violated by `2.0`, unitarity `1.06e−3` for a *lossless*
particle, ~31% median / 130% worst error vs COMSOL.

Replaced by SMARTIES' own `sparseTmatrix`, which applies the sign and returns the
matrix already in Nieminen (= `ott.utils.combined_index`) ordering.
**Requires SMARTIES on the MATLAB path.**

**Regression:** `sparseTmatrix` assembles the scattered-field T-matrix only, so
`TmatrixSmarties.simple(..., 'internal', true)` now errors instead of returning a
matrix. The code it replaced did return one, but transposed and wrong-signed in
the same way as the scattered case, so nothing that worked has been lost — only
something that appeared to. Restoring it means unfolding SMARTIES' `st4MR` blocks
with the sign and index order applied correctly.

## What is NOT changed, and why

**`wigner_rotation_matrix` is left alone.** It is exactly unitary (`|D†D−I| =
2e−15`) and equals the standard Condon–Shortley Wigner matrix Hermitian-
conjugated — verified against `treams` to **5.8e−16**:

```
|D_ott − D_treams† | = 5.796e-16
|D_ott − D_treams  | = 1.225e+00
```

So `rotate(R)` applies **`R⁻¹`**. That is a consistent convention (it is the
row-vector derivation the file's own header describes), not an error. Removing
the `D1 = D.'` transpose was tested and **breaks composition** — do not.

Before FIX 1 the phase and this convention combined to make `rotate(R)` apply
`Z·Rᵀ·Z` with `Z = diag(−1,−1,1)`, because `diag((−1)^m)` *is* the Wigner matrix
of a 180° z-rotation. **FIX 1 removes that mirror as a side effect**, leaving the
plain inverse convention.

**`translateXyz` is repaired as a side effect — measured, not assumed.**
`translateRtp`'s general branch is *rotate → translateZ → rotate back*, which
translates along `S⁻¹(ẑ)`. `translateZ` is exact; the rotation was the mirrored
one, so the composition moved the beam by `Z·d` instead of `d`. That produced the
notorious inconsistency — z-displacements matching one sign, x/y the other, and a
general displacement matching neither. Removing the mirror fixes it:

| d | `E(r−d)` | `E(r+d)` |
|---|---|---|
| (0, 0, +0.20) | 1.7e+00 | **1.0e−06** |
| (+0.20, 0, 0) | 6.3e−01 | **9.5e−07** |
| (0, +0.20, 0) | 1.2e+00 | **1.3e−06** |
| (+0.13, −0.09, +0.21) | 1.7e+00 | **1.2e−06** |
| (−0.05, +0.17, −0.12) | 3.3e−01 | **1.2e−06** |

All displacements now obey one rule, `translateXyz(d) → E(r + d)`.
(1e−06 is the finite-`Nmax` floor.)

**`BscPlane(θ,φ)` still propagates along `−r̂(θ,φ)`.** Independent convention,
unaffected by the phase — and note it disagrees with `BscPmGauss`, which
propagates along `+z`. Not a defect (`|F|` and `F·k̂` are both right), but a trap
if you assume the two beam classes share a sense.

## Consequences

| unchanged | changed (these were wrong) |
|---|---|
| on-axis beams, z-translations | oblique `BscPlane` |
| cross sections, all fields | displaced/rotated Gaussians (140–190%) |
| `axial_equilibrium`, `find_traps` | tilted non-axisymmetric particles |
| `F_z`, `T_z`, `\|F\|` anywhere | `F_x, F_y, T_x, T_y` with mixed `m` |

**Any script that compensated for the old rotation convention must be updated** —
notably `swOrientation.m` in the StandingWave project, whose transposed `Rz`/`Ry`
cancelled the old mirror.

---

# Appendix A — the root cause, at the Legendre level

`spharm` never adds the Condon–Shortley phase, and neither does the Legendre
routine underneath it. In the usual framing — *"the phase may be included either
in the associated Legendre polynomials or in the spherical harmonics, but it
should not be counted twice"* — **OTT counts it zero times.**

`+ott/+utils/legendrerow.m` does not call MATLAB's `legendre`. It rolls its own
recursion, citing **Holmes & Featherstone (2002)** and **Jekeli (2007)** — geodesy
papers. Geodesy conventionally omits Condon–Shortley, which is a quantum-mechanics
convention. Every term is built from positive square roots:

```matlab
Wnn = sqrt((2*n+1)/(4*pi)*prod(1-1/2./[1:n]))      % positive
pnm = a*ct.*pnm(jj+1,:) - b*st.^2.*pnm(jj+2,:);    % a, b > 0
pnm = pnm.*ST.^(M);                                % sin^m θ > 0
```

Measured signs at `θ = 0.7` (so `sinθ > 0`):

| n | m | `legendrerow` | `legendre` (has CS) | `legendre(...,'norm')` |
|---|---|---|---|---|
| 1 | 1 | +0.222573 | **−0.644218** | +0.557909 |
| 2 | 1 | +0.380654 | **−1.478175** | +0.954158 |
| 2 | 2 | +0.160310 | +1.245049 | +0.401838 |
| 3 | 3 | +0.111549 | **−4.010414** | +0.279613 |

`legendrerow` follows the CS-**free** pattern. `spharm` then applies `(−1)^m` only
to the `m < 0` rows — that is the `Y_{n,−m} = (−1)^m Y_{n,m}*` bookkeeping, not the
Condon–Shortley phase.

Nothing is wrong with that choice *in isolation*. The bug is that
`ott.forcetorque` implements Farsund's recursions, derived in the
quantum-mechanics convention. **A geodesy Legendre routine feeding a
quantum-mechanics force formula** — two communities' conventions meeting at a
boundary nobody tested.

### Trap when re-deriving this

MATLAB's `'norm'` option carries its own `(−1)^m`:

```
N_n^m(x) = (−1)^m · sqrt((n+½)(n−m)!/(n+m)!) · P_n^m(x)
```

which **cancels** the one inside `P_n^m`. So `legendre(...,'norm')` is CS-*free*
while plain `legendre(...)` is CS-*carrying*. Comparing `legendrerow` against
`'norm'` makes the two look identical and hides the phase completely. Compare
against the unnormalised form, or against the closed-form harmonics.

---

# Appendix B — where the `(−1)^m` comes from at all

**Nothing forces it.** It is a phase convention on the basis states. It is *not*
required by:

* the Legendre or Helmholtz differential equation — solutions are fixed only up
  to a multiplicative constant;
* normalisation — a unit-modulus phase does not change `∫|Y|² dΩ = 1`;
* orthogonality — a diagonal phase preserves it;
* the spherical Bessel functions — those are the *radial* part; this lives
  entirely in the angular part.

It is chosen so the **angular-momentum ladder operators have real, positive
matrix elements**:

```
L± |l,m⟩ = sqrt( l(l+1) − m(m±1) ) |l,m±1⟩ = sqrt( (l∓m)(l±m+1) ) |l,m±1⟩
```

Condon & Shortley (*The Theory of Atomic Spectra*, 1935) fixed the relative
phases of the `|l,m⟩` states precisely so this holds with no stray signs.
Everything built on angular-momentum algebra inherits that choice: Wigner 3j/6j
symbols, the Wigner–Eckart theorem, the Wigner d-matrices, and the standard
tables. In the formula the phase literally sits here:

```
P_l^m(x) = (−1)^m (1 − x²)^(m/2) dᵐ/dxᵐ P_l(x)
```

and the convention lets you move it to `Y` instead — but exactly once.

## Why this is fatal specifically for force, and harmless elsewhere

**Farsund's sums are angular-momentum ladder algebra.** The coefficient that
appears throughout them,

```
Λ_lm = sqrt( (l−m)(l+m+1) )
```

*is* the `L₊` matrix element above; the `l → l±1` coefficients are
Clebsch–Gordan factors for a rank-1 (vector) operator. Those numbers are only
correct in a basis that carries the CS phases.

Feed them a non-CS basis and every term that **moves `m` by ±1** picks up a
relative sign,

```
(−1)^m · (−1)^(m±1) = −1
```

while every term that **keeps `m` fixed** is immune,

```
(−1)^m · (−1)^m = +1
```

By the selection rules, `F_z` and `T_z` pair `m` with `m`, and
`F_x, F_y, T_x, T_y` pair `m` with `m±1`. That is exactly the observed pattern:
**axial components correct, transverse components sign-flipped** — and only when
the coefficients actually span several `m`, which is why on-axis beams on
axisymmetric particles never showed it.
