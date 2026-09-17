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

Measured on an oblate spheroid (`ellipsoid [0.25 0.25 0.12]`, `n_rel = 2.5`):
achirality violated by `2.0` — a 100% breach — and unitarity `1.06e−3` for a
*lossless* particle. After the fix, reproducible with
`swforce_compare/matlab/verify_tmat.m` in the BIC-Force project:

```
max |patched - direct SMARTIES| : 3.511e-16
achirality T12(-m) = -T12(+m)   : 0.000e+00
unitarity                       : 7.091e-05   (Nmax=6 truncation, not a defect)
```

Replaced by SMARTIES' own `sparseTmatrix`, which applies the sign and returns the
matrix already in Nieminen (= `ott.utils.combined_index`) ordering.
**Requires SMARTIES on the MATLAB path.**

**Independent confirmation.** `ott.TmatrixEbcm` solves the same spheroid by a
completely different method (null-field / extended boundary condition). Against
it, for `ellipsoid [0.25 0.25 0.12]`, `n_rel = 2.5`:

| | vs EBCM | unitarity (lossless) |
|---|---|---|
| unpatched | 9.98e−02 | 4.03e−03 |
| **patched** | **3.82e−03** | **6.95e−04** |

26x better agreement and 5.8x better unitarity. Note EBCM must be run at
`npts ≈ 50`: at `npts ≥ 100` its own system matrix goes singular
(`RCOND = 3e−27`) and it returns garbage, which is easy to mistake for a
disagreement.

**The internal-field option** (`'internal', true`) routes SMARTIES' `st4MR`
blocks through the same `sparseTmatrix` assembly. Three things are established
about it, and one is not:

* **Definition and normalisation — exact.** A spheroid with `a_x = a_z` *is* a
  sphere, and the internal Mie coefficients are analytic. Against
  `ott.TmatrixMie(..., 'internal', true)` at matched `Nmax`, the SMARTIES
  internal matrix agrees to **7e−16**, and the scattered one to **3.6e−16**.
  Deviation then grows smoothly with aspect ratio (1.7e−2 at 1.02, 8.4e−2 at
  1.10). So what OTT calls an internal T-matrix and what SMARTIES calls
  `R = Q⁻¹` are the same object, in the same normalisation.
* **Assembly — cross-checked.** Repairing OTT's original hand-rolled unfold
  (achirality sign restored, `meshgrid` → `ndgrid`) and comparing against the
  `sparseTmatrix` route gives `0.0000e+00` for the internal matrix and for the
  scattered one.
* **Not validated: a genuinely non-spherical internal matrix.** The sphere limit
  is *blind* to both defects FIX 2 repairs — for a sphere the T-matrix is
  diagonal, so a transposed block is unchanged and the `M12`/`M21` blocks
  carrying the `m < 0` sign are zero. Both builds pass it identically.

  And the obvious physical test is **unavailable**, for a reason worth recording.
  Tangential `E` must be continuous across the surface, but neither expansion
  converges there: the outgoing (scattered) series is valid only outside the
  circumscribing sphere `r > a`, the regular (internal) series only inside the
  inscribed sphere `r < c`, and a spheroid's surface spans both. Measured for
  `a = 0.25, c = 0.12`, the scattered field grows from `0.13` at `r = 1.0` to
  `5.6e+09` at `r = 0.06` — the Rayleigh hypothesis, not a defect. A boundary
  -condition check works only for a sphere, where the two radii coincide.

  That leaves `TmatrixEbcm`'s internal matrix as the only non-spherical
  reference, disagreeing by 16%; EBCM is the weaker method here (scattered
  unitarity 5.5e−3 against SMARTIES' 6.9e−4), so that figure more likely bounds
  EBCM than SMARTIES. Settling it needs an independent code — SMUTHI, DDA, or
  FEM internal fields. **Unresolved.**

## Validation against COMSOL

COMSOL model `force_torque_vs_angles.mph`: standing wave at `lam0 = 1550 nm`,
phase `ph_s = 45 deg`, `eps_p = 12`, with the wave axis and polarisation swept
over the particle. Force and torque from COMSOL's own Maxwell stress tensor on a
surrounding sphere. Two particles, so the two fixes can be told apart — a
**sphere** (r = 250 nm) uses analytic Mie and never touches `TmatrixSmarties`,
while the **spheroid** (357/357/250 nm) goes through it.

Median relative error against COMSOL:

| | `F_z` | `F_xy` | `T_z` | `T_xy` |
|---|---|---|---|---|
| sphere, unpatched | 0.49% | 0.58% | — | — |
| sphere, patched | 0.49% | 0.58% | — | — |
| spheroid, unpatched | 53.6% | 56.7% | — | 63.7% |
| **spheroid, patched** | **0.85%** | **0.49%** | — | **0.32%** |

`—` marks a component COMSOL itself returns as numerically zero, where no
relative error is defined: a lossless sphere in this beam feels no transverse
force and no torque at all, and `T_z` vanishes by symmetry for both particles.

**The spheroid discrepancy is FIX 2, not FIX 1.** On the sphere the two builds
are bit-identical — `max |unpatched − patched| = 0.000e+00` in both force and
torque — so the Condon–Shortley phase contributes nothing to the spheroid error.
Confirmed directly on a single configuration:

```
             sphere (analytic Mie)   spheroid (SMARTIES)
unpatched    Fz = -84.9768490        Fz =  -44.275
patched      Fz = -84.9768490        Fz = -111.056
```

`F_z` is identical to every digit on the sphere, exactly as the selection-rule
argument in Appendix B requires, and off by a factor 2.5 on the spheroid.

**This dataset does not exercise FIX 1 — at all.** Swapping the T-matrix for
`ott.TmatrixEbcm`, which bypasses FIX 2, and running it *unpatched* gives
0.77% / 0.54% / 0.42% — bit-identical to patched (`max |ΔF| = 0.000e+00` over all
114 orientations). The beam here is built from `BscPlane(0,0)` and
`BscPlane(pi,0)`, both on axis, then rotated; `diag((−1)^m)` is the Wigner matrix
of a 180° z-rotation, so on an on-axis beam the missing phase is a global phase
and cancels in the bilinear force. **The entire 53% error above is FIX 2.**

### FIX 1, validated separately — sphere, beam built at an angle

A sphere in a single plane wave must feel a force exactly along `k̂`. That is an
absolute symmetry statement needing no FEM, and a sphere uses `TmatrixMie`, so
FIX 2 cannot interfere. The beam must be *constructed* at the angle, not built
on axis and rotated.

| θ | unpatched `F·k̂/\|F\|` | patched |
|---|---|---|
| 0° | 1.00000 | 1.00000000 |
| 30° | 0.50000 | **1.00000000** |
| 45° | 0.00000 | **1.00000000** |
| 60° | −0.50000 | **1.00000000** |
| 75° | −0.86603 | **1.00000000** |

Unpatched, the force is perpendicular to the beam at 45° and points backwards at
75°. The measured ratio is `cos(2θ)` at every angle — the signature of
`F = diag(−1,−1,1)·F_true`. The transverse components flip sign while `F_z` is
bit-identical between builds (55.96, 45.69, 32.31, 16.72), which is exactly the
selection rule of Appendix B. Reproduce with
`swforce_compare/matlab/probe_fix1_oblique.m`.

Reproduce with `swforce_compare/matlab/ott_vs_comsol.m` (once per OTT build, per
shape) and `swforce_compare/verification/report_ott_vs_comsol.py`.

## What else changes, even though these files don't

Removing the phase changes the behaviour of code it feeds. Both items below are
things to check in existing scripts, not defences of code nobody questioned.

**Rotations.** `diag((−1)^m)` *is* the Wigner matrix of a 180° z-rotation, so
while the phase was missing, `rotate(R)` applied `Z·Rᵀ·Z` with `Z = diag(−1,−1,1)`
— a mirrored rotation. It now applies plain `Rᵀ`, i.e. `R⁻¹`, which is the
row-vector convention `wigner_rotation_matrix`'s own header describes.

`wigner_rotation_matrix` itself is deliberately **not** patched: it is unitary to
`2e−15` and matches `treams`' Wigner matrix Hermitian-conjugated to `5.8e−16`, so
the inverse convention is a choice, not a defect. Recorded here only because
removing its `D1 = D.'` transpose looks like the obvious fix and **breaks
composition** — don't.

**Translations.** `translateXyz` is fixed as a consequence. Its general branch is
*rotate → translateZ → rotate back*; `translateZ` is exact, but the rotation was
the mirrored one, so the composition moved the beam by `Z·d` instead of `d`. That
is the long-standing inconsistency where z-displacements matched one sign and x/y
the other:

Measured on a focused Gaussian, relative error against each candidate:

| d (units of λ) | unpatched `E(r−d)` | unpatched `E(r+d)` | patched `E(r+d)` |
|---|---|---|---|
| (0, 0, +0.20) | 1.8e+00 | **6.7e−15** | **6.7e−15** |
| (+0.20, 0, 0) | **5.2e−15** | 7.3e−01 | **5.3e−15** |
| (0, +0.20, 0) | **6.0e−15** | 7.8e−01 | **4.9e−15** |
| (+0.13, −0.09, +0.21) | 1.8e+00 | 5.8e−01 | **6.7e−15** |
| (−0.05, +0.17, −0.12) | 1.2e+00 | 6.9e−01 | **5.7e−15** |

Unpatched, z-displacements obeyed `E(r+d)`, x and y obeyed `E(r−d)` — the
opposite sign — and a general displacement matched **neither**. Patched, all five
obey `translateXyz(d) → E(r + d)` to machine precision. Reproduce with
`swforce_compare/matlab/probe_translate.m`.

**`BscPlane(θ,φ)` propagates along `−r̂(θ,φ)`**, unchanged and unrelated to the
phase, but note it disagrees with `BscPmGauss`, which goes along `+z`.

## What changes in existing results

The two fixes have different blast radii and must not be conflated.

**FIX 2 changes every result for a spheroid**, including `F_z` and `|F|`, because
the T-matrix itself was wrong — measured factor 2.5 on `F_z` above. Anything
built on `ott.Tmatrix.simple('ellipsoid', ...)` or `TmatrixSmarties` has to be
recomputed. It affects nothing else: spheres go through `TmatrixMie`, and other
shapes through their own classes.

**FIX 1 is narrower**, and the selection rules say exactly how narrow:

| unchanged by FIX 1 | changed by FIX 1 |
|---|---|
| `F_z`, `T_z` — terms pairing `m` with `m` | `F_x, F_y, T_x, T_y` — terms pairing `m` with `m ± 1` |
| on-axis beams on axisymmetric particles | oblique `BscPlane` |
| cross sections | rotated particles and tilted beams |
| `axial_equilibrium`, `find_traps` | `translateXyz` off the z axis (see below) |

Verified rather than asserted: on a sphere, where FIX 2 cannot apply, the two
builds agree to `0.000e+00` in both force and torque for this beam.

**Any script that compensated for the old rotation convention must be updated.**
The StandingWave project's `sw` package is the example: `sw.Dz` corrects the
force/torque output and `sw.toParticleFrame` rotates with `Dz*R*Dz`. Applied to
both, the two cancel and the workaround is inert; applied to only one, it now
introduces the very error it used to remove.

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

# Appendix B — why only the transverse components broke

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

