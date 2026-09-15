# Verification still to do

`FIXES.md` reports what was measured. This is what is **not** yet nailed down, and
exactly how to close each gap. The honest state of the evidence:

| claim | status |
|---|---|
| FIX 2 (T-matrix) is right | **done** — COMSOL, 114 orientations, 53.6% → 0.85% on `F_z` |
| FIX 2 reproduces SMARTIES | **done** — `3.5e-16`, achirality `0.0`, `verify_tmat.m` |
| FIX 1 leaves `F_z`, `T_z` alone | **done** — sphere, two builds bit-identical (`0.000e+00`) |
| `spharm` was missing `(−1)^m` | **done analytically** — closed forms, §1.1 |
| **FIX 1 gives the right transverse force/torque** | **NOT DONE** — see §1.3, §2 |

The gap matters: the existing COMSOL dataset cannot test FIX 1 at all. It uses a
sphere (no transverse force, no torque in that beam) and a spheroid (dominated by
FIX 2). Nothing in it exercises an `m ↔ m±1` observable. **FIX 1 is currently
justified by theory plus agreement with `swforce`, not by independent numerics.**

---

# 1. Analytic checks — no FEM needed

## 1.1 The phase itself (done, recorded here for completeness)

Condon & Shortley fix the relative phases of `|l,m⟩` so the ladder operators have
real, positive matrix elements:

```
L± |l,m⟩ = sqrt( (l∓m)(l±m+1) ) |l,m±1⟩
```

Nothing in the Helmholtz equation, the normalisation, or the orthogonality
requires it — a unit-modulus phase preserves all three. It is a basis convention,
and it must be applied **exactly once**: conventionally either in

```
P_l^m(x) = (−1)^m (1 − x²)^(m/2) dᵐ/dxᵐ P_l(x)
```

or in `Y_l^m`, never both. OTT applied it zero times (Appendix A of `FIXES.md`).

Check against closed forms, at any `(θ, φ)` with `sinθ > 0`:

```
Y_1^0 = sqrt(3/4π) cosθ
Y_1^1 = −sqrt(3/8π) sinθ e^{iφ}          <- the minus sign IS the phase
Y_2^1 = −sqrt(15/8π) sinθ cosθ e^{iφ}
Y_2^2 = +sqrt(15/32π) sin²θ e^{2iφ}
```

Odd `m` must come out negative relative to the CS-free form; even `m` unchanged.
**Trap:** MATLAB's `legendre(...,'norm')` carries its own `(−1)^m` that cancels
the one inside `P_l^m`, so comparing against `'norm'` hides the phase entirely.
Compare against the unnormalised `legendre(...)` or against the closed forms.

## 1.2 Ladder-operator test (to do)

Build the matrix of `L₊` in OTT's own angular basis by projecting
`L₊ Y_l^m` onto `Y_l^{m+1}` using OTT's `spharm` on a quadrature grid. Every
element must be **real and positive**, equal to `sqrt((l−m)(l+m+1))`. Before the
fix, alternate elements come out negative. This tests the phase in precisely the
algebra `ott.forcetorque` depends on, and needs nothing but `spharm`.

## 1.3 Force direction on a sphere (to do — the key analytic test of FIX 1)

For a **sphere** in a single plane wave, symmetry forces

```
F = khat * I * C_pr / c ,     C_pr = C_ext − g C_sca
```

exactly parallel to `khat`, with `C_pr` available in closed form from Mie theory
(Bohren & Huffman §4.6). At oblique incidence this has large transverse
components, so it is a direct test of `F_x, F_y` — and a sphere never touches
`TmatrixSmarties`, so it isolates FIX 1.

**Beware the circularity.** Building the beam with `BscPlane(θ,φ)` and then
checking `F ∥ khat` proves nothing on its own: the phase corrupts the beam
expansion too, so the force can come out parallel to a *mirrored* `khat` and look
self-consistent. The beam must first be anchored outside OTT:

1. evaluate `beam.emFieldXyz` at off-axis points (it is singular on the z axis);
2. compare against the analytic plane wave `E = ê E₀ exp(i k·r)` computed directly;
3. only once the field matches, test `F ∥ khat` and `|F|` against Mie `C_pr`.

Acceptance: `F·khat/|F| = 1` to truncation, and `|F|` matching Mie `C_pr` to the
same. Expected to fail on the unpatched build at oblique incidence and pass on the
patched one.

## 1.4 Torque anchor (to do)

A lossless sphere has zero torque, so it cannot test `T`. Use an **absorbing**
sphere (complex `index_particle`) in **circularly polarised** light, where
`T_z = C_abs I / ω` in closed form. That anchors `T_z` absolutely — though `T_z`
is `m`-diagonal and therefore immune to FIX 1, so it validates the torque
normalisation, not the phase. Transverse torque needs §2.2.

---

# 2. COMSOL verification — what to build

Existing model `force_torque_vs_angles.mph` covers FIX 2. These are the two new
systems needed for FIX 1. Both must avoid `TmatrixSmarties` so the two fixes stay
separable — i.e. **no spheroids**.

## 2.1 System A — sphere at oblique incidence (force only)

The minimal FIX 1 test, and the same physical case as §1.3, so the analytic and
FEM answers cross-check each other.

```
particle    sphere, r = 0.35 λ₀,  n_p = 4.0  (lossless),  medium vacuum
beam        single travelling plane wave, linear polarisation
angles      θ = 0° (control), 30°, 45°, 60°;  φ = 0°, 30°, 63°
```

`θ = 0` is the control: every build agrees there, which proves the COMSOL model is
set up correctly before the oblique numbers are believed. The discriminating
quantity is the **direction** of `F`, not its magnitude — the phase defect flips
`F_x, F_y`, so unpatched OTT should give a force that is not parallel to `k̂`,
while COMSOL and patched OTT both give `F ∥ k̂`.

`swforce_compare/matlab/comsol_validate.m` in the BIC-Force project already
implements exactly this and has never been run to completion — its output
`comsol_validate.json` does not exist. **Running it is the single highest-value
remaining task.** It needs `comsol mphserver -port 2036`.

## 2.2 System B — tilted cylinder (force *and* transverse torque)

A sphere gives no torque, so transverse torque needs a non-axisymmetric particle
that does not go through SMARTIES. OTT provides `TmatrixPm` (point-matching,
cylinders and other axisymmetric shapes) and `TmatrixDda`.

```
particle    cylinder, radius 0.3 λ₀, height 0.8 λ₀, n_p = 1.6 (lossless)
            built with ott.TmatrixPm  -- NOT TmatrixSmarties
beam        single plane wave, linear polarisation along x
tilt        cylinder axis tilted β = 0°, 15°, 30°, 45°, 60°, 90° from z,
            azimuth α = 0°, 45°, 90°
```

β = 0 and β = 90 with α = 0 are symmetry controls where transverse torque must
vanish identically. The tilted cases give `T_x, T_y` of order `|F| · h`, which is
where FIX 1 lives.

Rotate the **particle** in COMSOL and the **beam** in OTT, or vice versa, but
state which — the rotation convention is itself one of the things being tested
(`FIXES.md`, "What else changes"). Safest is to fix the particle in COMSOL and
sweep the beam direction, matching `vars_background_field.txt` from the existing
model.

## 2.3 COMSOL settings (same recipe as the existing model)

```
physics        Wave Optics (ewfd), scattered-field formulation
domain         particle + vacuum shell + spherical PML
force          intop over a sphere of radius R_int enclosing the particle:
                 F_i = ∫ ewfd.unTi dS
torque         T_x = ∫ (y·unTz − z·unTy) dS, cyclic
mesh           sweep λ/4, λ/6, λ/8 in the vacuum and confirm the answer
               has stopped moving -- report the converged value, not one mesh
R_int          sweep it too: in a lossless medium the result must be
               independent of the integration radius.  This is the cheapest
               and strongest check that the FEM model itself is sound.
```

## 2.4 Acceptance criteria

| quantity | expected |
|---|---|
| `θ = 0` control, all builds | agree with COMSOL to mesh convergence (~1%) |
| patched OTT, all cases | within a few % of COMSOL, force **and** torque |
| unpatched OTT, oblique/tilted | transverse components wrong; direction cosine < 1 |
| unpatched OTT, `F_z`, `T_z` | still correct — if these move, FIX 1 is not the whole story |

The last row is the falsifier. The selection-rule argument predicts the axial
components are untouched by the phase; if COMSOL shows unpatched `F_z` wrong on a
**sphere**, the explanation in `FIXES.md` Appendix B is incomplete.

---

# 3. Before submitting upstream

- [ ] §1.2 ladder-operator test
- [ ] §1.3 sphere force direction vs Mie `C_pr`, both builds
- [ ] §2.1 run `comsol_validate.m`, save `comsol_validate.json`
- [ ] §2.2 cylinder force + torque vs COMSOL, both builds
- [ ] restore the internal-field path (`st4MR`) or document it as out of scope —
      currently `TmatrixSmarties.simple(..., 'internal', true)` errors
- [ ] check `TmatrixPm` and `TmatrixDda` for the same `meshgrid` transpose bug
      that FIX 2 removed from `TmatrixSmarties`; the pattern may repeat
- [ ] run OTT's own test suite (`tests/`) on the patched build and report diffs,
      since some tests may encode the old convention
