# pionCT_kinematics_v2.C

ROOT macro for computing exclusive pion electroproduction kinematics at forward production.

## Reaction

```
e(k) + p(P)  →  e'(k') + π⁺(p_π) + n(P')
```

One-photon-exchange approximation. Target proton at rest in the lab. All
quantities computed at **θ\* = 0** (forward pion production in the hadronic
CM frame), which is the kinematic point of physical interest for pion-pole
dominated cross sections.

## Requirements

- ROOT (any recent version with `TLorentzVector` and `TMath`)
- No external libraries

## Usage

```bash
root -l -q "pionCT_kinematics_v2.C+(5.0, 10.7, \"W\", 2.42)"
root -l -q "pionCT_kinematics_v2.C+(5.0, 10.7, \"t\", -0.40)"
root -l -q "pionCT_kinematics_v2.C+()"
```

Or interactively:

```
root -l
.L pionCT_kinematics_v2.C+
pionCT_kinematics_v2(5.0, 10.7, "W", 2.42)
pionCT_kinematics_v2(5.0, 10.7, "t", -0.40)
pionCT_kinematics_v2()
```

### Arguments

| Argument | Type | Description |
|---|---|---|
| `Q2` | `Double_t` | Virtuality (GeV², positive). Set to 0 to print comparison table. |
| `Ebeam` | `Double_t` | Beam energy (GeV) |
| `mode` | `const char*` | `"W"` to specify W directly; `"t"` to find W from t_min; `"table"` for full table |
| `value` | `Double_t` | Value of W (GeV) or t_target (GeV², negative), depending on mode |

## Modes

### Mode A — specify W directly

```
pionCT_kinematics_v2(Q2, Ebeam, "W", W_value)
```

Computes all lab observables at the given W and Q², at θ\* = 0. Reports
`t_min` as an output.

**Example:**
```
pionCT_kinematics_v2(5.0, 10.7, "W", 2.42)
```

### Mode B — specify t_min

```
pionCT_kinematics_v2(Q2, Ebeam, "t", t_target)
```

Finds the W such that `t_min(W, Q², E) = t_target` using bisection, then
computes all observables at that W. This is the mode that reproduces the
reference kinematic table.

**Example:**
```
pionCT_kinematics_v2(5.0, 10.7, "t", -0.40)
```

### Mode C — comparison table

```
pionCT_kinematics_v2()
```

Prints the full computed table for Q² ∈ {5.0, 6.5, 7.5, 8.5} GeV² at
E = 10.7 GeV, t_min = −0.40 GeV², side by side with the paper reference
values.

## Output

Each single-point call prints:

```
=== Kinematics (mode: W) ===
  Input:          Q2 = 5.00000 GeV^2,  Ebeam = 10.70000 GeV
  W               = 2.42000 GeV
  nu              = 5.31620 GeV
  E'              = 5.38380 GeV
  theta_e         = 16.94200 deg
  |q|             = 5.71148 GeV
  beta_CM         = 0.85888
  p* (CM pion)    = 1.02272 GeV/c
  theta* (CM)     = 0.000 deg  [forward production]
  -----------------
  E_pi (lab)      = 5.15934 GeV
  p_pi (lab)      = 5.10279 GeV/c
  theta_pi (lab)  = 15.78540 deg
  t (= t_min)     = -0.39685 GeV^2
  k_pi            = 0.60869 GeV  [|q| - |p_pi|]
===========================
```

Followed by a scan of `t` vs `θ*` over the full CM angular range.

## Physics and definitions

| Symbol | Definition |
|---|---|
| Q² | Virtuality: Q² = 4EE′ sin²(θ_e/2), always positive |
| ν | Energy transfer: ν = E − E′ |
| \|**q**\| | Photon 3-momentum: √(ν² + Q²) |
| W | Hadronic invariant mass: W² = M_p² + 2M_p ν − Q² |
| p\* | CM pion momentum (two-body phase space formula) |
| θ\* | CM pion angle relative to **q** direction |
| t | Mandelstam variable: t = (q − p_π)², always negative |
| t_min | Value of t at θ\* = 0 (kinematic lower bound on \|t\|) |
| k_π | Longitudinal surplus: \|**q**\| − \|**p**_π\| |

### The Mandelstam variable t

`t` measures the four-momentum transfer from the virtual photon to the pion.
It is always negative. Its physical range at fixed W and Q² is:

```
t_min(W, Q², E)  ≤  t  ≤  0
```

`t_min` is a **hard kinematic boundary**: no event at a given (W, Q², E) can
have |t| smaller than |t_min|. It is reached when θ\* = 0 (pion goes
perfectly forward along **q** in the CM). As |t| increases the pion deviates
from the photon direction; the cross section falls steeply.

The approximate relation near forward production is:

```
−t  ≈  (p* sin θ*)²  =  k_T,CM²
```

so a cut on |t| is roughly a cut on the CM transverse momentum of the pion.

### Boost recipe

The hadronic CM frame moves along **q** with β_CM = |**q**| / (M_p + ν).
The boost is a **vector boost along q̂**, not along the beam axis. The
implementation decomposes the CM pion momentum into components parallel and
perpendicular to **q̂**, applies the standard Lorentz boost to the parallel
component, leaves the perpendicular component unchanged, then projects back
onto the lab (x, z) axes. No post-hoc rotation is needed or applied.

## Known issue in reference table

The reference kinematic table (HMS-electrons / SHMS-pions settings) has an
internal inconsistency in the Q² = 8.5 GeV² row:

| Q² | W (table) | t_min at that W | t (table) |
|---|---|---|---|
| 5.0 | 2.42 | −0.397 | −0.40 ✓ |
| 6.5 | 2.72 | −0.402 | −0.40 ✓ |
| 7.5 | 2.89 | −0.410 | −0.40 ✓ |
| **8.5** | **2.95** | **−0.462** | **−0.40 ✗** |

At W = 2.95 GeV, Q² = 8.5 GeV², E = 10.7 GeV, every produced pion has
|t| ≥ 0.462 GeV². The value t = −0.40 in the table is not reachable at
those kinematics. The physically correct W for t_min = −0.40 at Q² = 8.5
GeV² is **W ≈ 3.091 GeV**. Use Mode B to get the self-consistent values.

## Files

| File | Description |
|---|---|
| `pionCT_kinematics_v2.C` | This macro |
| `pion_kinematics.pdf` | Detailed derivation of all formulas with physics context |
| `pion_kinematics.tex` | LaTeX source for the above |

## History

`v2` is a rewrite of `pionCT_kinematics.C` (original ChatGPT-generated
version) correcting the following issues:

- **Boost direction**: original boosted along +z then applied `RotateY`;
  correct approach is a direct vector boost along **q̂**.
- **Primary input**: original scanned over E′ to match t exactly at an
  arbitrary W; correct approach uses W as the primary input and reports
  t_min as output.
- **k_π definition**: original reported `Pt()` (transverse momentum w.r.t.
  beam); correct definition is |**q**| − |**p**_π| (longitudinal surplus).
- **t_min vs t_cut**: the t column in the reference table is t_min (a
  kinematic output), not a cut value imposed as an input.
