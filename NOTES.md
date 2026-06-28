# Multibody Dynamics — Study Notes

> Reconstructed from this repo, mapped one-to-one to **Ahmed A. Shabana,
> *Computational Dynamics* / *Dynamics of Multibody Systems***.
> The code implements the **augmented formulation** (Lagrange multipliers +
> constraint Jacobian) for a **2D / planar** system using an **absolute
> Cartesian coordinate** description of each body.
>
> Notation here follows Shabana: bold-ish vectors written plainly, an overbar
> `ū` means "expressed in the body's local frame", a dot means time derivative,
> and a subscript `q` means partial derivative w.r.t. the generalized
> coordinates (i.e. a Jacobian).

---

## Why this field exists (read this first)

**The one-line point of computational multibody dynamics:** it turns "deriving
the equations of motion" from a one-off feat of human cleverness into a
*systematic algorithm a computer can run on any mechanism*. That shift — from
bespoke derivation to a general, repeatable procedure — is the whole reason the
discipline exists. Everything below is that procedure, applied to a 3-link toy.

Why that matters:

- **Hand derivation doesn't scale; the algorithm does.** You can derive a double
  pendulum by hand. A car suspension, a humanoid, a spacecraft with deploying
  panels — never. The *same* recipe here assembles `M, C_q, Q_e, Q_d` whether you
  have 3 bodies or 3,000. The triple pendulum is small enough to see the recipe
  you couldn't otherwise watch run.
- **Most real dynamics have no closed form.** This pendulum is *chaotic* — there
  is no formula `θ(t)`. Numerical integration is the only way to an answer, and
  computational dynamics is the craft of making that integration correct, stable,
  and efficient.
- **It's the engine inside every physics tool.** MuJoCo, Bullet, Drake,
  Pinocchio, Adams, Simscape, game engines — all implement these formulations.
  This repo is a tiny hand-built version of what they do internally. Shabana
  teaches *the engine behind the engine*.
- **It's the foundation of robotics.** Control, state estimation, planning, and
  model-based learning all need a *computable* equation of motion. You can't do
  model-based control of a system whose dynamics you can't compute.

> Classical mechanics (Lagrange, Hamilton) gives you the **principle**;
> computational dynamics makes that principle **executable** — general, scalable,
> numerically robust enough to run on real machines for real systems.

### Where any one equation sits — the layer map

Knowing which layer you're on prevents most conceptual confusion (e.g. "is the
augmented form *the* equation of motion?" — no, it's one *formulation* of it):

```
PRINCIPLE        Newton's 2nd law / Hamilton's principle / D'Alembert
   │             ── the actual physics; unique, not negotiable
   ▼
GOVERNING EQ.    rigid multibody → constrained Newton–Euler / Lagrange
  (per domain)   fluids → Navier–Stokes   heat → diffusion eq.   (all = F=ma specialized)
   │
   ▼
FORMULATION      augmented (this repo) | embedded/minimal | recursive O(n) |
   │             Kane | Hamiltonian | penalty ...  ── a CHOICE of how to express/solve
   ▼
MODEL            THIS triple pendulum: specific M, C_q, Q_e, Q_d
   │
   ▼
SCHEME           RK4, Newton–Raphson, the 15×15 augmented solve  ── how it runs on a machine
```

The **augmented formulation** this repo implements lives at the *formulation*
layer: one valid way to express the dynamics, chosen because it's the most
transparent to build from scratch. Its main siblings — all describing the *same*
motion, differing only in cost and bookkeeping:

| Formulation | Idea | Trade-off |
|---|---|---|
| **Augmented** (this repo) | keep all coords + explicit multipliers `λ` | uniform, easy to assemble, gives reactions free; larger, indefinite system |
| **Embedded / minimal** | eliminate `λ`, reduce to independent coords (the 3 angles) | small system; messier to derive, loses reaction forces |
| **Recursive O(n)** (Featherstone) | propagate body-to-body | linear-time; the robotics workhorse (see Section 13) |
| **Kane / Gibbs–Appell** | generalized speeds, skip constraint forces | compact; common in spacecraft/robotics |
| **Hamiltonian** | first-order in position + momentum | structure/energy-preserving integration |
| **Penalty** | replace hard constraints with stiff springs | no `λ`, trivial to code; stiff, approximate |

Section 13 then covers the *extensions* of this same machinery (3D, stabilization,
flexible bodies) — the "what to learn next" axis, orthogonal to the formulation
choice above.

---

## 0. The big picture

A constrained multibody system is described by:

1. **A set of generalized coordinates** `q` (here: absolute position + angle of
   every body — *not* joint angles).
2. **A set of algebraic constraint equations** `C(q, t) = 0` that encode the
   joints (these remove the redundant DOF).
3. **The equations of motion**, which for the augmented formulation are a
   coupled system of the dynamics (Newton–Euler) and the *second derivative*
   of the constraints, solved together for accelerations and reaction forces.

The simulation loop each timestep is:

```
position analysis (Newton–Raphson)   →  enforce C(q)=0
velocity analysis (linear solve)     →  enforce Cq·q̇ = 0
acceleration / dynamics (augmented)  →  solve for q̈ and λ
numerical integration (RK4)          →  advance independent coords to t+Δt
```

This is exactly the sequence in `main_tp.py:mainProg()`.

---

## 1. Body coordinates & the rotation (transformation) matrix

### Concept
Each planar rigid body `i` is fully located by **three** generalized coordinates:
the global position of its center of mass and its orientation angle.

```math
q_i = \begin{bmatrix} R_x^i \\ R_y^i \\ \theta^i \end{bmatrix}
```

For `nb` bodies the full coordinate vector has `n = 3·nb`. Here 3 links →
**n = 9**. (See `main_tp.py:77`, `link2index` in `calcModuleTP.py:46`.)

> Indexing convention in code: for link `i`,
> `x → 3(i-1)`, `y → 3(i-1)+1`, `theta → 3(i-1)+2`.

### The planar transformation matrix `A(θ)`
Rotates a vector from the body (local) frame into the global frame
(Shabana eq. for the planar direction cosine matrix):

```math
A(\theta) = \begin{bmatrix} \cos\theta & -\sin\theta \\ \sin\theta & \cos\theta \end{bmatrix}
```

→ `calcModuleTP.py:ATransformMatrix`

### Its derivative w.r.t. the angle, `A_θ`
Used everywhere in the Jacobian and the quadratic-velocity term:

```math
A_\theta = \frac{\partial A}{\partial \theta}
= \begin{bmatrix} -\sin\theta & -\cos\theta \\ \cos\theta & -\sin\theta \end{bmatrix}
```

→ `calcModuleTP.py:ATransformMatrixTHETA`

**Key identity used by the code (worth memorizing):**

```math
\frac{\partial^2 A}{\partial \theta^2} = -A(\theta)
```

This is *why* the quadratic-velocity vector `Qd` (Section 6) is built from
`A` itself, not `A_θ`.

---

## 2. Position of an arbitrary point on a body (local → global)

The global position of a material point `P` fixed on body `i` whose local
position is `ū^P` (constant in the body frame):

```math
r_i^P = R_i + A(\theta_i)\,\bar{u}_i^P
```

→ `calcModuleTP.py:local2global`

Here the *joint points* are defined as constant local vectors from each link's
COM to its ends (`main_tp.py:56–60`), e.g. for a rod of length `L` the COM is at
mid-length so the two endpoints are `ū = [0, ±L/2]ᵀ`.

---

## 3. Constraint equations `C(q) = 0`

### Concept
Joints are expressed as algebraic constraints. Each scalar equation removes one
DOF. With `n = 9` coordinates and `nc = 6` constraints the system has
`n − nc = 3` **degrees of freedom** — the three link angles, which are the
*independent* coordinates.

### Revolute (pin) joint
The standard planar revolute joint forces two points — one on each body — to
coincide. It supplies **2 scalar constraints**:

```math
C^{rev} = r_i^P - r_j^P = R_i + A_i\bar u_i^P - R_j - A_j\bar u_j^P = \mathbf{0}
```

→ `constraintModuleTP.py:revolutJoint`

### This model's three joints (`constraintModuleTP.py:constraintEquation`)
| Joint | Type | Constraint | Eqns |
|-------|------|-----------|------|
| A | pin to **ground** | `r_1^A = 0` (a fixed point in space) | 1–2 |
| B | revolute link1–link2 | `r_1^B = r_2^B` | 3–4 |
| C | revolute link2–link3 | `r_2^C = r_3^C` | 5–6 |

Total `nc = 6`. The ground pin is just a revolute joint where one "body" is the
fixed inertial frame, so its constraint collapses to `−r_1^A = 0`.

---

## 4. The constraint Jacobian `C_q`

### Concept
`C_q = ∂C/∂q` is the `nc × n` (here **6×9**) matrix of partial derivatives of
the constraints w.r.t. the coordinates. It is the single most important object
in the formulation — it appears in the position, velocity, acceleration, and
force computations.

For a revolute constraint the block structure for the two bodies is:

```math
C_q = \big[\; \dots,\ \underbrace{I_{2\times2}}_{\partial/\partial R_i},\ \underbrace{A_{\theta,i}\,\bar u_i^P}_{\partial/\partial \theta_i},\ \dots,\ \underbrace{-I_{2\times2}}_{\partial/\partial R_j},\ \underbrace{-A_{\theta,j}\,\bar u_j^P}_{\partial/\partial \theta_j},\ \dots \;\big]
```

→ `constraintModuleTP.py:jacobianMatrix` (builds the 6×9 matrix block by block).

### Coordinate partitioning (the heart of the DOF reduction)
Shabana splits `q` into **dependent** `q_d` and **independent** `q_i`
coordinates and correspondingly partitions the Jacobian:

```math
C_q\,q = \big[\,C_{q_d}\ \ C_{q_i}\,\big]\begin{bmatrix} q_d \\ q_i \end{bmatrix}
```

- **Independent** = the three link angles `θ₁, θ₂, θ₃` (the 3 DOF, integrated in time).
- **Dependent** = the six Cartesian positions `Rx, Ry` of each link (solved from constraints).

`C_{q_d}` is the square `6×6` block that must be invertible (this requires a
valid, non-singular choice of dependent coordinates).
→ partitioning in `jacobianMatrix` returns `Cq, Cq_dep, Cq_indep`.

---

## 5. Kinematic analysis (position, velocity, acceleration)

### (a) Position analysis — Newton–Raphson
Given the independent coordinates, the dependent ones must satisfy `C(q)=0`.
Because `C` is nonlinear in `θ`, solve iteratively:

```math
C_{q_d}\,\Delta q_d = -C(q) \qquad\Longrightarrow\qquad q_d \leftarrow q_d + \Delta q_d
```

Iterate until `‖Δq_d‖ < ε`.
→ `constraintModuleTP.py:positionAnalysis`, driven by the `while` loop in
`main_tp.py:108`.

### (b) Velocity analysis — linear
Differentiate `C(q)=0` once in time. For time-independent (scleronomic)
constraints:

```math
C_q\,\dot q = 0 \;\;\Longrightarrow\;\; C_{q_d}\dot q_d + C_{q_i}\dot q_i = 0
\;\;\Longrightarrow\;\; \dot q_d = \underbrace{-C_{q_d}^{-1} C_{q_i}}_{C_{di}}\;\dot q_i
```

→ `main_tp.py:124–129` (`Cdi = inv(-Cq_dep) · Cq_indep`).

### (c) Acceleration analysis
Differentiate once more. The accelerations satisfy:

```math
C_q\,\ddot q = Q_d
```

where `Q_d` is the **quadratic velocity vector** (Section 6). This is the lower
block row of the augmented system in Section 7.

---

## 6. The quadratic velocity vector `Q_d`

### Concept
When you differentiate `C_q q̇ = 0` again you get
`C_q q̈ + (C_q q̇)_q q̇ = 0`, so the right-hand side of the acceleration
equation is:

```math
Q_d = -\left(C_q \dot q\right)_q \dot q \quad(\text{plus time terms, zero here})
```

For a revolute joint, carrying out the differentiation and using the identity
`A_θθ = −A`, the messy term collapses to a clean form:

```math
Q_d = A(\theta_i)\,\dot\theta_i^{\,2}\,\bar u_i^P \;-\; A(\theta_j)\,\dot\theta_j^{\,2}\,\bar u_j^P
```

i.e. each body contributes a **centripetal** `ω²·(rotated local vector)` term.
→ `constraintModuleTP.py:QdCalc1` (ground joint, single body) and `QdCalc2`
(two-body revolute). Assembled in `main_tp.py:250–253`.

> This is the single neatest payoff of the `A_θθ = −A` identity: the acceleration
> RHS uses the plain rotation matrix `A`, scaled by angular velocity squared.

---

## 7. Equations of motion — the **augmented formulation**

### Concept
For a constrained system, Lagrange's equations with multipliers give:

```math
M\ddot q + C_q^{\mathsf T}\lambda = Q_e
```

where `λ` are the **Lagrange multipliers** (one per constraint) and `C_qᵀ λ`
is the **generalized constraint (reaction) force**. Coupling this with the
acceleration constraint `C_q q̈ = Q_d` gives a single linear system solved each
step for `[q̈, λ]`:

```math
\begin{bmatrix} M & C_q^{\mathsf T} \\ C_q & 0 \end{bmatrix}
\begin{bmatrix} \ddot q \\ \lambda \end{bmatrix}
=
\begin{bmatrix} Q_e \\ Q_d \end{bmatrix}
```

This is **Shabana's augmented (or "embedded-multiplier") equation** — the exact
matrix assembled in `main_tp.py:systemEquation` (the `massAugmented` matrix is
`(n+nc)×(n+nc) = 15×15`, solved by a direct inverse).

### The mass matrix `M`
For a planar rigid body the (constant, diagonal) mass matrix is:

```math
M_i = \begin{bmatrix} m_i & 0 & 0 \\ 0 & m_i & 0 \\ 0 & 0 & J_i \end{bmatrix}
```

Note: because absolute coordinates are taken at the **center of mass**, the
translation and rotation **decouple** → no inertia coupling terms, `M` is
constant and diagonal. The system `M` is the block-diagonal stack of all bodies.
→ `calcModuleTP.py:massMatrix`, assembled from `massVector` in `main_tp.py:70`.

### Rod moment of inertia (about the centroid)
```math
J = \tfrac{1}{12}\,m L^2
```
→ `calcModuleTP.py:inertiaRod`

---

## 8. Generalized forces `Q_e`

The external/applied generalized force vector. Contributions in this model:

### (a) Gravity
Acts on each body's `y` coordinate (since coordinates are at the COM, weight maps
straight into the `Ry` generalized force, no torque term):

```math
Q_e^{(y_i)} = -m_i g
```
→ `main_tp.py:230–232`

### (b) Torsional (rotational) spring at a joint
A spring between bodies `i` and `j` produces equal-and-opposite generalized
torques proportional to the relative angle:

```math
Q^{spring}_{\theta_i} = -k_r(\theta_i - \theta_j - \theta_0),\qquad
Q^{spring}_{\theta_j} = +k_r(\theta_i - \theta_j - \theta_0)
```
→ `forceModule.py:torSpring`

### (c) Rotational damper at a joint
Same structure, proportional to **relative angular velocity**:

```math
Q^{damp}_{\theta_i} = -c_r(\dot\theta_i - \dot\theta_j),\qquad
Q^{damp}_{\theta_j} = +c_r(\dot\theta_i - \dot\theta_j)
```
→ `forceModule.py:torDamp`

Each joint's spring + damper torques are summed into the `θ` slots of `Q_e`.
The middle link (link 2) receives contributions from **both** joints B and C
(`main_tp.py:246–248`).

---

## 9. Joint reaction forces from the multipliers

### Concept
A major practical reason to use the augmented formulation: the Lagrange
multipliers `λ` directly give the **joint reaction (constraint) forces** for
free, with no extra free-body diagrams. The generalized reaction force is:

```math
Q_c = -C_q^{\mathsf T}\lambda
```

The components of `λ` (solved as the lower `nc` entries of the augmented
solution vector) correspond to the force pairs enforcing each joint.
→ stored as `FReact_allTime` from `qiDotDot_lamda[n:n+nc]` in `main_tp.py:138`;
plotted as the y-axis joint reactions in figure 4.

---

## 10. Numerical integration — Runge–Kutta 4

### Concept
Only the **independent** coordinates (the 3 angles) are integrated in time;
the dependent Cartesian coordinates are recovered at every step by the
position/velocity analysis above. The 3 second-order ODEs `θ̈ = f(θ, θ̇)` are
written as 6 first-order ODEs (state `y = [θ, θ̇]`) and advanced with classic
RK4:

```math
y_{n+1} = y_n + \tfrac{\Delta t}{6}\,(k_1 + 2k_2 + 2k_3 + k_4)
```

where each `k` evaluates the full constrained dynamics (re-running `config` to
rebuild `C_q` at the trial state). → `main_tp.py:rungeKutta4_AtTimeNow`.

> Subtlety worth noting: at every RK stage the code rebuilds the Jacobian and
> re-solves the augmented system to extract `θ̈`. This is correct but expensive —
> it's the "honest" textbook way, with no shortcuts, which matches the repo's
> stated goal of building from scratch.

---

## 11. The complete per-timestep algorithm (as coded)

```
for each time step:
    # known: independent θ, θ̇ at t
    1. Position analysis  (Newton–Raphson):  C_qd Δq_d = −C  → dependent positions
    2. Velocity analysis  (linear):          q̇_d = −C_qd⁻¹ C_qi q̇_i
    3. Acceleration/dynamics (augmented):    [M Cqᵀ; Cq 0][q̈; λ] = [Qe; Qd]
    4. Store q, q̇, q̈, and reaction forces (from λ)
    5. Integrate independent θ, θ̇ to t+Δt with RK4
```

→ `main_tp.py:mainProg` lines 103–144.

---

## 12. Concept ↔ file cross-reference

| Concept (Shabana) | Equation | Code |
|---|---|---|
| Body coordinates `[Rx,Ry,θ]` | Section 1 | `link2index` |
| Rotation matrix `A(θ)` | Section 1 | `ATransformMatrix` |
| `A_θ`, identity `A_θθ=−A` | Section 1 | `ATransformMatrixTHETA` |
| Local→global point `r=R+Aū` | Section 2 | `local2global` |
| Revolute constraint `r_i=r_j` | Section 3 | `revolutJoint`, `constraintEquation` |
| Constraint Jacobian `C_q` | Section 4 | `jacobianMatrix` |
| Dependent/independent split | Section 4 | partition in `jacobianMatrix` / `mainProg` |
| Newton–Raphson position | Section 5a | `positionAnalysis` |
| Velocity `q̇_d = −C_qd⁻¹C_qi q̇_i` | Section 5b | `mainProg` (`Cdi`) |
| Quadratic velocity `Q_d` | Section 6 | `QdCalc1`, `QdCalc2` |
| Augmented EOM `[M Cqᵀ;Cq 0]` | Section 7 | `systemEquation` |
| Mass matrix `diag(m,m,J)` | Section 7 | `massMatrix` |
| Rod inertia `mL²/12` | Section 7 | `inertiaRod` |
| Gravity / spring / damper `Q_e` | Section 8 | `systemEquation`, `forceModule` |
| Reaction forces from `λ` | Section 9 | `mainProg` (`FReact_allTime`) |
| RK4 integration | Section 10 | `rungeKutta4_AtTimeNow` |

---

## 13. What to revisit next (toward "mastering robotics")

The repo nails the **planar absolute-coordinate augmented formulation**. Natural
extensions from Shabana / standard MBD, in roughly increasing difficulty:

- **3D rigid bodies**: orientation needs Euler parameters (quaternions) →
  the rotation matrix and the `A_θ`/`Q_d` machinery generalize but get heavier.
- **Baumgarte / coordinate-projection stabilization**: the pure augmented
  solve lets constraint error `C(q)` drift numerically; stabilization or
  index reduction keeps `‖C‖` bounded.
- **Recursive / minimal-coordinate methods** (Featherstone's Articulated Body
  Algorithm): O(n) joint-space dynamics — the formulation robotics control
  actually uses, vs. the O(n³) dense augmented solve here.
- **Flexible bodies** (the floating-frame-of-reference formulation): Shabana's
  later chapters, where `ū` is no longer constant.
- **Better integrators**: implicit / stiff solvers for stiff springs, and
  proper handling of the DAE (index-3) nature of the constrained system.

The conceptual through-line: this repo is the *dense, full-coordinate, DAE*
view of dynamics; robotics control favors the *minimal-coordinate, recursive*
view. Knowing both — and why they agree — is the actual mastery.

---

## 14. How much of Shabana does this repo actually cover?

Short version: **by page count, ~15–25% of the books. By the dependency tree of
concepts, this is the foundational ~20% that everything else is built on.**
That's the "vital few" — the part you cannot skip, fake, or shortcut.

### Why the 20% here *is* 80% of the concepts
The augmented-formulation pipeline below is the conceptual spine of all
computational multibody dynamics:

```
constraints C(q)=0  →  Jacobian C_q  →  augmented EOM [M Cqᵀ; Cq 0]
                    →  solve for q̈ and λ  →  integrate
```

Once this is built from scratch, the rest of Shabana is **extensions of these
same objects**, not new paradigms:

- **3D / spatial dynamics** → same formulation; `A(θ)` becomes a
  quaternion-based rotation matrix, inertia becomes a 3×3 tensor.
- **Other joint types** (prismatic, cam, gear, spatial) → same `C(q)=0`
  machinery, just different constraint rows and Jacobian blocks.
- **Flexible / deformable bodies** (floating frame of reference, much of
  book 2) → same EOM; the local vector `ū` stops being constant.

So this repo doesn't sit at "20% done" — it sits at roughly the **first 50–60%
of the conceptual dependency tree**, which happens to be a smaller fraction of
the page count.

### Coverage map

| Shabana topic | Status |
|---|---|
| Planar kinematics, `A(θ)`, local→global | ✅ Fully |
| Absolute Cartesian / generalized coordinates | ✅ Fully |
| Constraints, Jacobian `C_q`, dep/indep split | ✅ Fully |
| Position / velocity / acceleration analysis | ✅ Fully |
| Augmented formulation, Lagrange multipliers | ✅ Fully — the centerpiece |
| Mass matrix, generalized forces, reactions from λ | ✅ Fully |
| Analytical mechanics (virtual work, d'Alembert, Lagrange derivation) | ⚠️ Used, not derived |
| Spatial / 3D dynamics (quaternions, 3D inertia) | ❌ |
| Constraint stabilization (Baumgarte, DAE index, projection) | ❌ |
| Recursive / minimal-coordinate methods (Featherstone O(n)) | ❌ |
| Flexible bodies (floating frame of reference) — much of book 2 | ❌ |
| Joint variety (prismatic, cam, gears, spatial) | ❌ (revolute only) |

### The honest boundary
The ❌ rows are not "advanced footnotes" — they are the other 80% of the books,
and they are where the robotics-specific depth lives (O(n) dynamics for control,
3D for real manipulators, flexibility for real structures). But **none of them
make sense without Sections 1–11 above**, which is exactly why this repo is good proof
that the foundation is real. It's the 20% that unlocks reading the rest.
