# Model Derivation — Triple Pendulum (paper work *before* the code)

> This is the "do the math on paper first" worksheet for **this exact repo's
> model**. It is the companion to `NOTES.md`: where `NOTES.md` explains the
> general theory, this file applies it step-by-step to the specific 3-link
> system until every symbol has a concrete expression you can transcribe into
> Python. Section numbers in brackets like [Section 4] point back to `NOTES.md`.
>
> Workflow captured here:
> **sketch → coordinates → body points → constraints → Jacobian (by hand)
> → quadratic-velocity vector → mass matrix → forces → augmented EOM
> → solve plan → map symbols to code.**

---

## The destination: the augmented equation of motion

Everything in this worksheet exists to fill in **one matrix equation** — the
cardinal equation of constrained multibody dynamics (the *augmented* form).
Solve it at each instant for the accelerations `q̈` and the Lagrange
multipliers `λ`:

```math
\begin{bmatrix} M & C_q^{\mathsf T} \\ C_q & 0 \end{bmatrix}
\begin{bmatrix} \ddot{q} \\ \lambda \end{bmatrix}
=
\begin{bmatrix} Q_e \\ Q_d \end{bmatrix}
```

Read it as: *(top row)* Newton's law with the joint reaction forces `C_qᵀλ`
added; *(bottom row)* the requirement that accelerations keep the joints
satisfied (`C_q q̈ = Q_d`). Two unknown blocks, two equation blocks.

### Intuition — what this equation is really saying

Forget the matrix for a second. The augmented form is the answer to a very human
question: **"the links want to fall under gravity, but the joints won't let them
fly apart — so how does everything actually accelerate this instant?"** Two
things must be true at once, and the two rows are exactly those two truths.

**Top row — "obey Newton, *plus* whatever force the joints have to apply."**
A free body would just obey `M q̈ = Q_e` (mass × acceleration = applied force).
But our links aren't free — they're pinned. The pins must push/pull to hold the
joints together, and that extra push is `C_qᵀλ`. So the top row is literally
*Newton's law with an extra "joint force" term whose size we don't yet know*:

```math
\underbrace{M\ddot q}_{\text{mass}\times\text{accel}} = \underbrace{Q_e}_{\text{gravity, springs}} \;-\; \underbrace{C_q^{\mathsf T}\lambda}_{\text{force the joints exert}}
```

Think of `λ` as **"how hard each joint is pulling right now"** and `C_qᵀ` as
**"which direction that pull points."** `C_q` (the Jacobian) already encodes the
geometry of each joint — its transpose turns a pull-strength `λ` into an actual
force vector on the bodies. (Why must the joint force have the form `C_qᵀλ`?
Because a joint can only push *along the directions it constrains* — it does no
work sliding *along* the allowed motion. That "no free work" idea is D'Alembert's
principle; see `NOTES.md` Section 7.)

**Bottom row — "and the accelerations must keep the joints closed."**
The pins are welded shut: point B on link 1 *is* point B on link 2, forever. If
that's true for all time, it's true for position, velocity, *and* acceleration.
The position version (`C = 0`) and velocity version (`C_q q̇ = 0`) are already
guaranteed by how we start and step; the bottom row enforces the **acceleration**
version, `C_q q̈ = Q_d`. In words: *the accelerations aren't free to be anything —
they must be exactly the ones that keep the joints from pulling apart or
crushing together.* `Q_d` is just the leftover centripetal bookkeeping (the
`ω²r` terms) that comes from the joints rotating (Step 6).

**Why solve them together (the punchline).** Here's the chicken-and-egg that the
matrix resolves in one shot:

- You can't find the accelerations `q̈` until you know the joint forces `λ`
  (they appear in Newton's law).
- You can't find the joint forces `λ` until you know the accelerations `q̈`
  (the forces are *whatever it takes* to keep `C_q q̈ = Q_d`).

Neither can be computed first. The augmented matrix stacks both requirements and
solves for `q̈` **and** `λ` **simultaneously** — finding the one acceleration *and*
the one set of joint forces that are mutually consistent. That's the whole reason
it's a single bordered system instead of two separate steps.

**A one-line mental model:** *gravity proposes an acceleration; the joints
veto the parts that would tear them apart; `λ` is the exact strength of that
veto, and `q̈` is what's left after the veto is applied.* The top row is the
proposal, the bottom row is the veto, and solving them together is the
negotiation that happens every timestep.

Every step below is just **producing one of these ingredients** for our specific
3-link model:

| Block | What it is | Comes from |
|---|---|---|
| `M` | mass / inertia matrix | Step 7 |
| `C_q` | constraint Jacobian | Step 5 |
| `Q_e` | applied generalized forces (gravity, springs, dampers) | Step 8 |
| `Q_d` | quadratic-velocity vector (the `C_q q̈ = Q_d` right side) | Step 6 |
| `q̈`, `λ` | **the unknowns** — accelerations & reaction multipliers | Step 9 solve |

For our model the assembled system is **15×15** (`n + nc = 9 + 6`). Keep this
target in view: each section is handing you one piece of it. Steps 1–4 build the
scaffolding (`q`, points, constraints) the ingredients depend on; Steps 5–8 are
the ingredients themselves; Step 9 assembles and solves; Steps 10–11 wrap it in
time integration and code.

---

## Step 0 — Sketch and inventory the system

The model that the rest of this derivation is built on:

![Triple pendulum model](img/triplePendulum1.png)

*Joints A (pin to ground), B and C (revolute, each with a torsional spring `Kr`
and rotational damper `Cr`); body-fixed frames `XⁱYⁱ`; angles `θᵢ` from the
vertical; lengths `Lᵢ`; masses `mᵢ`. Every symbol below maps to this figure.*

Text schematic of the same thing:

```
            ┌── ground (fixed point, world origin O)
            ●  A          ← pin joint A (link1 ↔ ground)
            │
            │  link 1     (mass m1, length L1)
            │
            ●  B          ← revolute joint B (link1 ↔ link2)  + torsional spring krB, damper crB
            │
            │  link 2     (mass m2, length L2)
            │
            ●  C          ← revolute joint C (link2 ↔ link3)  + torsional spring krC, damper crC
            │
            │  link 3     (mass m3, length L3)
            │
            ○             ← free end
```

**Inventory:**
- **Bodies:** 3 rigid uniform rods.
- **Joints:** A (pin to ground), B (revolute), C (revolute) → 3 joints.
- **Force elements:** gravity (all bodies); torsional spring + rotational damper at B and at C.
- **Body reference point:** each rod's **center of mass** (mid-length) — chosen so the mass matrix is constant & diagonal (see `NOTES.md` Section 7 and the COM discussion).

**Expected DOF (sanity check before any math):**
`DOF = 3·(bodies) − (constraint equations)`.
A pin and each revolute supply 2 scalar constraints → `3·3 − 3·2 = 9 − 6 = 3`.
The 3 DOF will be the three link angles `θ₁, θ₂, θ₃`.

---

## Step 1 — Choose generalized coordinates [Section 1]

Absolute (Cartesian) coordinates, **3 per body** = COM position + orientation:

```math
q = \big[\, \underbrace{R_{x1},R_{y1},\theta_1}_{\text{body 1}},\ \underbrace{R_{x2},R_{y2},\theta_2}_{\text{body 2}},\ \underbrace{R_{x3},R_{y3},\theta_3}_{\text{body 3}} \,\big]^{\mathsf T}, \qquad n = 9
```

Angles `θ_i` measured from the global vertical (world +y), positive CCW.

---

## Step 2 — Locate the body-fixed joint points `ū` [Section 2]

Each rod's COM is at mid-length, so its two ends are `±L/2` along the body's
local y-axis. List every point that participates in a joint:

| Point | On body | Local vector `ū` | Meaning |
|---|---|---|---|
| A | 1 | `[0,  +L₁/2]` | top end of link 1 (to ground) |
| B | 1 | `[0,  −L₁/2]` | bottom end of link 1 |
| B | 2 | `[0,  +L₂/2]` | top end of link 2 |
| C | 2 | `[0,  −L₂/2]` | bottom end of link 2 |
| C | 3 | `[0,  +L₃/2]` | top end of link 3 |

> ⚠️ **Repo discrepancy to fix on paper vs. code:** in `main_tp.py:59`,
> `u_bar_2C` is coded as `[0, −L₁/2]` (uses L₁), but the correct model value is
> `[0, −L₂/2]`. Harmless while all `L=1`; wrong if `L₂ ≠ L₁`. The math below
> uses the correct `L₂/2`.

Shorthand for the rest of this document: `a = L₁/2`, `b = L₂/2`, `d = L₃/2`,
and `sᵢ = sin θᵢ`, `cᵢ = cos θᵢ`.

---

## Step 3 — Global position of each joint point [Section 2]

Using `rᵢᴾ = Rᵢ + A(θᵢ) ūᵢᴾ` with
`A(θ) = [[c, −s], [s, c]]`, work out each one (this is the algebra you'd grind on paper):

```math
A(\theta)\begin{bmatrix}0\\ \pm h\end{bmatrix} = \begin{bmatrix} \mp h\, s \\ \pm h\, c \end{bmatrix}
```

| Point | Global position `r` |
|---|---|
| `r₁ᴬ` | `[ R_{x1} − a s₁ ,  R_{y1} + a c₁ ]` |
| `r₁ᴮ` | `[ R_{x1} + a s₁ ,  R_{y1} − a c₁ ]` |
| `r₂ᴮ` | `[ R_{x2} − b s₂ ,  R_{y2} + b c₂ ]` |
| `r₂ᶜ` | `[ R_{x2} + b s₂ ,  R_{y2} − b c₂ ]` |
| `r₃ᶜ` | `[ R_{x3} − d s₃ ,  R_{y3} + d c₃ ]` |

### Why the `ū` vectors are `[0, ±h]` (and other "free" simplifications)

None of these simplifications change the physics — they're **choices of where to
put the body frame**, which is yours to make. Picking them well turns ugly
algebra into the clean form above. The rule of thumb: *the physics is fixed, but
the bookkeeping is yours — spend that freedom before you start differentiating.*

**1. Align the local `Yⁱ` axis along the rod ⇒ `ū = [0, ±h]`.**
We deliberately orient each body-fixed frame so its `Yⁱ` axis runs **down the
length of the linkage**. Then every point of interest (both ends, the COM) lies
on that axis, so its `x_local = 0` and the local vector collapses to just
`[0, ±h]` — a single nonzero number instead of two. That single zero is what
makes `A(θ)ū` reduce to the tidy `[∓h s, ±h c]` template, which in turn makes the
whole Jacobian (Step 5) sparse and hand-derivable. Had we tilted the local frame
arbitrarily, every `ū` would have two nonzero entries and twice the trig.

**2. Put the reference point `R` at the center of mass.**
Free to place anywhere on the body (see the anchor-choice discussion in `NOTES.md`), but the COM is
special: it **decouples translation from rotation**, making the mass matrix
constant and diagonal `diag(m, m, J)` (Step 7) with no velocity-dependent inertia
terms. Any other anchor is still valid physics but drags coupling terms into `M`.
For a uniform rod the COM is the geometric mid-point, which is *also* why the two
ends come out as the symmetric `±L/2` — choices (1) and (2) reinforce each other.

**3. Anchor the ground pin at the world origin ⇒ joint A is `r₁ᴬ = 0`.**
We're free to place the inertial frame anywhere; putting `O` exactly at the pin
makes the ground constraint `r₁ᴬ − 0 = 0` instead of `r₁ᴬ − (some constant) = 0`.
The constant vanishes, and so does its contribution to every derivative.

**4. Measure all angles `θᵢ` from the same global vertical.**
Using one common reference axis (not joint-relative angles) keeps each `A(θᵢ)`
independent of the others, so the Jacobian blocks don't chain-multiply. (This is
the absolute-coordinate choice; it trades more coordinates for simpler, decoupled
derivatives — see `NOTES.md` Section 13.)

**5. Uniform rod ⇒ closed-form inertia `J = mL²/12`.**
Assuming a slender uniform rod lets the moment of inertia be a one-line formula
instead of an integral, and centers it at the mid-point (reinforcing #2).

> The pattern across all five: *exploit the freedoms the formulation gives you
> (frame origin, frame orientation, inertial-frame placement, angle reference)
> to bury as many zeros and constants in the equations as possible — **before**
> you differentiate.* Every zero you place by choice is a term you never have to
> derive, code, or debug.

---

## Step 4 — Write the constraint equations `C(q) = 0` [Section 3]

### The canonical form first (where every equation below comes from)

Each joint is an instance of one **canonical constraint**. Derive the general
form once, then just substitute the specific bodies/points.

**Revolute (pin) joint — the master equation.** A revolute joint between bodies
`i` and `j` says: *the point P on body i and the point P on body j are the same
point in space.* That's a 2-D vector equation (⇒ 2 scalar constraints):

```math
\boxed{\;C^{rev}_{ij} = r_i^{P} - r_j^{P} = \big(R_i + A(\theta_i)\,\bar u_i^{P}\big) - \big(R_j + A(\theta_j)\,\bar u_j^{P}\big) = \mathbf{0}\;}
```

**Ground (fixed-pin) joint — the special case.** When one of the two bodies is
the inertial frame, its point is just a constant location `c` (no `R`, no `A`),
so the master equation degenerates to:

```math
C^{grd}_{i} = r_i^{P} - c = R_i + A(\theta_i)\,\bar u_i^{P} - c = \mathbf{0}
```

Everything in Step 4 is one of these two. To instantiate, you only plug in:
**(1)** which bodies, **(2)** their local points `ū` (from Step 2), and **(3)**
the `r = R + Aū` expansions (already done in Step 3). The substitution *is* the
derivation — there's no new calculus, just bookkeeping.

> General → specific recipe:
> 1. Write the canonical form for the joint type.
> 2. Replace `rᵢᴾ`, `rⱼᴾ` with their Step-3 expansions.
> 3. Split the 2-D vector equation into its x- and y-rows.

---

The ground pin fixes A to the world origin (`c = O = [0,0]`); each revolute makes
the two coincident points equal. Applying the recipe gives six scalar equations:

**Joint A** = `C^grd₁` with `P = A`, `c = 0` (`r₁ᴬ = 0`):
```math
\begin{aligned}
C_1 &: R_{x1} - a\,s_1 = 0\\
C_2 &: R_{y1} + a\,c_1 = 0
\end{aligned}
```

**Joint B** = `C^rev₁₂` with `P = B` (`r₁ᴮ − r₂ᴮ = 0`):
```math
\begin{aligned}
C_3 &: R_{x1} - R_{x2} + a\,s_1 + b\,s_2 = 0\\
C_4 &: R_{y1} - R_{y2} - a\,c_1 - b\,c_2 = 0
\end{aligned}
```

**Joint C** = `C^rev₂₃` with `P = C` (`r₂ᶜ − r₃ᶜ = 0`):
```math
\begin{aligned}
C_5 &: R_{x2} - R_{x3} + b\,s_2 + d\,s_3 = 0\\
C_6 &: R_{y2} - R_{y3} - b\,c_2 - d\,c_3 = 0
\end{aligned}
```

So `C = [C₁ … C₆]ᵀ`, `nc = 6`. DOF `= 9 − 6 = 3` ✓ (matches Step 0).

> **Sign note for code-matching:** the repo stores joint A as `−r₁ᴬ` rather than
> `+r₁ᴬ` (`constraintModuleTP.py:revolutJoint` / `constraintEquation`). That just
> negates rows 1–2 (and the matching `Q_d`/λ), which is physically identical.
> This worksheet keeps the natural `+` sign; flip rows 1–2 if you want a
> byte-for-byte match.

---

## Step 5 — Derive the Jacobian `C_q` by hand [Section 4]

`C_q = ∂C/∂q` is **6×9**. Differentiate each `Cᵢ` w.r.t. every coordinate.
Columns ordered `[R_{x1}, R_{y1}, θ₁ | R_{x2}, R_{y2}, θ₂ | R_{x3}, R_{y3}, θ₃]`.

Useful sub-result (this is the `A_θ ū` block from Section 4):
`∂(A(θ)ū)/∂θ = A_θ(θ) ū`, e.g. `∂r₁ᴬ/∂θ₁ = [−a c₁, −a s₁]`.

```math
C_q=\begin{bmatrix}
1 & 0 & -a c_1 & 0 & 0 & 0 & 0 & 0 & 0\\
0 & 1 & -a s_1 & 0 & 0 & 0 & 0 & 0 & 0\\
1 & 0 &  a c_1 & -1 & 0 & b c_2 & 0 & 0 & 0\\
0 & 1 &  a s_1 & 0 & -1 & b s_2 & 0 & 0 & 0\\
0 & 0 & 0 & 1 & 0 & b c_2 & -1 & 0 & d c_3\\
0 & 0 & 0 & 0 & 1 & b s_2 & 0 & -1 & d s_3
\end{bmatrix}
```

This is the matrix `constraintModuleTP.py:jacobianMatrix` builds block-by-block.
Each `±I₂` is a translation block; each trig column is a `±A_θ(θ) ū` block.

**Partition into dependent / independent [Section 4]:** pick the 3 angles as
independent. Then:
- `C_{q_i}` = columns `{θ₁, θ₂, θ₃}` (cols 3, 6, 9) → **6×3**
- `C_{q_d}` = the other 6 columns (the Cartesian `R`'s) → **6×6** (must be invertible)

```math
C_{q_d}=\begin{bmatrix}
1&0&0&0&0&0\\0&1&0&0&0&0\\
1&0&-1&0&0&0\\0&1&0&-1&0&0\\
0&0&1&0&-1&0\\0&0&0&1&0&-1
\end{bmatrix},\qquad
C_{q_i}=\begin{bmatrix}
-a c_1 & 0 & 0\\
-a s_1 & 0 & 0\\
 a c_1 & b c_2 & 0\\
 a s_1 & b s_2 & 0\\
 0 & b c_2 & d c_3\\
 0 & b s_2 & d s_3
\end{bmatrix}
```

(Note `C_{q_d}` is constant — a nice property of choosing the Cartesian coords as
dependent.)

---

## Step 6 — Velocity & acceleration kinematics; derive `Q_d` [Sections 5–6]

**Velocity** (differentiate `C=0` once): `C_q q̇ = 0`, so
`q̇_d = −C_{q_d}⁻¹ C_{q_i} q̇_i`.

**Acceleration** (differentiate again): `C_q q̈ = Q_d`, with the
quadratic-velocity vector built from the `A_θθ = −A` identity [Section 6]. For each joint
`Q_d = θ̇ᵢ² A(θᵢ) ūᵢ − θ̇ⱼ² A(θⱼ) ūⱼ` (single body for A). Plugging the points:

```math
\begin{aligned}
Q_d^{A} &= \dot\theta_1^{\,2}\,[\,-a s_1,\ a c_1\,]\\
Q_d^{B} &= \dot\theta_1^{\,2}\,[\,a s_1,\ -a c_1\,] - \dot\theta_2^{\,2}\,[\,-b s_2,\ b c_2\,]\\
Q_d^{C} &= \dot\theta_2^{\,2}\,[\,b s_2,\ -b c_2\,] - \dot\theta_3^{\,2}\,[\,-d s_3,\ d c_3\,]
\end{aligned}
```

Stack: `Q_d = [Q_d^A ; Q_d^B ; Q_d^C]` (6×1). Each entry is a centripetal
`ω²·(rotated arm)` term — no `θ̈` appears. → `QdCalc1`, `QdCalc2`.

---

## Step 7 — Mass matrix `M` [Section 7]

COM coordinates ⇒ translation/rotation decouple ⇒ constant diagonal **9×9**:

```math
M = \mathrm{diag}(\,m_1, m_1, J_1,\ m_2, m_2, J_2,\ m_3, m_3, J_3\,),
\qquad J_i = \tfrac{1}{12} m_i L_i^2
```

→ `massMatrix`, `inertiaRod`.

---

## Step 8 — Generalized forces `Q_e` (9×1) [Section 8]

Start at zero, then add each effect into the right slot.

**Gravity** → into the `R_y` slots only (force at COM, no torque):
```math
Q_e^{(R_{y1})} = -m_1 g,\quad Q_e^{(R_{y2})} = -m_2 g,\quad Q_e^{(R_{y3})} = -m_3 g
```

**Torsional springs + rotational dampers** → into the `θ` slots. With rest angle
`θ₀ = 0`, joint B couples (1,2) and joint C couples (2,3):
```math
\begin{aligned}
Q_e^{(\theta_1)} &= -k_{rB}(\theta_1-\theta_2) - c_{rB}(\dot\theta_1-\dot\theta_2)\\
Q_e^{(\theta_2)} &= +k_{rB}(\theta_1-\theta_2) - k_{rC}(\theta_2-\theta_3) + c_{rB}(\dot\theta_1-\dot\theta_2) - c_{rC}(\dot\theta_2-\dot\theta_3)\\
Q_e^{(\theta_3)} &= +k_{rC}(\theta_2-\theta_3) + c_{rC}(\dot\theta_2-\dot\theta_3)
\end{aligned}
```

The middle link gets contributions from **both** joints (the `±` reaction pairs).
All other entries of `Q_e` are 0. → `forceModule.torSpring`, `forceModule.torDamp`,
assembled in `systemEquation`.

---

## Step 9 — Assemble the augmented equation of motion [Section 7]

Combine dynamics (`M q̈ + C_qᵀ λ = Q_e`) with acceleration constraints
(`C_q q̈ = Q_d`) into one **(n+nc) × (n+nc) = 15×15** linear system:

```math
\begin{bmatrix} M & C_q^{\mathsf T} \\ C_q & 0 \end{bmatrix}
\begin{bmatrix} \ddot q \\ \lambda \end{bmatrix}
=
\begin{bmatrix} Q_e \\ Q_d \end{bmatrix}
```

Solve at each instant for the 9 accelerations `q̈` and 6 multipliers `λ`.
Joint reaction forces follow from `Q_c = −C_qᵀ λ` [Section 9]. → `systemEquation`.

---

## Step 10 — The time loop, derived step by step [Sections 5, 10–11]

Only `θ₁, θ₂, θ₃` (and their rates) are integrated; everything else is recovered
from the constraints each step. This section is the full anatomy of the `for`
loop in `mainProg()` — **for each step: the goal, the derivation of its formula
from scratch, and the code**. The four derivations are really just `C(q)=0`
seen at four levels: *position* (`C=0`), *velocity* (`Ċ=0`), *acceleration*
(`C̈=0`), and *time integration*.

**Carried state and outputs.** Three vectors persist across steps: `qi`
(9 positions), `qiDot` (9 velocities), `qiDotDot_lamda` (15 = 9 accelerations +
6 multipliers). Each step fills four output arrays — these *are* the solution:

| Output array | Holds | Source |
|---|---|---|
| `q_allTime` | positions (9) | recovered `qi` |
| `v_allTime` | velocities (9) | recovered `qiDot` |
| `a_allTime` | accelerations (9) | first 9 of the augmented solve |
| `FReact_allTime` | reaction multipliers `λ` (6) | last 6 of the augmented solve |

**What is known at the top of a step.** Only the 3 angles and their rates
(`θ`, `θ̇`) are authoritative — RK4 advanced them last step. The 6 dependent
Cartesian coordinates are stale and must be rebuilt. So every step is: *recover
everything from the angles → solve → push the angles forward.*

**The loop skeleton** (labels match the code comments in `mainProg`):

```
for each timestep t:
    (a) position analysis   → dependent POSITIONS  @ t    (Newton–Raphson loop)
    (b) velocity analysis   → dependent VELOCITIES @ t    (one linear solve)
    (c) augmented solve     → ACCELERATIONS q̈ + reactions λ @ t
    (d) store q, q̇, q̈, λ
    (e) RK4 integration     → independent θ, θ̇  @ t+Δt
```

Initial conditions: set `θ₁,θ₂,θ₃` (here `π/4, π/3, π/2`), velocities 0; run one
position analysis (a) once to get a consistent starting configuration.

---

### (a) Position analysis — recover dependent positions

**Goal.** Find the 6 dependent Cartesian coords that make the joints intact
(`C(q)=0`) for the current (fixed) angles.

**Derivation.** This is **Newton–Raphson**, the first-order-approximation idea.
We have a guess where the joints are slightly *open* by the violation vector
`C ≠ 0`. Take the first-order Taylor expansion of `C` and demand the corrected
point close the joints:

```math
C(q + \Delta q) \approx C(q) + C_q\,\Delta q = 0
\;\;\Longrightarrow\;\;
\Delta q = -\,C_q^{-1}\,C(q)
```

`C` is "how far off the joints are"; `C_qd⁻¹` converts that violation into "how
far to move the coordinates to close it"; the minus sign moves *opposite* to the
violation. It's the vector form of Newton's `x ← x − f/f'`.

**Why iterate until `‖Δq‖ < ε`.** The Taylor expansion is only first-order — a
straight-line estimate of a curved function — so for nonlinear `C` one step gets
closer but not exact. Recompute `C`, step again, until the correction drops below
tolerance `ε` (= "stopped moving → `C ≈ 0`"). It's `ε` not `0` because floating
point can't hit zero. (`max_iteration` caps the loop.) **Here it converges in one
step**, because with angles fixed the constraints are *linear* in the Cartesian
unknowns — `C_qd` from Step 5 is the constant `±1` matrix, no trig — so the
first-order step is exact.

**Code.** The single Newton step (`constraintModuleTP.py`):

```python
def positionAnalysis(constraintVector, jacobianMatrix, qi):
    inverse_jacobian = np.linalg.inv(jacobianMatrix)            # C_qd⁻¹
    delta_qi = - np.matmul(inverse_jacobian, constraintVector)  # Δq = −C_qd⁻¹ · C
    delta_qi_norm = np.linalg.norm(delta_qi)                    # ‖Δq‖  (for the ε test)
    qi = qi + delta_qi                                          # q ← q + Δq
    return qi, delta_qi_norm
```

The `until ‖Δq‖ < ε` loop wrapping it (`main_tp.py`):

```python
delta_qDep_norm = 1
while delta_qDep_norm > epsilon:                               # repeat until converged
    Cq, Cq_dep, Cq_indep, constraintVect = config(qi)         # rebuild C and C_q at current guess
    q_dep = np.concatenate((qi[0:2], qi[3:5], qi[6:8]))       # the 6 dependent positions (Rx,Ry)
    q_depNew, delta_qDep_norm = conMod.positionAnalysis(      # one Newton step → new q_dep, ‖Δq‖
        constraintVect, Cq_dep, q_dep)
    if (delta_qDep_norm < epsilon) or (count > max_iteration): # converged, or safety cap
        break
qi[0:2], qi[3:5], qi[6:8] = q_depNew[0:2], q_depNew[2:4], q_depNew[4:6]  # store back
```

`constraintVect`=`C`, `Cq_dep`=`C_qd`, `epsilon`=`ε`. The slices `[0:2],[3:5],
[6:8]` are each link's `(Rx,Ry)`; the angle indices `2,5,8` are held fixed.

---

### (b) Velocity analysis — recover dependent velocities

**Goal.** Given the 3 angular velocities, find the 6 dependent Cartesian
velocities that keep the joints intact *as they move*.

**Derivation.** "Joints stay closed for all time" means `C(q(t)) = 0` for every
`t`. Differentiate once in time (chain rule). For our time-independent
constraints there is no explicit `∂C/∂t` term, so:

```math
\frac{d}{dt}\,C(q(t)) = C_q\,\dot q = 0
```

Now split the coordinates into dependent and independent columns
(`C_q q̇ = C_qd q̇_dep + C_qi q̇_indep`) and solve for the unknown dependent part:

```math
C_{qd}\,\dot q_{dep} + C_{qi}\,\dot q_{indep} = 0
\;\;\Longrightarrow\;\;
\dot q_{dep} = -\,C_{qd}^{-1}\,C_{qi}\,\dot q_{indep}
```

This is **linear and exact** — one solve, no iteration — because `C_qd`, `C_qi`
are already known from the positions found in (a), and `q̇_indep` is the known
angular velocity. (Same `C_q q̇ = 0` appears in Section 5; it is also the velocity
form whose *derivative* gives the acceleration constraint in (c).)

**Code** (`main_tp.py`):

```python
qDot_indep = np.concatenate((qiDot[2:3], qiDot[5:6], qiDot[8:9]))  # the 3 angular rates θ̇
Cdi      = np.dot(np.linalg.inv(-Cq_dep), Cq_indep)               # Cdi = −C_qd⁻¹ · C_qi
qDot_dep = np.dot(Cdi, qDot_indep)                               # q̇_dep = Cdi · θ̇
qiDot[0:2], qiDot[3:5], qiDot[6:8] = \
    qDot_dep[0:2], qDot_dep[2:4], qDot_dep[4:6]                  # store back into qiDot
```

`Cdi` is literally the matrix `−C_qd⁻¹ C_qi`.

---

### (c) Acceleration & reactions — the augmented solve

**Goal.** Find the 9 accelerations `q̈` and the 6 reaction multipliers `λ` at
this instant.

**Derivation — part 1, the acceleration constraint.** Differentiate the velocity
constraint `C_q q̇ = 0` once more in time (product rule):

```math
\frac{d}{dt}\big(C_q\,\dot q\big) = C_q\,\ddot q + \dot C_q\,\dot q = 0
\;\;\Longrightarrow\;\;
C_q\,\ddot q = -\,\dot C_q\,\dot q \;\equiv\; Q_d
```

So the accelerations are *not free* — they must satisfy `C_q q̈ = Q_d`, where
`Q_d` (the **quadratic-velocity vector**) collects the leftover `q̇`-terms. For
this model the only time-varying part of `C_q` is each `A_θ(θ)·ū` block. Its time
derivative is:

```math
\frac{d}{dt}\big[A_\theta(\theta)\,\bar u\big]
= A_{\theta\theta}\,\dot\theta\,\bar u
= -A(\theta)\,\dot\theta\,\bar u
```

using the `A_θθ = −A` identity (Step 6). Working `−Ċ_q q̇` through for a joint
point gives the centripetal term `θ̇²·A(θ)·ū` — exactly `QdCalc1` (single body)
and `QdCalc2` = `θ̇ᵢ²·Aᵢ·ūᵢ − θ̇ⱼ²·Aⱼ·ūⱼ` (two bodies).

**Derivation — part 2, stack with the dynamics.** The dynamics (Step 9 /
"destination") give `M q̈ + C_qᵀ λ = Q_e`. Stacking it with the acceleration
constraint above yields one linear system in the unknowns `q̈` and `λ`:

```math
\begin{bmatrix} M & C_q^{\mathsf T} \\ C_q & 0 \end{bmatrix}
\begin{bmatrix} \ddot q \\ \lambda \end{bmatrix}
=
\begin{bmatrix} Q_e \\ Q_d \end{bmatrix}
```

15×15 here (`9 + 6`). Invert → `q̈` (first 9) and `λ` (last 6).

**Code** (`main_tp.py`, trimmed):

```python
def systemEquation(t, Cq, qi, qiDot):
    # --- augmented matrix  [ M  Cqᵀ ; Cq  0 ] ---
    massAugmented[0:massSize, 0:massSize]      = mass_Matrix       # M    (top-left)
    massAugmented[0:massSize, massSize:matDim] = np.transpose(Cq)  # Cqᵀ  (top-right)
    massAugmented[massSize:matDim, 0:massSize] = Cq                # Cq   (bottom-left)
    #                                            bottom-right stays 0

    # --- right-hand side  [ Qe ; Qd ] ---
    Qe = np.zeros((massSize, 1))
    Qe[l2i(1,"y")] = -mass1*gravity                    # gravity → Ry slots
    Qe[l2i(2,"y")] = -mass2*gravity
    Qe[l2i(3,"y")] = -mass3*gravity
    Qe[l2i(2,"theta")] = QSpring2B + QSpring2C + QDamp2B + QDamp2C  # middle link: BOTH joints

    Qd1 = conMod.QdCalc1(qi, qiDot, u_bar_1A, 1)               # joint A
    Qd2 = conMod.QdCalc2(qi, qiDot, u_bar_1B, u_bar_2B, 1, 2)  # joint B
    Qd3 = conMod.QdCalc2(qi, qiDot, u_bar_2C, u_bar_3C, 2, 3)  # joint C
    Qd  = np.concatenate((-Qd1, Qd2, Qd3))                     # quadratic-velocity vector

    QeAug = np.concatenate((Qe, Qd))                           # [ Qe ; Qd ]  (15×1)
    qiDotDot_lamda = np.dot(np.linalg.inv(massAugmented), QeAug)  # solve → [ q̈ ; λ ]
    return qiDotDot_lamda
```

The `Q_d` centripetal helper itself:

```python
def QdCalc1(qi, qiDot, u_bar_iP, i):                  # constraintModuleTP.py (joint A)
    id = link2index(i, "theta")
    return np.square(float(qiDot[id])) * np.dot(A_i(qi[id]), u_bar_iP)   # θ̇² · A(θ)·ū
```

---

### (d) Store

`q_allTime ← qi`, `v_allTime ← qiDot`, `a_allTime ← qiDotDot_lamda[0:9]`,
`FReact_allTime ← qiDotDot_lamda[9:15]`. The full instantaneous answer
(position, velocity, acceleration, reactions) is now recorded.

---

### (e) Integrate — RK4

**Goal.** Advance the independent state `y = [θ, θ̇]` from `t` to `t+Δt`.

**Derivation.** The 3 angles obey a first-order ODE in the stacked state
`y = [θ, θ̇]`:

```math
\dot y = f(y) = \begin{bmatrix} \dot\theta \\ \ddot\theta(y) \end{bmatrix}
```

where `θ̇` is just the lower half of `y`, and `θ̈` comes from the augmented solve
(c) evaluated at `y`. Classical 4th-order Runge–Kutta advances it with a weighted
average of four slope evaluations:

```math
k_1 = f(y),\;\; k_2 = f\!\big(y+\tfrac{\Delta t}{2}k_1\big),\;\;
k_3 = f\!\big(y+\tfrac{\Delta t}{2}k_2\big),\;\; k_4 = f\!\big(y+\Delta t\,k_3\big)
```
```math
y(t+\Delta t) = y + \frac{\Delta t}{6}\big(k_1 + 2k_2 + 2k_3 + k_4\big)
```

Because each `f` evaluation needs `θ̈`, **RK4 calls the augmented solve (c) four
times per step** (re-forming `C_q` at each trial state) — separate from the one
solve in (c) that gets *recorded*.

**Code** (`main_tp.py`, first stage shown; k2/k3/k4 repeat at the trial points):

```python
def rungeKutta4_AtTimeNow(qi, qiDot, systemFunction, stepSize, timeNow):
    x    = [θ₁, θ₂, θ₃]                       # independent positions  (from qi)
    xDot = [θ̇₁, θ̇₂, θ̇₃]                      # independent velocities (from qiDot)
    y    = np.concatenate((x, xDot))          # the 6-vector state

    Cq, _, _, _ = config(qi)                   # rebuild C_q at this state
    f_1 = systemFunction(t1, Cq, qi, qiDot)    # augmented solve → accelerations
    for x in range(numberOfDOF):
        k1[x]             = y[x + numberOfDOF] # dθ/dt = θ̇  (top half = f(y))
        k1[x+numberOfDOF] = f_1[l2i(x+1,"theta")]   # dθ̇/dt = θ̈ (from the solve)
    # ... k2, k3, k4: write trial θ,θ̇ into qi/qiDot, config(), systemFunction() again ...

    yNew = y + stepSize*(k1 + 2*k2 + 2*k3 + k4)/6   # weighted RK4 step
    # write yNew back into the θ, θ̇ slots of qi, qiDot
    return qi, qiDot
```

Only the angle/rate slots of `qi`/`qiDot` are updated; the Cartesian slots get
rebuilt by (a) and (b) at the start of the next step.

---

### The whole loop at a glance

```
known: θ, θ̇  ──(a)──► dependent positions  ──(b)──► dependent velocities
        │                                                      │
        └──────────────(c) augmented solve ◄───────────────────┘
                          │  → accelerations q̈ + reactions λ
                          ▼
              (d) store  →  (e) RK4 advances θ, θ̇ to t+Δt  → repeat
```

**The insight worth keeping.** For *this* model, `M`, `C_q`, `Q_e`, `Q_d` depend
**only on the angles and angular velocities** — never on the Cartesian
coordinates (`M` is constant because reference points are at the COM; `C_q` is
`A_θ(θ)·ū`; forces are gravity + torsional terms). So the 3 independent angular
DOF form a **self-contained ODE**: RK4 marching `(θ, θ̇)` is fully correct on its
own. Steps (a) and (b) exist mainly to **reconstruct the Cartesian trajectories
for output and reaction forces**, not to drive the integration. (Add a
*translational* spring between two points and that breaks — Cartesian positions
would then feed back into `Q_e`.)

So the "we get position, velocity, acceleration" memory is the **output** view;
the mechanism is: *integrate the few real DOF (angles), reconstruct the rest from
geometry, and solve the augmented system for accelerations and joint forces at
each instant.*

---

## Step 11 — Map every symbol to code (the transcription table)

This is the final paper-to-Python step: each derived object becomes one
function or array.

| Paper object (this file) | Code |
|---|---|
| `q` layout, `link→index` map | `link2index` (`calcModuleTP.py`) |
| `A(θ)`, `A_θ(θ)` | `ATransformMatrix`, `ATransformMatrixTHETA` |
| `r = R + Aū` (Step 3) | `local2global` |
| `ū` joint points (Step 2) | `u_bar_*` globals (`main_tp.py:56–60`) |
| Constraints `C` (Step 4) | `constraintEquation`, `revolutJoint` |
| Jacobian `C_q`, partitions (Step 5) | `jacobianMatrix` |
| Newton–Raphson position (Step 10.1) | `positionAnalysis` |
| Velocity solve (Step 10.2) | `Cdi` block in `mainProg` |
| `Q_d` (Step 6) | `QdCalc1`, `QdCalc2` |
| `M`, `Jᵢ` (Step 7) | `massMatrix`, `inertiaRod` |
| `Q_e` gravity/spring/damper (Step 8) | `systemEquation`, `forceModule` |
| Augmented 15×15 solve (Step 9) | `systemEquation` |
| Reactions `−C_qᵀλ` | `FReact_allTime` in `mainProg` |
| RK4 (Step 10.4) | `rungeKutta4_AtTimeNow` |

---

## The method in one breath

> Sketch the system and count DOF → pick coordinates → locate body-fixed points
> → turn each joint into algebraic constraints → differentiate them by hand for
> the Jacobian → differentiate twice for `Q_d` → write down `M` and `Q_e` →
> stack everything into the augmented matrix → decide the per-step solve order →
> and only then translate each symbol into a function. The code is the *last*
> step, not the first — exactly the order you did it in.
