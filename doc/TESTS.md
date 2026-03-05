# Tests

## Overview

The test suite consists of four programs in `bin/`. <br>
Run them interactively *via*:

```bash
make test
```

This presents a menu to run a single test or all of them at once. <br>
Individual executables can also be invoked directly (see below).

---

## test_demo

**Purpose:** End-to-end demonstration of the full FDPACK workflow <br>
on a bivariate analytic test function.

```bash
bin/test_demo.x
```

**Test function:**

```
f(x, y) = sin(a*x) * exp(-b*y^2),   a=3, b=2
```

evaluated at `x = y = 0.5`.

**Steps performed:**

1. Noise estimation (`get_noise_table`)
2. Offset selection — three configurations:
   - central differences, `approx=.false.`
   - central differences, `approx=.true.`
   - forward differences, `approx=.false.`
3. Numerical gradient (`get_num_grad`, central)
4. Numerical Hessian (`get_num_hess`, central)
5. Spectral decomposition (`diag_num_hess`) — <br>
   only if built with `USE_LAPACK=1`

**Output explanation:**

*Step 1* prints the estimated noise level per variable and the resulting <br>
min/max offset bounds. <br>
For a smooth, noise-free function, the noise floor is at machine epsilon <br>
(~2.2e-16).

*Step 2* prints the function-call count and the optimal offsets for each <br>
gradient and Hessian entry. <br>
In `approx=.false.` mode (17 calls), full Richardson extrapolation is <br>
used for all entries. <br>
In `approx=.true.` mode (9 calls), off-diagonal Hessian offsets are <br>
estimated from cached gradient data, cutting the cost roughly in half. <br>
In forward-difference mode (5 calls), Hessian offsets are inherited <br>
directly from the gradient offsets — no separate optimization <br>
is performed.

*Steps 3–4* print a three-column table:

```
            numerical      analytic   rel. err
  g(1)=    1.2871E-01    1.2871E-01    1.2E-10
```

Relative errors at the `1e-10` to `1e-7` level are expected for central <br>
differences on a smooth function. <br>
The `H(2,2)` entry is analytically zero at this point; it is printed as <br>
`n/a` and the numerical value should be near machine epsilon.

*Step 5* (LAPACK only) prints eigenvalues, eigenvectors, and a <br>
definiteness label (positive definite / negative definite / indefinite).

**Expected result:** All relative errors below `1e-6`.

---

## test_func_noise

**Purpose:** Stress-test of `get_noise_table` across a wide range of <br>
function types and noise levels.

```bash
bin/test_func_noise.x
```

**Setup:** 150 runs, involving 8 univariate function classes with <br>
randomly drawn parameters and noise levels between ~1e-15 and ~1e-2. <br>
Artificial Gaussian noise with amplitude `p_nl` is added to each <br>
evaluation.

| Class      | Description              |
|------------|--------------------------|
| `g_osc`    | oscillatory              |
| `g_peak`   | sharp peak (Lorentzian)  |
| `g_gauss`  | Gaussian                 |
| `g_corner` | corner (exponential)     |
| `g_cont`   | continuous piecewise exp |
| `g_disc`   | discontinuous step       |
| `rosen`    | Rosenbrock               |
| `sing`     | singular (power-law)     |

**Definition of output columns:**

- Column 1 [`ID`]: <br>
  ◦ run index (1–150)
- Column 2 [`func`]: <br>
  ◦ function class name
- Column 3 [`nu(true)`]: <br>
  ◦ true noise level: `max(p_nl, eps * max(1, |f|))`
- Column 4 [`nu(est)`]: <br>
  ◦ noise level estimated by `get_noise_table`
- Column 5 [`ratio`]: <br>
  ◦ `max(nu(est)/nu(true), nu(true)/nu(est))` <br>
  ◦ 1.0 is perfect, values up to ~10 are acceptable
- Column 6 [`#func`]: <br>
  ◦ number of function evaluations used by `get_noise_table`

**Output explanation:**

A value of `ratio` close to 1 means an accurate estimate. <br>
Outliers with large `ratio` (e.g. ID 123) are edge cases where the <br> 
true noise is indistinguishable from roundoff; the estimate is not <br> 
meaningful in such cases.

Because parameters and noise levels are drawn randomly, the summary <br>
statistics vary between runs. <br>

*Typical values:*

- median `ratio`: 1–5
- max `ratio`: up to a few hundred (edge cases where the noise floor <br>
  is indistinguishable from roundoff)
- average function-call count: ~5–6 per variable
- maximum function-call count: ≤ 7 per variable

**Expected result:**

- median `ratio` around 2
- maximum below ~500
- average cost ≤ 7 evaluations per variable

---

## test_grad_offsets

**Purpose:** Stress-test of `get_fdoffsets` and `get_num_grad` across <br>
50 univariate functions, for both central and forward difference <br>
schemes.

```bash
bin/test_grad_offsets.x
```

**Setup:** 50 functions with randomly drawn parameters, evaluated at <br>
`x = 0.5`. <br>
The maximum offset is set adaptively per scheme:

- central:  `hmax = 1e5 * eps^(1/3) * max(1, |x|)`
- forward:  `hmax = 1e5 * sqrt(eps)  * max(1, |x|)`

**Definition of output columns:**

- Column 1 [`ID`]: <br>
  ◦ run index (1–100)
- Column 2 [`scheme`]: <br>
  ◦ FD scheme (`central`/`forward`)
- Column 3 [`func`]: <br>
  ◦ function ID (`f_01`–`f_50`)
- Column 4 [`ratio(r)`]: <br>
  ◦ selected offset / geometric-mean starting offset <br>
  ◦ indicator of how much the optimizer moved from the initial guess
- Column 5 [`g_ana`]: <br>
  ◦ analytic gradient value
- Column 6 [`omega`]: <br>
  ◦ optimal offset selected by `get_fdoffsets`
- Column 7 [`err`]: <br>
  ◦ relative error: `|g_num - g_ana| / |g_ana|` <br>
    (0.0 if `|g_ana| < eps`)
- Column 8 [`#func`]: <br>
  ◦ total evaluations: `get_noise_table` + `get_fdoffsets`

**Output explanation:**

Central differences achieve relative errors in the `1e-14` to `1e-9` <br>
range for well-behaved functions. <br>
Forward differences are ~5 orders of magnitude less accurate, as <br>
expected from the lower-order scheme. <br>
Outliers (*e.g.* ID 26: near-zero gradient, ID 97: function with <br>
near-singular behavior) drive the maximum error up, but are not <br>
representative of typical use. <br>
Cost is 3 evaluations (forward) or 5 (central).

Since parameters are drawn randomly, the summary statistics vary <br>
between runs. <br>

*Typical values:*

- median relative error: ~1e-10 (central) and ~1e-6 (forward)
- max relative error: up to ~1e-3 for near-singular functions <br>
  or near-zero gradients
- average function-call count: 4 per run
- maximum function-call count: 5

**Expected result:**

- central case: median relative error below `1e-8`
- forward case: median relative error below `1e-5`

---

## test_hess_offsets

**Purpose:** Stress-test of `get_fdoffsets` and `get_num_hess` across <br>
50 bivariate functions, for central differences in exact and approximate <br>
modes.

```bash
bin/test_hess_offsets.x           # exact off-diagonal offsets
bin/test_hess_offsets.x approx    # approximate off-diagonal offsets
```

**Setup:** 50 two-variable functions with randomly drawn parameters, <br>
evaluated at `x = y = 0.5`. <br>
Both central and forward differences are tested in the default run; <br>
only central differences are tested in `approx` mode. <br>
Three Hessian entries are tracked: `H(1,1)`, `H(2,1)`, and `H(2,2)`.

**Definition of output columns:**

- Column 1 [`ID`]: <br>
  ◦ 1–100 in default mode <br>
  ◦ 1–50 in `approx` mode
- Column 2 [`scheme`]: <br>
  ◦ FD scheme (`central`/`forward`)
- Column 3 [`func`]: <br>
  ◦ function ID (`f_01`–`f_50`)
- Columns 4–6 [`r(1,1)`, `r(2,1)`, `r(2,2)`]: <br>
  ◦ offset ratio (selected / geometric-mean starting offset) <br>
  ◦ for each Hessian entry
- Columns 7–9 [`h_ana(1,1)`, `h_ana(2,1)`, `h_ana(2,2)`]: <br>
  ◦ analytic Hessian values
- Columns 10–12 [`omega(1,1)`, `omega(2,1)`, `omega(2,2)`]: <br>
  ◦ optimal offsets selected by `get_fdoffsets` <br>
  ◦ for off-diagonal entries this is the second <br>
  ◦ (column-direction) offset
- Columns 13–15 [`err(1,1)`, `err(2,1)`, `err(2,2)`]: <br>
  ◦ relative errors of the Hessian entries: <br>
    `|h_num - h_ana| / |h_ana|` (0.0 if `|h_ana| < eps`) <br>
- Column 16 [`#func`]: <br>
  ◦ total evaluations: `get_noise_table` + `get_fdoffsets`

**Output explanation:**

In the default central-difference case, 17 function evaluations are <br>
used per function (full Richardson extrapolation for all Hessian <br>
entries). <br>
Forward Hessians cost only 5 evaluations but are ~4–5 orders of <br>
magnitude less accurate. <br>
In `approx` mode the cost drops to 9 evaluations; accuracy is slightly <br>
worse than in the full central mode, but well within `1e-5` for most <br>
cases.

Because parameters are drawn randomly, the summary statistics vary <br>
between runs. <br>

*Typical values:*

~ Default mode:

- median relative error: ~1e-7 (central) and ~1e-5 (forward)
- maximum relative error: up to ~1e-2 for challenging entries
- average function-call count: ~11 per run (central) and ~5 (forward)
- maximum function-call count: 17

~ `approx` mode:

- median relative error: comparable to full central mode
- maximum relative error: ~1e-5
- function-call count: 9 (fixed)

**Expected result:**

- central case: median relative error below `1e-6`
- forward case: median relative error below `1e-4`
- `approx` mode: median comparable to full central mode <br>
  at ~60% of the cost
