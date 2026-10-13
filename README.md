# PDE Homework 2 — Numerical Solution of the 2D Stokes Equations

Numerical solution of the two-dimensional Stokes equations on a staggered (MAC) grid, solved as a sparse saddle-point linear system in MATLAB. Accuracy is verified against a manufactured solution, demonstrating second-order convergence for velocity and pressure. An additional experiment models flow through an idealized 2D porous microstructure (circular obstacles) using Brinkman penalization.

## Accessing the Report

**Option 1 — Pre-built PDF (quickest):**

The compiled report is included in this repository as [`report.pdf`](report.pdf). View it directly on GitHub by clicking the file, or download it:

```bash
curl -O https://github.com/lhenriquez02/PDE_HW2/raw/main/report.pdf
```

**Option 2 — Rebuild from source:**

```bash
git clone https://github.com/lhenriquez02/PDE_HW2.git
cd PDE_HW2
latexmk -pdf report.tex   # or: pdflatex report.tex
```

This requires a LaTeX distribution with `pdflatex` and the packages loaded in the preamble (`graphicx`, `float`, `amsmath`, `hyperref`, etc.).

## Repository Contents

| File | Description |
|---|---|
| `report.tex` | LaTeX source of the report |
| `report.pdf` | Compiled report |
| `code.m` | MATLAB solver: MAC-grid Stokes solve, manufactured-solution error study, and generation of the velocity and error figures |
| `porousBrinkman.m` | MATLAB solver for the additional experiment: flow through a porous microstructure via Brinkman penalization |
| `stokes_velocity.pdf` | Velocity magnitude plot with streamlines (Figure 1 in the report) |
| `stokes_errors.pdf` | Relative L² error vs. grid resolution, with first/second-order reference slopes (Figure 2 in the report) |
| `porous_velocity.pdf` | Velocity magnitude and streamlines for the porous medium experiment (Figure 3 in the report) |

## Running the Code

Both scripts are self-contained MATLAB files (no toolboxes required):

```matlab
code.m              % Stokes solve + convergence study, saves stokes_*.pdf
porousBrinkman.m    % Brinkman-penalized porous flow, saves porous_velocity.pdf
```

Each script saves its figures as vector PDFs next to the script file. `code.m` also prints the observed convergence orders for `U`, `V`, and `P` to the console.

## Problem Summary

- **Governing equations:** μΔ**u** − ∇P = **b** with ∇·**u** = 0 on the unit square.
- **Discretization:** staggered MAC grid — pressure at cell centers, `U` on vertical faces, `V` on horizontal faces — assembled into a sparse saddle-point system and solved with MATLAB's backslash. The pressure null space is handled with a one-point constraint.
- **Boundary conditions (main problem):** periodic in x; no-slip top/bottom for `U` (U = 0); prescribed vertical wall velocity V = −3.5.
- **Verification:** manufactured solutions with relative L² error norms showing ≈ O(N⁻²) convergence for all variables.
- **Additional experiment:** solid grains are approximated with Brinkman penalization (large drag coefficient αₛ in solid regions), driven by a constant vertical body force with fully periodic boundaries.
