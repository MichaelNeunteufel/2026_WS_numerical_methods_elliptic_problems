# Lecture notebooks

Interactive examples for *Numerical Methods for Elliptic Problems*. Use the [course JupyterLite site](https://michaelneunteufel.github.io/2026_WS_numerical_methods_elliptic_problems/lab/index.html?path=index.ipynb) or a local JupyterLab installation, see [Getting started](../getting-started/README.md).

Run code cells in order. Before working through the examples, use the [NGSolve and WebGUI check](00_ngsolve_test.ipynb) to verify that a triangular mesh is displayed.

## Lecture 1 — October 7, 2026

[Motivation examples from Chapter 1](VO_2026_10_07/index.ipynb):

- [Mass and heat transfer](VO_2026_10_07/01_mass_and_heat_transfer.ipynb)
- [Minimal surfaces and membranes](VO_2026_10_07/02_minimal_surface_and_membrane.ipynb)
- [The biharmonic plate](VO_2026_10_07/03_biharmonic_plate.ipynb)
- [Linear elasticity](VO_2026_10_07/04_linear_elasticity.ipynb)
- [Incompressible flow](VO_2026_10_07/05_incompressible_flow.ipynb)
- [Magneto-quasistatics](VO_2026_10_07/06_magneto_quasistatics.ipynb)

## Lecture 3 — October 14, 2026

[Laplace problem discretization with linear finite elements](VO_2026_10_14/01_simple_laplace.ipynb)

## Lecture 4 — October 21, 2026

[Friedrichs, Poincaré, and trace constants](VO_2026_10_21/01_sobolev_constants.ipynb)

Experiments with shrinking Dirichlet boundary parts, volume and boundary measurements, a single fixed finite element vertex, and the trace constant on thin rectangles.

## October 29, 2026 — Boundary data

[Boundary data: lifting, mixed conditions, and the Neumann kernel](VO_2026_10_29/01_boundary_data.ipynb)

Examples from Chapter 4: nonhomogeneous Dirichlet data and liftings, combined Dirichlet–Neumann–Robin conditions, and compatibility and mean constraints for pure Neumann problems. The convection examples compare a coarse Galerkin breakdown with a divergence-free field and show how strong convection creates an outflow boundary layer.

## November 5, 2026 — Finite element approximation and error estimates

- [Interpolation and Galerkin approximation](VO_2026_11_05/01_interpolation_and_galerkin.ipynb): nodal interpolation, best approximation in energy, and first- and second-order error decay.
- [Why does shape regularity matter?](VO_2026_11_05/02_shape_regularity.ipynb): interpolation on a flattening triangle and the different geometric assumptions in the energy and L2 estimates.
- [Convergence near a reentrant corner](VO_2026_11_05/03_corner_singularity.ipynb): an exact singular solution on an L-shaped domain and the roles of primal and dual regularity.
- [Inverse estimates and maximum-norm errors](VO_2026_11_05/04_inverse_estimates.ipynb): inverse norm ratios for a shrinking nodal hat and logarithmic maximum-error growth on two uniform triangulations.

Each notebook is self-contained and can be run in a fresh kernel.

The examples accompany the [lecture notes](../script/README.md). Related tasks are listed under [Exercises](../exercises/README.md).

[Course overview](../index.ipynb)
