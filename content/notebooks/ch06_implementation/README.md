# Finite element implementation

Companion material for Chapter 6 of the [lecture notes](../../script/numerical_methods.pdf).

These notebooks compile a C++ extension and must be run locally. Native C++ compilation is unavailable in JupyterLite. See [Getting started](../../getting-started/README.md) for installation instructions.

Download the complete `ch06_implementation/` directory, including `teaching_utils.py` and all files in `cpp/`. Keep this directory structure when running the examples.

Run the notebooks from this directory, each in its own fresh Python kernel:

- [00_linear_triangle_space.ipynb](00_linear_triangle_space.ipynb): a complete P1 triangle and boundary segment, vertex numbering, all four element matrices and load vectors from the manuscript, explicit scatter-add compared with C++ assembly, and nonzero Dirichlet data.
- [01_nodal_elements.ipynb](01_nodal_elements.ipynb): P1/Q1 functions, nodal/hierarchical P2 coordinates, integrated Legendre edges and Dubiner bubbles, interactive WebGUI views on one triangle and a triangle/quad patch, shared-edge numbering and orientation, Pp/Qp enrichment for a smooth manufactured solution, and P1 tetrahedra.
- [02_assembly.ipynb](02_assembly.ipynb): reference-to-physical gradients, local quadrature, numerical comparison of three assembly paths, mixed boundary data, quadrature and actual errors.

The two element/space pairs are:

- `P1TriangleElement` / `P1TriangleFESpace`: an elementary triangle with explicit values and gradients, direct inheritance from NGSolve's `FESpace`, and one coefficient per mesh vertex. `P1SegmentElement` supplies the boundary trace. This path needs neither automatic differentiation nor edge orientation.
- `HighOrderElement<ET>` / `HighOrderFESpace`: hierarchical Pp triangles and Qp quadrilaterals for every uniform integer degree p>=1 in 2D, with matching segment traces. Here `ET` selects the reference topology, for example `ET_TRIG` for triangles. Shape formulas are templated for automatic differentiation. The element also derives from NGSolve's `VertexOrientedFE<ET>`. The space supplies global vertex numbers with `SetVertexNumbers`. Edge modes use the endpoints from `GetVertexOrientedEdge`, and triangle interior modes use `GetVertexOrientedFace(0)` to order their barycentric coordinates. Thus global numbering lists edge modes by increasing degree. Cell-interior modes are marked local, separately from shared vertex and edge coefficients. Three-dimensional examples remain P1 tetrahedra.

The high-order files are [myHighOrderElement.hpp](cpp/myHighOrderElement.hpp) and [myHighOrderFESpace.hpp](cpp/myHighOrderFESpace.hpp). Triangle interior modes multiply the reference Dubiner basis by the triangle bubble. Their global vertex ordering is optional for cell-local bubbles in 2D, but illustrates the convention needed when triangular modes become shared face traces of high-order tetrahedra. Quadrilateral modes are tensor products of endpoint and integrated Legendre functions. Homogeneous Legendre recurrences avoid division by collapsed barycentric coordinates.

The introductory `P1TriangleFESpace` is fixed to degree one on 2D triangular meshes. Its explicit local functions are in [myP1TriangleElement.hpp](cpp/myP1TriangleElement.hpp) and [myP1SegmentElement.hpp](cpp/myP1SegmentElement.hpp), and its vertex-only numbering is in [myP1TriangleFESpace.hpp](cpp/myP1TriangleFESpace.hpp). The notebooks list reference vertices and edge-midpoint coordinates explicitly for the basis-coordinate comparison. The finite elements returned by `fes.GetFE(...)` provide reference values and gradients through `CalcShape` and `CalcDShape`. `HighOrderShapes` exposes values and gradients of the high-order family for the basis illustrations.

The shared [myIntegrator.hpp](cpp/myIntegrator.hpp) contains one class per weak-form term: `LaplaceIntegrator<D>`, `MassIntegrator<D>`, `BoundaryMassIntegrator<D>`, `SourceIntegrator<D>` and `BoundarySourceIntegrator<D>`. Their Python factories remain `Laplace`, `Mass`, `BoundaryMass`, `Source` and `BoundarySource`. The module also supplies serial assembly. Geometry, sparse storage, bindings and direct solvers come from NGSolve.

Requirements: NGSolve/Netgen with development headers, a working C++ compiler, NumPy, SciPy, Matplotlib, ipywidgets and Jupyter.

After editing a C++ header, restart the kernel before recompilation. The runtime compiler can reuse a loaded module. A header edit alone does not reliably trigger recompilation.
