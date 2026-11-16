# Exercises

Exercise sheets and supporting files for *Numerical Methods for Elliptic Equations*, winter semester 2026/27.

## Exercise sheet 1

[Exercise sheet 1 (PDF)](ex01/ueb1.pdf)

The sheet also includes a voluntary local NGSolve installation task. Its supporting files are:

- [C++ extension test notebook](ex01/CppExtensionTest.ipynb)
- [C++ source: mymodule.cpp](ex01/mymodule.cpp)

Download both files into the same directory. Follow the [installation guide (PDF)](../getting-started/ngsolve-installation.pdf), activate your local environment, and start JupyterLab from that directory.

## Exercise sheet 2

[Exercise sheet 2 (PDF)](ex02/ueb2.pdf)

The numerical exercises use these starter notebooks:

- [Poisson problem on a Pacman domain](ex02/pacman_domain.ipynb)
- [Heat conduction with a material inclusion](ex02/heat_conduction.ipynb)

Work through the Pacman notebook first as the heat-conduction notebook builds on its explanations.

## Exercise sheet 3

[Exercise sheet 3 (PDF)](ex03/ueb3.pdf)

The numerical exercise uses this starter notebook:

- [Exterior electrostatics](ex03/exterior_electrostatics.ipynb)

## Exercise sheet 4

[Exercise sheet 4 (PDF)](ex04/ueb4.pdf)

The programming exercises use these starter notebooks and C++ files:

- [Anisotropic diffusion](ex04/tensor_diffusion.ipynb)
- [Cubic triangles and edge orientation](ex04/cubic_triangles.ipynb)
- [C++ starter files](ex04/implementation_ueb4/)

These exercises require a local NGSolve/Netgen installation with development headers and a compatible C++ compiler. Follow [Getting started](../getting-started/README.md).

Keep the repository layout when downloading: both notebooks use [teaching_utils.py](../notebooks/ch06_implementation/teaching_utils.py) from the [Chapter 6 material](../notebooks/ch06_implementation/README.md). Download or clone the repository, keep all files in `ex04/implementation_ueb4/` together, and start JupyterLab from `content/exercises/ex04/`. Run each notebook in its own fresh kernel and restart the kernel after editing C++ headers.

[Course overview](../index.ipynb) · [Lecture notes](../script/README.md) · [Lecture notebooks](../notebooks/README.md)
