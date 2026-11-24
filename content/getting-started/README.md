# Getting started

These materials accompany *Numerical Methods for Elliptic Equations*, winter semester 2026/27 at Johannes Kepler University Linz.

## Work locally

Use your local NGSolve installation for the course notebooks and programming exercises. The [Chapter 6 implementation notebooks](../notebooks/ch06_implementation/README.md) and the programming tasks on [Exercise sheet 4](../exercises/README.md#exercise-sheet-4) compile native C++ code and require NGSolve/Netgen with development headers and a compatible C++ compiler. The C++ programming tasks are part of the regular course exercises.

The [NGSolve installation guide (PDF)](ngsolve-installation.pdf) explains installation and checks for Python, WebGUI, and the C++ compiler. Run the [C++ extension test](../exercises/ex01/CppExtensionTest.ipynb) to check your compiler setup.

**Windows.** Visual Studio 2022 Community, Professional, Enterprise, and the standalone Build Tools all support the course exercises with the *Desktop development with C++* workload installed. Start JupyterLab from your installation’s *x64 Native Tools Command Prompt* and run the extension test with the course Python kernel.

For Chapter 6, download the complete `ch06_implementation/` directory with `teaching_utils.py` and `cpp/`, and start JupyterLab from that directory. Each notebook should use its own fresh Python kernel.

For Exercise sheet 4, keep `content/exercises/ex04/` and `content/notebooks/ch06_implementation/` in their repository locations, including the C++ starter directory and `teaching_utils.py`. Start JupyterLab from `content/exercises/ex04/` and restart the kernel after changing C++ headers.

Start with the [NGSolve and WebGUI check](../notebooks/00_ngsolve_test.ipynb). Then explore the [lecture notebooks](../notebooks/README.md) and [exercise sheets](../exercises/README.md).

## Work in the browser

The [course JupyterLite site](https://michaelneunteufel.github.io/2026_WS_numerical_methods_elliptic_problems/lab/index.html?path=index.ipynb) provides browser access to many lecture examples and exercise notebooks. Tasks that compile C++ extensions must be run locally.

Download your work regularly. Browser files are stored in the current browser profile and are not synchronized between devices.

[Course overview](../index.ipynb)
