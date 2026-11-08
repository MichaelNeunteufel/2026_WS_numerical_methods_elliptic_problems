"""Small mesh/plot helpers for the Chapter 6 teaching notebooks."""

from pathlib import Path
from itertools import combinations
import numpy as np
from ngsolve import Mesh, CF, Integrate, grad, sqrt, x, y, z
from ngsolve.fem import CompilePythonModule
from netgen.meshing import Mesh as NetgenMesh, MeshPoint, Pnt, Element2D, Element1D


def load_module():
    return CompilePythonModule(
        Path(__file__).resolve().parent / "cpp/mymodule.cpp",
        init_function_name="teaching_fem",
    )


def square_mesh(n=4, quadrilateral=False, distorted=False):
    """A pure triangular or pure quadrilateral mesh, with four named sides."""
    if n < 1:
        raise ValueError("n must be positive")
    mesh = NetgenMesh(dim=2)
    domain = mesh.AddRegion("domain", dim=2)
    sides = {
        name: mesh.AddRegion(name, dim=1) for name in ["left", "right", "bottom", "top"]
    }
    points = {}
    for j in range(n + 1):
        for i in range(n + 1):
            xx, yy = i / n, j / n
            if distorted and 0 < i < n and 0 < j < n:
                displacement = np.sin(np.pi * i / n) * np.sin(np.pi * j / n) / n
                xx += 0.16 * displacement
                yy -= 0.11 * displacement
            points[i, j] = mesh.Add(MeshPoint(Pnt(xx, yy, 0)))
    for j in range(n):
        for i in range(n):
            bl, br = points[i, j], points[i + 1, j]
            tr, tl = points[i + 1, j + 1], points[i, j + 1]
            if quadrilateral:
                mesh.Add(Element2D(domain, [bl, br, tr, tl]))
            elif (i + j) % 2 == 0:
                mesh.Add(Element2D(domain, [bl, br, tr]))
                mesh.Add(Element2D(domain, [bl, tr, tl]))
            else:
                mesh.Add(Element2D(domain, [bl, br, tl]))
                mesh.Add(Element2D(domain, [br, tr, tl]))
    for i in range(n):
        mesh.Add(Element1D([points[i, 0], points[i + 1, 0]], index=sides["bottom"]))
        mesh.Add(Element1D([points[i + 1, n], points[i, n]], index=sides["top"]))
    for j in range(n):
        mesh.Add(Element1D([points[0, j + 1], points[0, j]], index=sides["left"]))
        mesh.Add(Element1D([points[n, j], points[n, j + 1]], index=sides["right"]))
    return Mesh(mesh)


def max_diameter(mesh):
    coordinates = np.array([vertex.point for vertex in mesh.vertices])
    return max(
        np.linalg.norm(coordinates[a.nr] - coordinates[b.nr])
        for element in mesh.Elements()
        for a, b in combinations(element.vertices, 2)
    )


def errors(mesh, uh, exact, quadrature_order=12):
    exact = CF(exact)
    difference = exact - uh
    exact_gradient = CF(tuple(exact.Diff(t) for t in (x, y, z)[: mesh.dim]))
    gradient = exact_gradient - grad(uh)
    l2_squared = Integrate(difference**2, mesh, order=quadrature_order)
    gradient_squared = Integrate(gradient * gradient, mesh, order=quadrature_order)
    return dict(L2=sqrt(l2_squared), H1=sqrt(l2_squared + gradient_squared))
