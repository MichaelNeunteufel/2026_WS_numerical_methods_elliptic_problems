#ifndef TEACHING_MY_HIGH_ORDER_FESPACE_HPP
#define TEACHING_MY_HIGH_ORDER_FESPACE_HPP

#include <comp.hpp>
#include <cmath>
#include "myHighOrderElement.hpp"
namespace teaching
{
  using namespace ngcomp;

  // Read the NGSolve order option, defaulting to degree one.
  inline int HighOrderDegree(const Flags &flags)
  {
    double degree = flags.GetNumFlag("order", 1);
    if (!std::isfinite(degree) || degree < 1 || degree != std::floor(degree))
      throw Exception("The polynomial degree must be a positive integer.");
    return int(degree);
  }

  // FESpace owns global numbering and evaluators. HighOrderElement owns local basis formulas.
  class HighOrderFESpace : public FESpace
  {
    int p;
    int nv = 0, edge_dofs = 0;
    // Cell c owns the half-open DOF range [first_cell_dof[c], first_cell_dof[c+1]).
    Array<int> first_cell_dof;

  public:
    HighOrderFESpace(shared_ptr<MeshAccess> mesh, const Flags &flags)
        : FESpace(mesh, flags), p(HighOrderDegree(flags))
    {
      // GetDimension() counts field components. ma->GetDimension() is the geometric dimension.
      if (IsComplex() || GetDimension() != 1)
        throw Exception("The teaching space is real and scalar.");
      int dim = ma->GetDimension();
      if (dim != 2 && dim != 3)
        throw Exception("Use a 2D or 3D mesh.");
      if (dim == 3 && p != 1)
        throw Exception("Use positive degree in 2D; only P1 tetrahedra in 3D.");
      type = "teachinghighorder";
      if (dim == 2)
      {
        // Supply symbolic value, gradient and boundary-trace operators for this scalar space.
        // The differential-operator template argument is the physical dimension, not the degree.
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<2>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<2>>>();
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundary<2>>>();
      }
      else
      {
        evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<3>>>();
        flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<3>>>();
        evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundary<3>>>();
      }
    }
    string GetClassName() const override { return "HighOrderFESpace"; }
    static DocInfo GetDocu()
    {
      auto doc = FESpace::GetDocu();
      doc.short_docu = "Hierarchical Pp triangles/Qp quads; uniform positive degree; P1 tetrahedra.";
      return doc;
    }
    // Rebuild entity-based numbering and coupling types for the current mesh.
    void Update() override
    {
      FESpace::Update();
      for (auto active : definedon[VOL])
        if (!active)
          throw Exception("The teaching space is defined on the full mesh.");

      // Global blocks: one DOF per vertex, p-1 per edge, then the cell-interior DOFs.
      nv = ma->GetNV();
      edge_dofs = p - 1;
      int next = nv + edge_dofs * ma->GetNEdges();
      first_cell_dof.SetSize(ma->GetNE(VOL) + 1);
      for (auto ei : ma->Elements(VOL))
      {
        auto et = ma->GetElement(ei).GetType();
        HighOrderNDof(et, p);
        first_cell_dof[ei.Nr()] = next;
        if (et == ET_TRIG)
          next += (p - 1) * (p - 2) / 2;
        if (et == ET_QUAD)
          next += (p - 1) * (p - 1);
      }

      // The final offset closes the DOF range of the last cell.
      first_cell_dof.Last() = next;
      SetNDof(next);
      ctofdof.SetSize(next);
      // INTERFACE_DOF keeps global/shared vertex and edge coefficients in the interface system.
      // GetDofNrs creates the sharing. Dirichlet constraints use the separate GetFreeDofs() mask.
      ctofdof = INTERFACE_DOF;
      // LOCAL_DOF marks cell bubbles: they can be eliminated by static condensation.
      for (int cell = 0; cell < ma->GetNE(VOL); cell++)
        for (int d = first_cell_dof[cell]; d < first_cell_dof[cell + 1]; d++)
          ctofdof[d] = LOCAL_DOF;
    }
    // Map local shapes to global coefficients: vertices, edges, then cell interiors.
    // A shared edge has one global block, reused by its neighboring elements.
    void GetDofNrs(ElementId ei, Array<DofId> &dnums) const override
    {
      dnums.SetSize0();
      auto element = ma->GetElement(ei);
      auto vertices = element.Vertices();
      for (auto vertex : vertices)
        dnums.Append(vertex);
      if (edge_dofs)
      {
        auto edges = element.Edges();
        for (int e = 0; e < edges.Size(); e++)
        {
          for (int k = 0; k < edge_dofs; k++)
            dnums.Append(nv + edge_dofs * edges[e] + k);
        }
      }
      if (ei.VB() == VOL)
        for (int d = first_cell_dof[ei.Nr()]; d < first_cell_dof[ei.Nr() + 1]; d++)
          dnums.Append(d);
    }
    // Construct the local reference element in caller-provided storage, using its mesh topology.
    FiniteElement &GetFE(ElementId ei, Allocator &alloc) const override
    {
      auto et = ma->GetElement(ei).GetType();
      FiniteElement *fe;
      switch (et)
      {
      case ET_SEGM:
        fe = new (alloc) HighOrderElement<ET_SEGM>(p);
        break;
      case ET_TRIG:
        fe = new (alloc) HighOrderElement<ET_TRIG>(p);
        break;
      case ET_QUAD:
        fe = new (alloc) HighOrderElement<ET_QUAD>(p);
        break;
      case ET_TET:
        fe = new (alloc) HighOrderElement<ET_TET>(p);
        break;
      default:
        throw Exception("Point finite elements are not needed.");
      }
      Array<int> vertices;
      for (auto v : ma->GetElement(ei).Vertices())
        vertices.Append(v);
      // Global vertex IDs make edge-mode orientations agree between adjacent elements.
      fe->SetVertexNumbers(vertices);
      return *fe;
    }
    // Entity queries let NGSolve locate DOFs, for example when imposing Dirichlet data.
    void GetVertexDofNrs(int vertex, Array<DofId> &dnums) const override
    {
      dnums.SetSize(1);
      dnums[0] = vertex;
    }
    // Only edge-interior modes belong to this block. Endpoint DOFs belong to the vertices.
    void GetEdgeDofNrs(int edge, Array<DofId> &dnums) const override
    {
      dnums.SetSize(edge_dofs);
      for (int k = 0; k < edge_dofs; k++)
        dnums[k] = nv + edge_dofs * edge + k;
    }
    // The supported P1 tetrahedra have no face modes. 2D cell bubbles are queried below.
    void GetFaceDofNrs(int, Array<DofId> &dnums) const override { dnums.SetSize0(); }
    void GetInnerDofNrs(int cell, Array<DofId> &dnums) const override
    {
      dnums.SetSize0();
      for (int d = first_cell_dof[cell]; d < first_cell_dof[cell + 1]; d++)
        dnums.Append(d);
    }
  };
  // Make this space available to the NGSolve factory under the given type name.
  static RegisterFESpace<HighOrderFESpace> register_high_order_space("teachinghighorder");
} // namespace teaching

#endif // TEACHING_MY_HIGH_ORDER_FESPACE_HPP
