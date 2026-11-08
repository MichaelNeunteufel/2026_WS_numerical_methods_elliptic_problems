#ifndef TEACHING_MY_HIGH_ORDER_FESPACE_HPP
#define TEACHING_MY_HIGH_ORDER_FESPACE_HPP

#include <comp.hpp>
#include <cmath>
#include "myHighOrderElement.hpp"
namespace teaching
{
  using namespace ngcomp;

  inline int HighOrderDegree(const Flags &flags)
  {
    double degree = flags.GetNumFlag("order", 1);
    if (!std::isfinite(degree) || degree < 1 || degree != std::floor(degree))
      throw Exception("The polynomial degree must be a positive integer.");
    return int(degree);
  }

  class HighOrderFESpace : public FESpace
  {
    int p;
    int nv = 0, edge_dofs = 0;
    Array<int> first_cell_dof;

  public:
    HighOrderFESpace(shared_ptr<MeshAccess> mesh, const Flags &flags)
        : FESpace(mesh, flags), p(HighOrderDegree(flags))
    {
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
    void Update() override
    {
      FESpace::Update();
      for (auto active : definedon[VOL])
        if (!active)
          throw Exception("The teaching space is defined on the full mesh.");
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

      first_cell_dof.Last() = next;
      SetNDof(next);
      ctofdof.SetSize(next);
      ctofdof = INTERFACE_DOF;
      for (int cell = 0; cell < ma->GetNE(VOL); cell++)
        for (int d = first_cell_dof[cell]; d < first_cell_dof[cell + 1]; d++)
          ctofdof[d] = LOCAL_DOF;
    }
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
      fe->SetVertexNumbers(vertices);
      return *fe;
    }
    void GetVertexDofNrs(int vertex, Array<DofId> &dnums) const override
    {
      dnums.SetSize(1);
      dnums[0] = vertex;
    }
    void GetEdgeDofNrs(int edge, Array<DofId> &dnums) const override
    {
      dnums.SetSize(edge_dofs);
      for (int k = 0; k < edge_dofs; k++)
        dnums[k] = nv + edge_dofs * edge + k;
    }
    void GetFaceDofNrs(int, Array<DofId> &dnums) const override { dnums.SetSize0(); }
    void GetInnerDofNrs(int cell, Array<DofId> &dnums) const override
    {
      dnums.SetSize0();
      for (int d = first_cell_dof[cell]; d < first_cell_dof[cell + 1]; d++)
        dnums.Append(d);
    }
  };
  static RegisterFESpace<HighOrderFESpace> register_high_order_space("teachinghighorder");
} // namespace teaching

#endif // TEACHING_MY_HIGH_ORDER_FESPACE_HPP
