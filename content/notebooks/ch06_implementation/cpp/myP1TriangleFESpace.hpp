#ifndef TEACHING_MY_P1_TRIANGLE_FESPACE_HPP
#define TEACHING_MY_P1_TRIANGLE_FESPACE_HPP

#include <comp.hpp>
#include "myP1TriangleElement.hpp"
#include "myP1SegmentElement.hpp"

namespace teaching
{
  using namespace ngcomp;

  // Continuous piecewise affine functions on a two-dimensional triangular mesh
  class P1TriangleFESpace : public FESpace
  {
  public:
    P1TriangleFESpace(shared_ptr<MeshAccess> mesh, const Flags &flags)
        : FESpace(mesh, flags)
    {
      if (ma->GetDimension() != 2)
        throw Exception("P1TriangleFESpace requires a two-dimensional triangular mesh.");
      if (IsComplex() || GetDimension() != 1)
        throw Exception("P1TriangleFESpace is real and scalar.");
      if (flags.GetNumFlag("order", 1) != 1)
        throw Exception("P1TriangleFESpace has fixed polynomial degree one.");
      type = "teachingp1triangle";

      evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<2>>>();
      flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<2>>>();
      evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundary<2>>>();
    }

    string GetClassName() const override { return "P1TriangleFESpace"; }
    static DocInfo GetDocu()
    {
      auto doc = FESpace::GetDocu();
      doc.short_docu = "P1 triangles with one nodal coefficient per mesh vertex.";
      return doc;
    }

    void Update() override
    {
      FESpace::Update();
      for (auto active : definedon[VOL])
        if (!active)
          throw Exception("P1TriangleFESpace is defined on the full mesh.");
      for (auto ei : ma->Elements(VOL))
        if (ma->GetElement(ei).GetType() != ET_TRIG)
          throw Exception("P1TriangleFESpace supports triangular volume elements only.");

      // A global degree of freedom has the same number as its mesh vertex
      SetNDof(ma->GetNV());
      ctofdof.SetSize(GetNDof());
      ctofdof = INTERFACE_DOF;
    }

    void GetDofNrs(ElementId ei, Array<DofId> &dnums) const override
    {
      dnums.SetSize0();
      for (auto vertex : ma->GetElement(ei).Vertices())
        dnums.Append(vertex);
    }

    FiniteElement &GetFE(ElementId ei, Allocator &alloc) const override
    {
      switch (ma->GetElement(ei).GetType())
      {
      case ET_TRIG:
        return *new (alloc) P1TriangleElement;
      case ET_SEGM:
        return *new (alloc) P1SegmentElement;
      default:
        throw Exception("P1TriangleFESpace supplies triangles and boundary segments.");
      }
    }

    // Used by FESpace::UpdateFreeDofs through vertex-node DOF queries
    void GetVertexDofNrs(int vertex, Array<DofId> &dnums) const override
    {
      dnums.SetSize(1);
      dnums[0] = vertex;
    }

    // P1 has no edge DOFs; extend for P2/P3 triangles and Q2 quadrilaterals
    void GetEdgeDofNrs(int edge, Array<DofId> &dnums) const override
    {
      dnums.SetSize0();
    }

    // Face DOFs are used in 3D; use GetInnerDofNrs for 2D cell interiors
    void GetFaceDofNrs(int face, Array<DofId> &dnums) const override
    {
      dnums.SetSize0();
    }

    // P1 has no interior DOFs; extend for P3 triangles and Q2 quadrilaterals
    void GetInnerDofNrs(int cell, Array<DofId> &dnums) const override
    {
      dnums.SetSize0();
    }
  };

  static RegisterFESpace<P1TriangleFESpace> register_p1_triangle_space("teachingp1triangle");
} // namespace teaching

#endif // TEACHING_MY_P1_TRIANGLE_FESPACE_HPP
