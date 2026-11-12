#ifndef TEACHING_MY_P1_TRIANGLE_FESPACE_HPP
#define TEACHING_MY_P1_TRIANGLE_FESPACE_HPP

#include <comp.hpp>
#include "myP1TriangleElement.hpp"
#include "myP1SegmentElement.hpp"

namespace teaching
{
  using namespace ngcomp;

  // Continuous piecewise affine functions on a two-dimensional triangular mesh
  // FESpace connects local reference elements to global DOFs and symbolic operators.
  class P1TriangleFESpace : public FESpace
  {
  public:
    // MeshAccess supplies mesh topology/geometry. Flags contains options such as dirichlet.
    P1TriangleFESpace(shared_ptr<MeshAccess> mesh, const Flags &flags)
        : FESpace(mesh, flags)
    {
      // ma is the mesh; FESpace::GetDimension() below counts field components (one for scalars).
      if (ma->GetDimension() != 2)
        throw Exception("P1TriangleFESpace requires a two-dimensional triangular mesh.");
      if (IsComplex() || GetDimension() != 1)
        throw Exception("P1TriangleFESpace is real and scalar.");
      if (flags.GetNumFlag("order", 1) != 1)
        throw Exception("P1TriangleFESpace has fixed polynomial degree one.");
      type = "teachingp1triangle";

      // Register value and gradient evaluators for symbolic forms. <2> is the spatial dimension.
      // VOL acts on cells. The BND identity evaluator supplies the boundary trace.
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

    // Rebuild global numbering when the space is updated, for example after refinement.
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
      // INTERFACE_DOF: global/shared vertex coefficients remain in the condensed system.
      // Sharing comes from GetDofNrs. Dirichlet constraints use the separate GetFreeDofs() mask.
      ctofdof = INTERFACE_DOF;
    }

    // Return global DOF numbers in the same order as the local shape functions.
    // Neighboring elements reuse vertex numbers, so they share the corresponding coefficients.
    void GetDofNrs(ElementId ei, Array<DofId> &dnums) const override
    {
      dnums.SetSize0();
      for (auto vertex : ma->GetElement(ei).Vertices())
        dnums.Append(vertex);
    }

    // ElementId identifies a volume or boundary element and its mesh number.
    // Placement new uses the supplied allocator; the returned reference has that lifetime.
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

    // Node queries identify DOFs on mesh entities, including Dirichlet boundary vertices.
    // NGSolve uses them to build the free-DOF mask (Python: fes.FreeDofs()).
    void GetVertexDofNrs(int vertex, Array<DofId> &dnums) const override
    {
      dnums.SetSize(1);
      dnums[0] = vertex;
    }

    // P1 has no edge DOFs. Extend for P2/P3 triangles and Q2 quadrilaterals
    void GetEdgeDofNrs(int edge, Array<DofId> &dnums) const override
    {
      dnums.SetSize0();
    }

    // Face DOFs are used in 3D. Use GetInnerDofNrs for 2D cell interiors
    void GetFaceDofNrs(int face, Array<DofId> &dnums) const override
    {
      dnums.SetSize0();
    }

    // P1 has no interior DOFs. Extend for P3 triangles and Q2 quadrilaterals
    void GetInnerDofNrs(int cell, Array<DofId> &dnums) const override
    {
      dnums.SetSize0();
    }
  };

  // Register this class under the name used by the NGSolve FESpace factory.
  static RegisterFESpace<P1TriangleFESpace> register_p1_triangle_space("teachingp1triangle");
} // namespace teaching

#endif // TEACHING_MY_P1_TRIANGLE_FESPACE_HPP
