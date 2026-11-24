#ifndef TEACHING_MY_P3_TRIANGLE_FESPACE_HPP
#define TEACHING_MY_P3_TRIANGLE_FESPACE_HPP
#include <comp.hpp>
#include "myP3TriangleElement.hpp"
#include "myP3SegmentElement.hpp"

namespace teaching
{
  using namespace ngcomp;

  // Fixed-degree extension of P1TriangleFESpace; all orientation is in GetDofNrs.
  class P3TriangleFESpace : public FESpace
  {
  public:
    P3TriangleFESpace(shared_ptr<MeshAccess> mesh, const Flags &flags)
        : FESpace(mesh, flags)
    {
      if (ma->GetDimension() != 2)
        throw Exception("P3TriangleFESpace requires a two-dimensional triangular mesh.");
      if (IsComplex() || GetDimension() != 1)
        throw Exception("P3TriangleFESpace is real and scalar.");
      if (flags.GetNumFlag("order", 3) != 3)
        throw Exception("P3TriangleFESpace has fixed polynomial degree three.");
      order = 3;
      type = "teachingueb4p3triangle";
      evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpId<2>>>();
      flux_evaluator[VOL] = make_shared<T_DifferentialOperator<DiffOpGradient<2>>>();
      evaluator[BND] = make_shared<T_DifferentialOperator<DiffOpIdBoundary<2>>>();
    }
    string GetClassName() const override { return "P3TriangleFESpace"; }
    static DocInfo GetDocu()
    {
      auto doc = FESpace::GetDocu();
      doc.short_docu = "Nodal P3 triangles with explicit edge-DOF permutations.";
      return doc;
    }
    void Update() override
    {
      FESpace::Update();
      const auto &active_domains = DefinedOn(VOL);
      for (auto domain : Range(active_domains.Size()))
        if (!active_domains[domain])
          throw Exception("P3TriangleFESpace is defined on the full mesh.");
      for (auto ei : ma->Elements(VOL))
        if (ma->GetElement(ei).GetType() != ET_TRIG)
          throw Exception("P3TriangleFESpace supports triangular volume elements only.");
      // BEGIN EXERCISE_P3_NUMBERING
      // TODO: implement this part of the exercise.
      throw Exception("TODO: EXERCISE_P3_NUMBERING");
      // END EXERCISE_P3_NUMBERING
      // Supplied NGSolve bookkeeping: shared coefficients and cell interiors.
      ctofdof.SetSize(GetNDof());
      ctofdof = INTERFACE_DOF;
      int first_inner = ma->GetNV() + 2 * ma->GetNEdges();
      for (int cell = 0; cell < ma->GetNE(VOL); cell++)
        ctofdof[first_inner + cell] = LOCAL_DOF;
    }
    void GetDofNrs(ElementId ei, Array<DofId> &dnums) const override
    {
      auto element = ma->GetElement(ei);
      auto vertices = element.Vertices();
      auto edges = element.Edges();
      auto topology = ElementTopology::GetEdges(element.GetType());
      dnums.SetSize0();
      for (auto vertex : vertices)
        dnums.Append(vertex);
      for (int e = 0; e < edges.Size(); e++)
      {
        int a = topology[e][0], b = topology[e][1];
        // Global edge direction: smaller to larger global vertex number.
        // The same rule is needed for volume triangles AND boundary segments.
        // BEGIN EXERCISE_EDGE_ORIENTATION
        // TODO: implement this part of the exercise.
        throw Exception("TODO: EXERCISE_EDGE_ORIENTATION");
        // END EXERCISE_EDGE_ORIENTATION
      }
      if (ei.VB() == VOL)
      {
        Array<DofId> inner;
        GetInnerDofNrs(ei.Nr(), inner);
        for (auto dof : inner)
          dnums.Append(dof);
      }
    }
    FiniteElement &GetFE(ElementId ei, Allocator &alloc) const override
    {
      switch (ma->GetElement(ei).GetType())
      {
      case ET_TRIG:
        return *new (alloc) P3TriangleElement;
      case ET_SEGM:
        return *new (alloc) P3SegmentElement;
      default:
        throw Exception("P3TriangleFESpace supplies triangles and boundary segments.");
      }
    }
    // Entity queries are used, in particular, to mark Dirichlet DOFs.
    void GetVertexDofNrs(int vertex, Array<DofId> &dnums) const override
    {
      dnums.SetSize(1);
      dnums[0] = vertex;
    }
    void GetEdgeDofNrs(int edge, Array<DofId> &dnums) const override
    {
      dnums.SetSize(2);
      dnums[0] = ma->GetNV() + 2 * edge;
      dnums[1] = ma->GetNV() + 2 * edge + 1;
    }
    void GetFaceDofNrs(int, Array<DofId> &dnums) const override { dnums.SetSize0(); }
    void GetInnerDofNrs(int cell, Array<DofId> &dnums) const override
    {
      // BEGIN EXERCISE_P3_NUMBERING
      // TODO: implement this part of the exercise.
      throw Exception("TODO: EXERCISE_P3_NUMBERING");
      // END EXERCISE_P3_NUMBERING
    }
  };
  static RegisterFESpace<P3TriangleFESpace> register_p3_triangle_space("teachingueb4p3triangle");
} // namespace teaching

#endif // TEACHING_MY_P3_TRIANGLE_FESPACE_HPP
