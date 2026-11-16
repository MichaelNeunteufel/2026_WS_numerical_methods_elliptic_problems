#ifndef TEACHING_MY_P3_TRIANGLE_ELEMENT_HPP
#define TEACHING_MY_P3_TRIANGLE_ELEMENT_HPP
#include <fem.hpp>

namespace teaching
{
  using namespace ngfem;
  using namespace ngbla;
  using namespace ngstd;

  // Fixed reference basis: vertices, two nodes per topology edge, barycenter.
  // No vertex numbers enter this element; the FESpace handles edge orientation.
  class P3TriangleElement : public ScalarFiniteElement<2>
  {
  public:
    // <2> is the reference dimension; (10,3) means ten DOFs, polynomial degree three.
    P3TriangleElement() : ScalarFiniteElement<2>(10, 3) {}
    ELEMENT_TYPE ElementType() const override { return ET_TRIG; }

    void CalcShape(const IntegrationPoint &ip, BareSliceVector<> shape) const override
    {
      double lam[3] = {ip(0), ip(1), 1 - ip(0) - ip(1)};
      auto edges = ElementTopology::GetEdges(ET_TRIG);
      // BEGIN EXERCISE_P3_SHAPES
      // TODO: implement this part of the exercise.
      throw Exception("TODO: EXERCISE_P3_SHAPES");
      // END EXERCISE_P3_SHAPES
    }

    // Rows are basis functions; columns are reference derivatives d/dxi,d/deta.
    void CalcDShape(const IntegrationPoint &ip, BareSliceMatrix<> dshape) const override
    {
      double lam[3] = {ip(0), ip(1), 1 - ip(0) - ip(1)};
      double dlam[3][2] = {{1, 0}, {0, 1}, {-1, -1}};
      auto edges = ElementTopology::GetEdges(ET_TRIG);
      // BEGIN EXERCISE_P3_DSHAPES
      // TODO: implement this part of the exercise.
      throw Exception("TODO: EXERCISE_P3_DSHAPES");
      // END EXERCISE_P3_DSHAPES
    }
  };
} // namespace teaching

#endif // TEACHING_MY_P3_TRIANGLE_ELEMENT_HPP
