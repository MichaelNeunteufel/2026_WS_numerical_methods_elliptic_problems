#ifndef TEACHING_MY_P1_TRIANGLE_ELEMENT_HPP
#define TEACHING_MY_P1_TRIANGLE_ELEMENT_HPP

#include <fem.hpp>

namespace teaching
{
  using namespace ngfem;
  using namespace ngbla;
  using namespace ngstd;

  // Reference vertices: (1,0), (0,1), (0,0)
  // The template argument <2> is the reference dimension: a triangle is 2D.
  class P1TriangleElement : public ScalarFiniteElement<2>
  {
  public:
    // Base constructor (ndof, order): three local basis functions, degree one.
    P1TriangleElement() : ScalarFiniteElement<2>(3, 1) {}
    // ET_TRIG identifies the reference topology used by NGSolve.
    ELEMENT_TYPE ElementType() const override { return ET_TRIG; }

    // ip contains reference coordinates. Fill the caller's length-three view.
    // Entry alpha is the value of the reference basis function alpha at ip.
    void CalcShape(const IntegrationPoint &ip, BareSliceVector<> shape) const override
    {
      shape(0) = ip(0);
      shape(1) = ip(1);
      shape(2) = 1 - ip(0) - ip(1);
    }

    // Fill a 3-by-2 matrix: basis functions in rows, reference derivatives in columns.
    void CalcDShape(const IntegrationPoint &, BareSliceMatrix<> dshape) const override
    {
      // These affine functions have constant gradients, so ip is unused.
      dshape(0, 0) = 1.0;
      dshape(0, 1) = 0.0;
      dshape(1, 0) = 0.0;
      dshape(1, 1) = 1.0;
      dshape(2, 0) = -1.0;
      dshape(2, 1) = -1.0;
    }
  };

} // namespace teaching

#endif // TEACHING_MY_P1_TRIANGLE_ELEMENT_HPP
