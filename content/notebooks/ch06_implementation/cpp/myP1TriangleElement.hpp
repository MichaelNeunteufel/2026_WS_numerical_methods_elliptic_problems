#ifndef TEACHING_MY_P1_TRIANGLE_ELEMENT_HPP
#define TEACHING_MY_P1_TRIANGLE_ELEMENT_HPP

#include <fem.hpp>

namespace teaching
{
  using namespace ngfem;
  using namespace ngbla;
  using namespace ngstd;

  // Reference vertices: (1,0), (0,1), (0,0)
  class P1TriangleElement : public ScalarFiniteElement<2>
  {
  public:
    P1TriangleElement() : ScalarFiniteElement<2>(3, 1) {}
    ELEMENT_TYPE ElementType() const override { return ET_TRIG; }

    void CalcShape(const IntegrationPoint &ip, BareSliceVector<> shape) const override
    {
      shape(0) = ip(0);
      shape(1) = ip(1);
      shape(2) = 1 - ip(0) - ip(1);
    }

    void CalcDShape(const IntegrationPoint &, BareSliceMatrix<> dshape) const override
    {
      // Each row is the reference gradient of one basis function
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
