#ifndef TEACHING_MY_P1_SEGMENT_ELEMENT_HPP
#define TEACHING_MY_P1_SEGMENT_ELEMENT_HPP

#include <fem.hpp>

namespace teaching
{
  using namespace ngfem;
  using namespace ngbla;

  // Boundary reference endpoints: 1, 0
  class P1SegmentElement : public ScalarFiniteElement<1>
  {
  public:
    P1SegmentElement() : ScalarFiniteElement<1>(2, 1) {}
    ELEMENT_TYPE ElementType() const override { return ET_SEGM; }

    void CalcShape(const IntegrationPoint &ip, BareSliceVector<> shape) const override
    {
      shape(0) = ip(0);
      shape(1) = 1 - ip(0);
    }

    void CalcDShape(const IntegrationPoint &, BareSliceMatrix<> dshape) const override
    {
      dshape(0, 0) = 1.0;
      dshape(1, 0) = -1.0;
    }
  };
} // namespace teaching

#endif // TEACHING_MY_P1_SEGMENT_ELEMENT_HPP
