#ifndef TEACHING_MY_P1_SEGMENT_ELEMENT_HPP
#define TEACHING_MY_P1_SEGMENT_ELEMENT_HPP

#include <fem.hpp>

namespace teaching
{
  using namespace ngfem;
  using namespace ngbla;

  // Boundary reference endpoints: 1, 0
  // <1> is the reference dimension, even for a boundary edge in a 2D mesh.
  class P1SegmentElement : public ScalarFiniteElement<1>
  {
  public:
    // Base constructor (ndof, order): two endpoint coefficients, degree one.
    P1SegmentElement() : ScalarFiniteElement<1>(2, 1) {}
    // ET_SEGM selects the one-dimensional reference segment.
    ELEMENT_TYPE ElementType() const override { return ET_SEGM; }

    // Fill the supplied view with basis values at the reference coordinate ip(0).
    void CalcShape(const IntegrationPoint &ip, BareSliceVector<> shape) const override
    {
      shape(0) = ip(0);
      shape(1) = 1 - ip(0);
    }

    // Each row contains one derivative along the reference segment, not physical space.
    void CalcDShape(const IntegrationPoint &, BareSliceMatrix<> dshape) const override
    {
      dshape(0, 0) = 1.0;
      dshape(1, 0) = -1.0;
    }
  };
} // namespace teaching

#endif // TEACHING_MY_P1_SEGMENT_ELEMENT_HPP
