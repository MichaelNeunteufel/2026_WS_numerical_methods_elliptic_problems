#ifndef TEACHING_MY_P3_SEGMENT_ELEMENT_HPP
#define TEACHING_MY_P3_SEGMENT_ELEMENT_HPP
#include <fem.hpp>

namespace teaching
{
  using namespace ngfem;
  using namespace ngbla;

  // Supplied boundary trace element: endpoints 1,0, then nodes 2/3,1/3.
  // This is the trace of the triangle basis along its local edge direction.
  class P3SegmentElement : public ScalarFiniteElement<1>
  {
    static constexpr double nodes[4] = {1, 0, 2.0 / 3, 1.0 / 3};

  public:
    P3SegmentElement() : ScalarFiniteElement<1>(4, 3) {}
    ELEMENT_TYPE ElementType() const override { return ET_SEGM; }
    void CalcShape(const IntegrationPoint &ip, BareSliceVector<> shape) const override
    {
      for (int i = 0; i < 4; i++)
      {
        shape(i) = 1;
        for (int j = 0; j < 4; j++)
          if (j != i)
            shape(i) *= (ip(0) - nodes[j]) / (nodes[i] - nodes[j]);
      }
    }
    void CalcDShape(const IntegrationPoint &ip, BareSliceMatrix<> dshape) const override
    {
      for (int i = 0; i < 4; i++)
      {
        dshape(i, 0) = 0;
        for (int k = 0; k < 4; k++)
          if (k != i)
          {
            double term = 1 / (nodes[i] - nodes[k]);
            for (int j = 0; j < 4; j++)
              if (j != i && j != k)
                term *= (ip(0) - nodes[j]) / (nodes[i] - nodes[j]);
            dshape(i, 0) += term;
          }
      }
    }
  };
} // namespace teaching

#endif // TEACHING_MY_P3_SEGMENT_ELEMENT_HPP
