#ifndef TEACHING_MY_INTEGRATOR_HPP
#define TEACHING_MY_INTEGRATOR_HPP
#include <comp.hpp>

namespace teaching
{
  using namespace ngcomp;

  // Volume integrator for integral_T (K grad phi_j) . grad phi_i.
  // The input is a real symmetric 2-by-2 matrix coefficient function.
  class TensorDiffusionIntegrator : public BilinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> coefficient;
    int quadrature_order;

  public:
    TensorDiffusionIntegrator(shared_ptr<CoefficientFunction> coef, int q)
        : coefficient(coef), quadrature_order(q)
    {
      if (!coef || coef->IsComplex() || coef->Dimension() != 4)
        throw Exception("TensorLaplace requires a real 2x2 coefficient.");
    }
    string Name() const override { return "TeachingTensorLaplace"; }
    int DimElement() const override { return 2; }
    int DimSpace() const override { return 2; }
    VorB VB() const override { return VOL; }
    xbool IsSymmetric() const override { return true; }
    void CalcElementMatrix(const FiniteElement &base,
                           const ElementTransformation &trafo,
                           FlatMatrix<double> matrix, LocalHeap &lh) const override
    {
      auto &fe = dynamic_cast<const BaseScalarFiniteElement &>(base);
      int n = fe.GetNDof();
      matrix = 0;
      FlatMatrix<> reference_gradients(n, 2, lh), gradients(n, 2, lh);
      FlatVector<> values(4, lh);
      FlatMatrix<> tensor(2, 2, lh);
      // Default exact for affine elements and polynomial K of degree at most four.
      int q = quadrature_order >= 0 ? quadrature_order : 2 * fe.Order() + 2;
      IntegrationRule rule(fe.ElementType(), q);
      for (auto &ip : rule)
      {
        MappedIntegrationPoint<2, 2> mip(ip, trafo);
        fe.CalcDShape(ip, reference_gradients);
        coefficient->Evaluate(mip, values);

        // vector to matrix conversion for the tensor coefficient
        for (int i = 0; i < 2; i++)
          for (int j = 0; j < 2; j++)
            tensor(i, j) = values(2 * i + j);

        double weight = ip.Weight() * mip.GetMeasure();

        // Derivative rows transform by J^{-1}, while column gradients use J^{-T}
        // BEGIN EXERCISE_TENSOR_DIFFUSION
        // TODO: implement this part of the exercise.
        throw Exception("TODO: EXERCISE_TENSOR_DIFFUSION");
        // END EXERCISE_TENSOR_DIFFUSION
      }
    }
  };
} // namespace teaching

#endif // TEACHING_MY_INTEGRATOR_HPP
