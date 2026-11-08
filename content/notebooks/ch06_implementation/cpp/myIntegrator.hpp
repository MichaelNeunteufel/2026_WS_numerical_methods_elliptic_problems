#ifndef TEACHING_MY_INTEGRATOR_HPP
#define TEACHING_MY_INTEGRATOR_HPP

#include <comp.hpp>

namespace teaching
{
  using namespace ngcomp;

  // Each class implements one weak-form term. D is the physical dimension.
  template <int D>
  class LaplaceIntegrator : public BilinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> coefficient;
    int quadrature_order;

  public:
    LaplaceIntegrator(shared_ptr<CoefficientFunction> coef, int q)
        : coefficient(coef), quadrature_order(q)
    {
      if (!coef || coef->Dimension() != 1)
        throw Exception("Use a scalar coefficient.");
    }
    string Name() const override { return "TeachingLaplace"; }
    int DimElement() const override { return D; }
    int DimSpace() const override { return D; }
    VorB VB() const override { return VOL; }
    xbool IsSymmetric() const override { return true; }
    void CalcElementMatrix(const FiniteElement &base,
                           const ElementTransformation &trafo,
                           FlatMatrix<double> matrix, LocalHeap &lh) const override
    {
      auto &fe = dynamic_cast<const BaseScalarFiniteElement &>(base);
      int n = fe.GetNDof();
      matrix = 0.0;
      FlatMatrix<> reference_gradients(n, D, lh), gradients(n, D, lh);
      IntegrationRule rule(fe.ElementType(), quadrature_order >= 0 ? quadrature_order : 2 * fe.Order() + 2);
      for (auto &ip : rule)
      {
        MappedIntegrationPoint<D, D> mip(ip, trafo);
        double weight = ip.Weight() * mip.GetMeasure() * coefficient->Evaluate(mip);
        fe.CalcDShape(ip, reference_gradients);
        gradients = reference_gradients * mip.GetJacobianInverse();
        matrix += weight * gradients * Trans(gradients);
      }
    }
  };

  template <int D>
  class MassIntegrator : public BilinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> coefficient;
    int quadrature_order;

  public:
    MassIntegrator(shared_ptr<CoefficientFunction> coef, int q)
        : coefficient(coef), quadrature_order(q)
    {
      if (!coef || coef->Dimension() != 1)
        throw Exception("Use a scalar coefficient.");
    }
    string Name() const override { return "TeachingMass"; }
    int DimElement() const override { return D; }
    int DimSpace() const override { return D; }
    VorB VB() const override { return VOL; }
    xbool IsSymmetric() const override { return true; }
    void CalcElementMatrix(const FiniteElement &base,
                           const ElementTransformation &trafo,
                           FlatMatrix<double> matrix, LocalHeap &lh) const override
    {
      auto &fe = dynamic_cast<const BaseScalarFiniteElement &>(base);
      int n = fe.GetNDof();
      FlatVector<> shape(n, lh);
      matrix = 0.0;
      IntegrationRule rule(fe.ElementType(), quadrature_order >= 0 ? quadrature_order : 2 * fe.Order() + 2);
      for (auto &ip : rule)
      {
        MappedIntegrationPoint<D, D> mip(ip, trafo);
        double weight = ip.Weight() * mip.GetMeasure() * coefficient->Evaluate(mip);
        fe.CalcShape(ip, shape);
        matrix += weight * shape * Trans(shape);
      }
    }
  };

  template <int D>
  class BoundaryMassIntegrator : public BilinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> coefficient;
    int quadrature_order;

  public:
    BoundaryMassIntegrator(shared_ptr<CoefficientFunction> coef, int q)
        : coefficient(coef), quadrature_order(q)
    {
      if (!coef || coef->Dimension() != 1)
        throw Exception("Use a scalar coefficient.");
    }
    string Name() const override { return "TeachingBoundaryMass"; }
    int DimElement() const override { return D - 1; }
    int DimSpace() const override { return D; }
    VorB VB() const override { return BND; }
    xbool IsSymmetric() const override { return true; }
    void CalcElementMatrix(const FiniteElement &base,
                           const ElementTransformation &trafo,
                           FlatMatrix<double> matrix, LocalHeap &lh) const override
    {
      auto &fe = dynamic_cast<const BaseScalarFiniteElement &>(base);
      int n = fe.GetNDof();
      FlatVector<> shape(n, lh);
      matrix = 0.0;
      IntegrationRule rule(fe.ElementType(), quadrature_order >= 0 ? quadrature_order : 2 * fe.Order() + 2);
      for (auto &ip : rule)
      {
        MappedIntegrationPoint<D - 1, D> mip(ip, trafo);
        double weight = ip.Weight() * mip.GetMeasure() * coefficient->Evaluate(mip);
        fe.CalcShape(ip, shape);
        matrix += weight * shape * Trans(shape);
      }
    }
  };

  template <int D>
  class SourceIntegrator : public LinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> coefficient;
    int quadrature_order;

  public:
    SourceIntegrator(shared_ptr<CoefficientFunction> coef, int q)
        : coefficient(coef), quadrature_order(q)
    {
      if (!coef || coef->Dimension() != 1)
        throw Exception("Use a scalar coefficient.");
    }
    string Name() const override { return "TeachingSource"; }
    int DimElement() const override { return D; }
    int DimSpace() const override { return D; }
    VorB VB() const override { return VOL; }
    void CalcElementVector(const FiniteElement &base,
                           const ElementTransformation &trafo,
                           FlatVector<double> vector, LocalHeap &lh) const override
    {
      auto &fe = dynamic_cast<const BaseScalarFiniteElement &>(base);
      FlatVector<> shape(fe.GetNDof(), lh);
      vector = 0.0;
      IntegrationRule rule(fe.ElementType(), quadrature_order >= 0 ? quadrature_order : 2 * fe.Order() + 2);
      for (auto &ip : rule)
      {
        MappedIntegrationPoint<D, D> mip(ip, trafo);
        fe.CalcShape(ip, shape);
        vector += ip.Weight() * mip.GetMeasure() * coefficient->Evaluate(mip) * shape;
      }
    }
  };

  template <int D>
  class BoundarySourceIntegrator : public LinearFormIntegrator
  {
    shared_ptr<CoefficientFunction> coefficient;
    int quadrature_order;

  public:
    BoundarySourceIntegrator(shared_ptr<CoefficientFunction> coef, int q)
        : coefficient(coef), quadrature_order(q)
    {
      if (!coef || coef->Dimension() != 1)
        throw Exception("Use a scalar coefficient.");
    }
    string Name() const override { return "TeachingBoundarySource"; }
    int DimElement() const override { return D - 1; }
    int DimSpace() const override { return D; }
    VorB VB() const override { return BND; }
    void CalcElementVector(const FiniteElement &base,
                           const ElementTransformation &trafo,
                           FlatVector<double> vector, LocalHeap &lh) const override
    {
      auto &fe = dynamic_cast<const BaseScalarFiniteElement &>(base);
      FlatVector<> shape(fe.GetNDof(), lh);
      vector = 0.0;
      IntegrationRule rule(fe.ElementType(), quadrature_order >= 0 ? quadrature_order : 2 * fe.Order() + 2);
      for (auto &ip : rule)
      {
        MappedIntegrationPoint<D - 1, D> mip(ip, trafo);
        fe.CalcShape(ip, shape);
        vector += ip.Weight() * mip.GetMeasure() * coefficient->Evaluate(mip) * shape;
      }
    }
  };
} // namespace teaching

#endif // TEACHING_MY_INTEGRATOR_HPP
