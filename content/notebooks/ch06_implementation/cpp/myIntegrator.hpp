#ifndef TEACHING_MY_INTEGRATOR_HPP
#define TEACHING_MY_INTEGRATOR_HPP

#include <comp.hpp>

namespace teaching
{
  using namespace ngcomp;

  // NGSolve calls CalcElementMatrix for bilinear forms and CalcElementVector for linear forms.
  // The finite element supplies reference shapes. trafo supplies the physical geometry.
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
    // DimElement is the reference dimension; DimSpace is the physical embedding dimension.
    int DimElement() const override { return D; }
    int DimSpace() const override { return D; }
    // VOL selects cell integration. BND in the boundary classes selects boundary integration.
    VorB VB() const override { return VOL; }
    xbool IsSymmetric() const override { return true; }
    // Fill the supplied local matrix in basis order; lh provides temporary storage.
    void CalcElementMatrix(const FiniteElement &base,
                           const ElementTransformation &trafo,
                           FlatMatrix<double> matrix, LocalHeap &lh) const override
    {
      // Use the scalar-FE interface to access reference values and derivatives.
      auto &fe = dynamic_cast<const BaseScalarFiniteElement &>(base);
      int n = fe.GetNDof();
      matrix = 0.0;
      FlatMatrix<> reference_gradients(n, D, lh), gradients(n, D, lh);
      // The rule order is polynomial exactness, not the number of quadrature points.
      IntegrationRule rule(fe.ElementType(), quadrature_order >= 0 ? quadrature_order : 2 * fe.Order() + 2);
      for (auto &ip : rule)
      {
        // Template arguments are reference and physical dimensions. mip evaluates the map at ip.
        MappedIntegrationPoint<D, D> mip(ip, trafo);
        // ip.Weight() is the reference weight. GetMeasure() is |det J| for a volume cell.
        // Evaluate(mip) samples the coefficient at the mapped physical point.
        double weight = ip.Weight() * mip.GetMeasure() * coefficient->Evaluate(mip);
        fe.CalcDShape(ip, reference_gradients);
        // Gradients are rows: physical derivatives equal reference derivatives times J^{-1}.
        gradients = reference_gradients * mip.GetJacobianInverse();
        matrix += weight * gradients * Trans(gradients);
      }
    }
  };

  // Volume mass term: integrate coefficient * phi_i * phi_j.
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
        // Shape values are evaluated on the reference element. Geometry enters through the weight.
        fe.CalcShape(ip, shape);
        matrix += weight * shape * Trans(shape);
      }
    }
  };

  // Boundary mass term: use the boundary trace element and boundary measure.
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
    // A boundary reference element has dimension D-1 and is embedded in physical dimension D.
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
        // Here GetMeasure() is the arclength/surface factor sqrt(det(J^T J)), not |det J|.
        MappedIntegrationPoint<D - 1, D> mip(ip, trafo);
        double weight = ip.Weight() * mip.GetMeasure() * coefficient->Evaluate(mip);
        fe.CalcShape(ip, shape);
        matrix += weight * shape * Trans(shape);
      }
    }
  };

  // Volume load term: integrate the physical source times each test function.
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
    // Fill one entry per local test function. The assembler adds these to global DOFs.
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

  // Boundary load term: evaluate the prescribed data at mapped boundary points.
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
