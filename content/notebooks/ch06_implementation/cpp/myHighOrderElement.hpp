#ifndef TEACHING_MY_HIGH_ORDER_ELEMENT_HPP
#define TEACHING_MY_HIGH_ORDER_ELEMENT_HPP

#include <fem.hpp>
#include <vector>
namespace teaching
{
  using namespace ngfem;
  using namespace ngbla;
  using namespace ngstd;

  // One uniform degree. Arbitrary-order segments, Pp triangles and Qp quads.
  inline int HighOrderNDof(ELEMENT_TYPE et, int p)
  {
    if (p < 1)
      throw Exception("The polynomial degree must be a positive integer.");
    switch (et)
    {
    case ET_SEGM:
      return p + 1;
    case ET_TRIG:
      return (p + 1) * (p + 2) / 2;
    case ET_QUAD:
      return (p + 1) * (p + 1);
    case ET_TET:
      if (p == 1)
        return 4;
      throw Exception("The three-dimensional teaching example uses P1 tetrahedra.");
    default:
      throw Exception("Use segments, triangles, quadrilaterals or P1 tetrahedra.");
    }
  }

  // Scaled Legendre polynomials: S_n(d,s)=s^n P_n(d/s).
  // The recurrence remains polynomial when s=0. No division by s is used.
  template <class T>
  std::vector<T> ScaledLegendre(int p, T d, T s)
  {
    std::vector<T> value(p + 1);
    value[0] = 1.0;
    if (p)
      value[1] = d;
    for (int n = 1; n < p; n++)
      value[n + 1] = ((2 * n + 1) * d * value[n] - n * s * s * value[n - 1]) / (n + 1);
    return value;
  }

  // Integrated Legendre polynomials
  template <class T>
  T IntegratedMode(const std::vector<T> &scaled, int n, T s)
  {
    return (scaled[n] - s * s * scaled[n - 2]) / (2 * n - 1);
  }

  // Jacobi polynomials P_n^(alpha,0), with the standard P_0=1 normalization.
  template <class T>
  std::vector<T> JacobiModes(int p, double alpha, T t)
  {
    std::vector<T> value(p + 1);
    value[0] = 1.0;
    if (p)
      value[1] = 0.5 * (alpha + (alpha + 2) * t);
    for (int n = 1; n < p; n++)
    {
      double c = 2 * n + alpha;
      T numerator =
          (c + 1) * (c * (c + 2) * t + alpha * alpha) * value[n] - 2 * (n + alpha) * n * (c + 2) * value[n - 1];
      double denominator = 2 * (n + 1) * (n + alpha + 1) * c;
      value[n + 1] = numerator / denominator;
    }
    return value;
  }

  // ET selects the reference topology. ET_trait<ET>::DIM supplies its dimension.
  // VertexOrientedFE gives edges and triangular faces a common global vertex ordering.
  template <ELEMENT_TYPE ET>
  class HighOrderElement : public ScalarFiniteElement<ET_trait<ET>::DIM>,
                           public VertexOrientedFE<ET>
  {
    static constexpr int D = ET_trait<ET>::DIM;
    int p;

    // One formula for T=double (values) and T=AutoDiff<D> (values and reference derivatives).
    template <class T>
    void CalcShapes(T x, T y, T z, BareSliceVector<T> shape) const
    {
      if constexpr (ET == ET_SEGM)
      {
        T lam[2] = {x, 1 - x};
        shape(0) = lam[0];
        shape(1) = lam[1];
        if (p == 1)
          return;
        // The returned local endpoint indices are ordered using the supplied global vertex numbers.
        auto edge = this->GetVertexOrientedEdge(0);
        auto pol = ScaledLegendre(p, lam[edge[1]] - lam[edge[0]], T(1.0));
        for (int n = 2; n <= p; n++)
          shape(n) = IntegratedMode(pol, n, T(1.0));
        return;
      }
      if constexpr (ET == ET_TRIG)
      {
        T lam[3] = {x, y, 1 - x - y};
        for (int i = 0; i < 3; i++)
          shape(i) = lam[i];
        if (p == 1)
          return;
        // Local order: vertices, edge modes in increasing degree, then cell-interior modes.
        int k = 3;
        for (int e = 0; e < 3; e++)
        {
          auto edge = this->GetVertexOrientedEdge(e);
          int a = edge[0], b = edge[1];
          T s = lam[a] + lam[b], d = lam[b] - lam[a];
          auto pol = ScaledLegendre(p, d, s);
          for (int n = 2; n <= p; n++)
            shape(k++) = IntegratedMode(pol, n, s);
        }
        if (p >= 3)
        {
          // Use the same vertex ordering as triangular faces of high-order tetrahedra.
          // Optional for 2D cell bubbles, this makes future shared 3D face modes consistent.
          auto face = this->GetVertexOrientedFace(0);
          T s = lam[face[0]] + lam[face[1]];
          T bubble = lam[face[0]] * lam[face[1]] * lam[face[2]];
          auto leg = ScaledLegendre(p - 3, lam[face[1]] - lam[face[0]], s);
          // jac[i][j] is P_j^(2*i+1,0). Each i needs its own Jacobi parameter.
          // The cubic bubble leaves i+j <= p-3, hence a maximal j of p-3-i.
          std::vector<std::vector<T>> jac;
          for (int i = 0; i <= p - 3; i++)
            jac.push_back(JacobiModes(p - 3 - i, 2 * i + 1.0, 2 * lam[face[2]] - 1));
          // Here degree=i+j, so j=degree-i. The bubble adds three to the total degree.
          for (int degree = 0; degree <= p - 3; degree++)
            for (int i = 0; i <= degree; i++)
              shape(k++) = bubble * leg[i] * jac[i][degree - i];
        }
        return;
      }
      if constexpr (ET == ET_QUAD)
      {
        shape(0) = (1 - x) * (1 - y);
        shape(1) = x * (1 - y);
        shape(2) = x * y;
        shape(3) = (1 - x) * y;
        if (p == 1)
          return;
        auto vertices = ElementTopology::GetVertices(ET_QUAD);
        int k = 4;
        for (int e = 0; e < 4; e++)
        {
          auto edge = this->GetVertexOrientedEdge(e);
          int a = edge[0], b = edge[1];
          T t = (vertices[b][0] - vertices[a][0]) * (2 * x - 1.0) + (vertices[b][1] - vertices[a][1]) * (2 * y - 1.0);
          T transverse = shape(a) + shape(b);
          auto pol = ScaledLegendre(p, t, T(1.0));
          for (int n = 2; n <= p; n++)
            shape(k++) = transverse * IntegratedMode(pol, n, T(1.0));
        }
        auto px = ScaledLegendre(p, 2 * x - 1, T(1.0)), py = ScaledLegendre(p, 2 * y - 1, T(1.0));
        for (int i = 2; i <= p; i++)
          for (int j = 2; j <= p; j++)
            shape(k++) = IntegratedMode(px, i, T(1.0)) * IntegratedMode(py, j, T(1.0));
        return;
      }
      if constexpr (ET == ET_TET)
      {
        shape(0) = x;
        shape(1) = y;
        shape(2) = z;
        shape(3) = 1 - x - y - z;
      }
    }

  public:
    // ScalarFiniteElement<D> takes the local DOF count and polynomial order.
    HighOrderElement(int ap, const std::vector<int> &numbers = {})
        : ScalarFiniteElement<D>(HighOrderNDof(ET, ap), ap), p(ap)
    {
      int nv = ElementTopology::GetNVertices(ET);
      if (!numbers.empty() && int(numbers.size()) != nv)
        throw Exception("Wrong vertex count.");
      for (int i = 0; i < nv; i++)
        this->SetVertexNumber(i, numbers.empty() ? i : numbers[i]);
    }

    // The FESpace supplies global vertex IDs through this orientation interface.
    using VertexOrientedFE<ET>::SetVertexNumbers;
    HighOrderElement *SetVertexNumbers(FlatArray<int> numbers) override
    {
      if (numbers.Size() != ElementTopology::GetNVertices(ET))
        throw Exception("Wrong vertex count.");
      VertexOrientedFE<ET>::SetVertexNumbers(numbers);
      return this;
    }

    ELEMENT_TYPE ElementType() const override { return ET; }
    // NGSolve requests all local basis values at one reference integration point.
    void CalcShape(const IntegrationPoint &ip, BareSliceVector<> shape) const override
    {
      CalcShapes<double>(ip(0), D > 1 ? ip(1) : 0.0, D > 2 ? ip(2) : 0.0, shape);
    }
    // Return ndof-by-D reference gradients, with one basis function per row.
    void CalcDShape(const IntegrationPoint &ip, BareSliceMatrix<> dshape) const override
    {
      // AutoDiff<D> propagates a value and D derivatives. Coordinate j seeds derivative j.
      AutoDiff<D> coord[3] = {0.0, 0.0, 0.0};
      for (int j = 0; j < D; j++)
        coord[j] = AutoDiff<D>(ip(j), j);
      Vector<AutoDiff<D>> value(this->GetNDof());
      CalcShapes<AutoDiff<D>>(coord[0], coord[1], coord[2], value);
      for (int i = 0; i < this->GetNDof(); i++)
        for (int j = 0; j < D; j++)
          dshape(i, j) = value(i).DValue(j);
    }
  };
} // namespace teaching

#endif // TEACHING_MY_HIGH_ORDER_ELEMENT_HPP
