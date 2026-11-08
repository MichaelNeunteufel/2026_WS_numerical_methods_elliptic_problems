#include <comp.hpp>
#include <python_comp.hpp>
#include <optional>
#include <string>
#include <vector>
#include "myP1TriangleFESpace.hpp"
#include "myHighOrderFESpace.hpp"
#include "myIntegrator.hpp"
#include "myAssembling.hpp"

namespace teaching
{
  // Convert the element names accepted by the Python evaluation helpers.
  ELEMENT_TYPE ElementTypeFromName(const std::string &name)
  {
    if (name == "segment")
      return ET_SEGM;
    if (name == "triangle")
      return ET_TRIG;
    if (name == "quad")
      return ET_QUAD;
    if (name == "tetrahedron")
      return ET_TET;
    throw Exception("Use segment, triangle, quad or tetrahedron.");
  }

  template <class I>
  shared_ptr<I> Restrict(shared_ptr<I> integrator, std::optional<ngcomp::Region> region)
  {
    if (region)
    {
      if (region->VB() != integrator->VB())
        throw ngstd::Exception("Region has the wrong codimension.");
      integrator->SetDefinedOn(region->Mask());
    }
    return integrator;
  }
}

extern "C" NGCORE_API_EXPORT void teaching_fem(py::object &res)
{
  py::module::import("ngsolve");
  static py::module::module_def def;
  py::module m = py::module::create_extension_module("", "", &def);
  using namespace teaching;

  ExportFESpace<teaching::P1TriangleFESpace>(m, "P1TriangleFESpace");
  ExportFESpace<teaching::HighOrderFESpace>(m, "HighOrderFESpace");

  m.def("Laplace", [](shared_ptr<CoefficientFunction> coef, int dim, int q, std::optional<Region> region) -> shared_ptr<BilinearFormIntegrator>
        {
    if (dim==2) return Restrict<BilinearFormIntegrator>(make_shared<teaching::LaplaceIntegrator<2>>(coef,q),region);
    if (dim==3) return Restrict<BilinearFormIntegrator>(make_shared<teaching::LaplaceIntegrator<3>>(coef,q),region);
    throw Exception("Use dimension 2 or 3."); }, py::arg("coefficient"), py::arg("dim") = 2, py::arg("quadrature_order") = -1, py::arg("definedon") = py::none());

  m.def("Mass", [](shared_ptr<CoefficientFunction> coef, int dim, int q, std::optional<Region> region) -> shared_ptr<BilinearFormIntegrator>
        {
    if (dim==2) return Restrict<BilinearFormIntegrator>(make_shared<teaching::MassIntegrator<2>>(coef,q),region);
    if (dim==3) return Restrict<BilinearFormIntegrator>(make_shared<teaching::MassIntegrator<3>>(coef,q),region);
    throw Exception("Use dimension 2 or 3."); }, py::arg("coefficient"), py::arg("dim") = 2, py::arg("quadrature_order") = -1, py::arg("definedon") = py::none());

  m.def("BoundaryMass", [](shared_ptr<CoefficientFunction> coef, int dim, int q, std::optional<Region> region) -> shared_ptr<BilinearFormIntegrator>
        {
    if (dim==2) return Restrict<BilinearFormIntegrator>(make_shared<teaching::BoundaryMassIntegrator<2>>(coef,q),region);
    if (dim==3) return Restrict<BilinearFormIntegrator>(make_shared<teaching::BoundaryMassIntegrator<3>>(coef,q),region);
    throw Exception("Use dimension 2 or 3."); }, py::arg("coefficient"), py::arg("dim") = 2, py::arg("quadrature_order") = -1, py::arg("definedon") = py::none());

  m.def("Source", [](shared_ptr<CoefficientFunction> coef, int dim, int q, std::optional<Region> region) -> shared_ptr<LinearFormIntegrator>
        {
    if (dim==2) return Restrict<LinearFormIntegrator>(make_shared<teaching::SourceIntegrator<2>>(coef,q),region);
    if (dim==3) return Restrict<LinearFormIntegrator>(make_shared<teaching::SourceIntegrator<3>>(coef,q),region);
    throw Exception("Use dimension 2 or 3."); }, py::arg("coefficient"), py::arg("dim") = 2, py::arg("quadrature_order") = -1, py::arg("definedon") = py::none());

  m.def("BoundarySource", [](shared_ptr<CoefficientFunction> coef, int dim, int q, std::optional<Region> region) -> shared_ptr<LinearFormIntegrator>
        {
    if (dim==2) return Restrict<LinearFormIntegrator>(make_shared<teaching::BoundarySourceIntegrator<2>>(coef,q),region);
    if (dim==3) return Restrict<LinearFormIntegrator>(make_shared<teaching::BoundarySourceIntegrator<3>>(coef,q),region);
    throw Exception("Use dimension 2 or 3."); }, py::arg("coefficient"), py::arg("dim") = 2, py::arg("quadrature_order") = -1, py::arg("definedon") = py::none());

  m.def("AssembleMatrix", &AssembleMatrix, py::arg("fes"), py::arg("integrators"));
  m.def("AssembleVector", &AssembleVector, py::arg("fes"), py::arg("integrators"));

  m.def("HighOrderShapes", [](string name, int p, std::vector<double> point, std::vector<int> vertices)
        {
    auto et=ElementTypeFromName(name);
    int dim=ElementTopology::GetSpaceDim(et),n=HighOrderNDof(et,p);
    if (point.size()!=size_t(dim)) throw Exception("Reference point has the wrong dimension.");
    Vector<> values(n); Matrix<> gradients(n,dim);
    IntegrationPoint ip(point[0],dim>1 ? point[1] : 0.,dim>2 ? point[2] : 0.,0.);
    if (et==ET_SEGM) {HighOrderElement<ET_SEGM> fe(p,vertices);fe.CalcShape(ip,values);fe.CalcDShape(ip,gradients);}
    if (et==ET_TRIG) {HighOrderElement<ET_TRIG> fe(p,vertices);fe.CalcShape(ip,values);fe.CalcDShape(ip,gradients);}
    if (et==ET_QUAD) {HighOrderElement<ET_QUAD> fe(p,vertices);fe.CalcShape(ip,values);fe.CalcDShape(ip,gradients);}
    if (et==ET_TET) {HighOrderElement<ET_TET> fe(p,vertices);fe.CalcShape(ip,values);fe.CalcDShape(ip,gradients);}
    std::vector<double> value(n);
    std::vector<std::vector<double>> derivative(n,std::vector<double>(dim));
    for (int i=0;i<n;i++) { value[i]=values(i); for(int j=0;j<dim;j++) derivative[i][j]=gradients(i,j); }
    return py::make_tuple(value,derivative); }, py::arg("element"), py::arg("order"), py::arg("point"), py::arg("vertex_numbers") = std::vector<int>{});

  res = m;
}
