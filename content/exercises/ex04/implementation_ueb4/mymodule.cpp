#include <comp.hpp>
#include <python_comp.hpp>
#include <string>
#include <vector>
#include "myP1TriangleFESpace.hpp"
#include "myP3TriangleFESpace.hpp"
#include "myIntegrator.hpp"

namespace teaching
{
  // The helpers expose the same fixed local ordering as the finite elements.
  std::vector<std::vector<double>> ReferenceNodes(const std::string &kind, int degree)
  {
    if (degree != 1 && degree != 3)
      throw Exception("Use degree 1 or 3.");
    if (kind == "segment")
    {
      if (degree == 1)
        return {{1}, {0}};
      return {{1}, {0}, {2.0 / 3}, {1.0 / 3}};
    }
    if (kind != "triangle")
      throw Exception("Use triangle or segment.");
    std::vector<std::vector<double>> nodes = {{1, 0}, {0, 1}, {0, 0}};
    if (degree == 3)
    {
      auto edges = ElementTopology::GetEdges(ET_TRIG);
      for (int e = 0; e < 3; e++)
      {
        auto a = nodes[edges[e][0]], b = nodes[edges[e][1]];
        nodes.push_back({(2 * a[0] + b[0]) / 3, (2 * a[1] + b[1]) / 3});
        nodes.push_back({(a[0] + 2 * b[0]) / 3, (a[1] + 2 * b[1]) / 3});
      }
      nodes.push_back({1.0 / 3, 1.0 / 3});
    }
    return nodes;
  }
}

// CompilePythonModule loads this named entry point; extern C preserves its symbol.
extern "C" NGCORE_API_EXPORT void teaching_fem(py::object &res)
{
  // Import NGSolve's Python types before exporting subclasses and integrators.
  py::module::import("ngsolve");
  static py::module::module_def def;
  py::module m = py::module::create_extension_module("", "", &def);
  using namespace teaching;
  // These bindings supply keyword-to-Flags conversion and initialize each space.
  ExportFESpace<P1TriangleFESpace>(m, "P1TriangleFESpace");
  ExportFESpace<P3TriangleFESpace>(m, "P3TriangleFESpace");
  // The registered base return type makes the result usable in BilinearForm += ... .
  m.def("TensorLaplace", [](shared_ptr<CoefficientFunction> coef, int q) -> shared_ptr<BilinearFormIntegrator>
        { return make_shared<TensorDiffusionIntegrator>(coef, q); }, py::arg("coefficient"), py::arg("quadrature_order") = -1);
  m.def("ReferenceNodes", &ReferenceNodes, py::arg("element"), py::arg("degree"));
  m.def("ReferenceShapes", [](std::string kind, int degree, std::vector<double> point)
        {
    auto nodes = ReferenceNodes(kind,degree); // validates kind and degree
    int dim = kind == "triangle" ? 2 : 1;
    if (point.size() != size_t(dim)) throw Exception("Reference point has the wrong dimension.");
    Vector<> values(nodes.size());
    Matrix<> gradients(nodes.size(),dim);
    IntegrationPoint ip(point[0], dim == 2 ? point[1] : 0, 0, 0);
    if (kind == "triangle" && degree == 1)
    { P1TriangleElement fe; fe.CalcShape(ip,values); fe.CalcDShape(ip,gradients); }
    if (kind == "triangle" && degree == 3)
    { P3TriangleElement fe; fe.CalcShape(ip,values); fe.CalcDShape(ip,gradients); }
    if (kind == "segment" && degree == 1)
    { P1SegmentElement fe; fe.CalcShape(ip,values); fe.CalcDShape(ip,gradients); }
    if (kind == "segment" && degree == 3)
    { P3SegmentElement fe; fe.CalcShape(ip,values); fe.CalcDShape(ip,gradients); }
    std::vector<double> shape(nodes.size());
    std::vector<std::vector<double>> derivative(nodes.size(),std::vector<double>(dim));
    for (size_t i=0; i<nodes.size(); i++)
    {
      shape[i] = values(i);
      for (int j=0; j<dim; j++) derivative[i][j] = gradients(i,j);
    }
    // Return owned Python data rather than views into the temporary NGSolve arrays.
    return py::make_tuple(shape,derivative); }, py::arg("element"), py::arg("degree"), py::arg("point"));
  res = m;
}
