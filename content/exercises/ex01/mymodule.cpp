#include <comp.hpp>
#include <python_comp.hpp>

extern "C" NGCORE_API_EXPORT void mymodule(py::object &res)
{
  cout << "called mymodule" << endl;
  // import ngsolve such that python base classes are defined
  py::module::import("ngsolve");

  static py::module::module_def def;
  py::module m = py::module::create_extension_module("", "", &def);
  m.def("answer", []() { return 42; });
  m.def("ndof", [](std::shared_ptr<ngcomp::FESpace> fes) {
    return fes->GetNDof();
  });
  res = m;
}
