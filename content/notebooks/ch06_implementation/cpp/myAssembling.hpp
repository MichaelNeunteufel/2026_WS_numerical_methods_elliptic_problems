#ifndef TEACHING_MY_ASSEMBLING_HPP
#define TEACHING_MY_ASSEMBLING_HPP

#include <comp.hpp>
#include <vector>

namespace teaching
{
  using namespace ngcomp;

  inline void CheckAssemblySpace(const shared_ptr<FESpace> &fes)
  {
    if (fes->IsComplex() || fes->GetDimension() != 1 || fes->IsParallel())
      throw Exception("The teaching assembler uses a real scalar, non-MPI space.");
  }

  inline shared_ptr<BaseSparseMatrix> AssembleMatrix(shared_ptr<FESpace> fes,
                                                     const std::vector<shared_ptr<BilinearFormIntegrator>> &integrators)
  {
    CheckAssemblySpace(fes);
    auto mesh = fes->GetMeshAccess();
    Array<DofId> dnums;
    Array<int> counts(mesh->GetNE(VOL));
    for (auto ei : mesh->Elements(VOL))
    {
      fes->GetDofNrs(ei, dnums);
      counts[ei.Nr()] = dnums.Size();
    }
    Table<int> element_dofs(counts);
    for (auto ei : mesh->Elements(VOL))
    {
      fes->GetDofNrs(ei, dnums);
      element_dofs[ei.Nr()] = dnums;
    }
    auto matrix = make_shared<SparseMatrix<double>>(fes->GetNDof(), fes->GetNDof(), element_dofs, element_dofs, false);
    matrix->SetZero();
    LocalHeap heap(1000000, "teaching assembly");
    for (auto vb : {VOL, BND})
      for (auto ei : mesh->Elements(vb))
      {
        HeapReset reset(heap);
        auto element = mesh->GetElement(ei);
        const FiniteElement &fe = fes->GetFE(ei, heap);
        const ElementTransformation &trafo = mesh->GetTrafo(ei, heap);
        fes->GetDofNrs(ei, dnums);
        FlatMatrix<> local(fe.GetNDof(), heap);
        for (auto integrator : integrators)
          if (integrator->VB() == vb && integrator->DefinedOn(element.GetIndex()) && integrator->DefinedOnElement(ei.Nr()))
          {
            integrator->CalcElementMatrix(fe, trafo, local, heap);
            matrix->AddElementMatrix(dnums, dnums, local);
          }
      }
    return matrix;
  }

  inline shared_ptr<BaseVector> AssembleVector(shared_ptr<FESpace> fes,
                                               const std::vector<shared_ptr<LinearFormIntegrator>> &integrators)
  {
    CheckAssemblySpace(fes);
    auto vector = make_shared<VVector<double>>(fes->GetNDof());
    *vector = 0.0;
    auto mesh = fes->GetMeshAccess();
    LocalHeap heap(1000000, "teaching load assembly");
    Array<DofId> dnums;
    for (auto vb : {VOL, BND})
      for (auto ei : mesh->Elements(vb))
      {
        HeapReset reset(heap);
        auto element = mesh->GetElement(ei);
        const FiniteElement &fe = fes->GetFE(ei, heap);
        const ElementTransformation &trafo = mesh->GetTrafo(ei, heap);
        fes->GetDofNrs(ei, dnums);
        FlatVector<> local(fe.GetNDof(), heap);
        for (auto integrator : integrators)
          if (integrator->VB() == vb && integrator->DefinedOn(element.GetIndex()) && integrator->DefinedOnElement(ei.Nr()))
          {
            integrator->CalcElementVector(fe, trafo, local, heap);
            vector->AddIndirect(dnums, local);
          }
      }
    return vector;
  }
} // namespace teaching

#endif // TEACHING_MY_ASSEMBLING_HPP
