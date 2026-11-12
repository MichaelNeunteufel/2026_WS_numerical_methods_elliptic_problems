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

  // Assemble the full global matrix; Dirichlet DOFs are handled later via the free-DOF mask.
  // FESpace provides basis functions and DOF maps, while MeshAccess provides geometry.
  inline shared_ptr<BaseSparseMatrix> AssembleMatrix(shared_ptr<FESpace> fes,
                                                     const std::vector<shared_ptr<BilinearFormIntegrator>> &integrators)
  {
    CheckAssemblySpace(fes);
    auto mesh = fes->GetMeshAccess();
    Array<DofId> dnums;
    // Build the sparsity graph: all pairs of DOFs belonging to one cell may couple.
    Array<int> counts(mesh->GetNE(VOL));
    for (auto ei : mesh->Elements(VOL))
    {
      fes->GetDofNrs(ei, dnums);
      counts[ei.Nr()] = dnums.Size();
    }
    // A Table stores variable-length DOF lists, one row per volume element.
    Table<int> element_dofs(counts);
    for (auto ei : mesh->Elements(VOL))
    {
      fes->GetDofNrs(ei, dnums);
      element_dofs[ei.Nr()] = dnums;
    }
    // Allocate the global matrix graph first. False requests full rather than symmetric storage.
    auto matrix = make_shared<SparseMatrix<double>>(fes->GetNDof(), fes->GetNDof(), element_dofs, element_dofs, false);
    matrix->SetZero();
    // LocalHeap owns temporary element objects and matrix/vector storage.
    LocalHeap heap(1000000, "teaching assembly");
    for (auto vb : {VOL, BND})
      for (auto ei : mesh->Elements(vb))
      {
        // Release this iteration's temporary heap storage when leaving the element scope.
        HeapReset reset(heap);
        auto element = mesh->GetElement(ei);
        // GetFE returns the reference basis. GetTrafo returns its map to the physical element.
        const FiniteElement &fe = fes->GetFE(ei, heap);
        const ElementTransformation &trafo = mesh->GetTrafo(ei, heap);
        // dnums maps each local basis index to a global matrix/vector index.
        fes->GetDofNrs(ei, dnums);
        // A square local matrix with one row/column per element DOF, allocated from the heap.
        FlatMatrix<> local(fe.GetNDof(), heap);
        // Respect the integration category (VOL/BND) and any region/element restrictions.
        for (auto integrator : integrators)
          if (integrator->VB() == vb && integrator->DefinedOn(element.GetIndex()) && integrator->DefinedOnElement(ei.Nr()))
          {
            integrator->CalcElementMatrix(fe, trafo, local, heap);
            // Scatter-add: contributions from neighboring elements accumulate at shared DOFs.
            matrix->AddElementMatrix(dnums, dnums, local);
          }
      }
    return matrix;
  }

  // Linear-form counterpart: assemble one global right-hand-side entry per DOF.
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
            // Add each local load entry to its global DOF. Shared entries accumulate.
            vector->AddIndirect(dnums, local);
          }
      }
    return vector;
  }
} // namespace teaching

#endif // TEACHING_MY_ASSEMBLING_HPP
