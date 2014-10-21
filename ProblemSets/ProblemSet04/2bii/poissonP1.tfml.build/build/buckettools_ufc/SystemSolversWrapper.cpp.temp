
#include "SystemSolversWrapper.h"
#include "BoostTypes.h"
#include <dolfin.h>
#include "PoissonSNES.h"
#include "Errorsolver.h"

namespace buckettools
{
  // A function to return a functionspace from a system given a mesh (defaults to first solver in system as they should all be the same).
  FunctionSpace_ptr ufc_fetch_functionspace(const std::string &systemname,                                             Mesh_ptr mesh, PythonPeriodicMap_ptr periodicmap, MeshFunction_size_t_ptr facetdomains, 
                                            const std::vector<std::size_t> &masterids, const std::vector<std::size_t> &slaveids)
  {
    FunctionSpace_ptr functionspace;
    if (systemname ==  "Poisson")
    {
      // All solvers within a system should return the same functionspace so just take the first one
      if (periodicmap && facetdomains)
      {
        functionspace.reset( new PoissonSNES::FunctionSpace(mesh, periodicmap, facetdomains, masterids, slaveids) );
      }
      else
      {
        functionspace.reset( new PoissonSNES::FunctionSpace(mesh) );
      }
    }
    else if (systemname ==  "Error")
    {
      // All solvers within a system should return the same functionspace so just take the first one
      if (periodicmap && facetdomains)
      {
        functionspace.reset( new Errorsolver::FunctionSpace(mesh, periodicmap, facetdomains, masterids, slaveids) );
      }
      else
      {
        functionspace.reset( new Errorsolver::FunctionSpace(mesh) );
      }
    }
    else
    {
      dolfin::error("Unknown systemname in ufc_fetch_functionspace");
    }
    return functionspace;
  }

  // A function to return a functionspace from a system given a mesh and a solvername.
  FunctionSpace_ptr ufc_fetch_functionspace(const std::string &systemname, const std::string &solvername, 
                                            Mesh_ptr mesh, PythonPeriodicMap_ptr periodicmap, MeshFunction_size_t_ptr facetdomains, 
                                            const std::vector<std::size_t> &masterids, const std::vector<std::size_t> &slaveids)
  {
    FunctionSpace_ptr functionspace;
    if (systemname ==  "Poisson")
    {
      // All solvers within a system should return the same functionspace
      if (solvername ==  "SNES")
      {
        if (periodicmap && facetdomains)
        {
          functionspace.reset(new PoissonSNES::FunctionSpace(mesh, periodicmap, facetdomains, masterids, slaveids) );
        }
        else
        {
          functionspace.reset(new PoissonSNES::FunctionSpace(mesh));
        }
      }
      else
      {
        dolfin::error("Unknown solvername in ufc_fetch_functionspace");
      }
    }
    else if (systemname ==  "Error")
    {
      // All solvers within a system should return the same functionspace
      if (solvername ==  "solver")
      {
        if (periodicmap && facetdomains)
        {
          functionspace.reset(new Errorsolver::FunctionSpace(mesh, periodicmap, facetdomains, masterids, slaveids) );
        }
        else
        {
          functionspace.reset(new Errorsolver::FunctionSpace(mesh));
        }
      }
      else
      {
        dolfin::error("Unknown solvername in ufc_fetch_functionspace");
      }
    }
    else
    {
      dolfin::error("Unknown systemname in ufc_fetch_functionspace");
    }
    return functionspace;
  }

  // A function to return a functionspace (for a coefficient) from a system given a mesh, a solvername and a uflsymbol.
  FunctionSpace_ptr ufc_fetch_coefficientspace_from_solver(const std::string &systemname, const std::string &solvername, const std::string &uflsymbol, 
                                                           Mesh_ptr mesh, PythonPeriodicMap_ptr periodicmap, MeshFunction_size_t_ptr facetdomains, 
                                                           const std::vector<std::size_t> &masterids, const std::vector<std::size_t> &slaveids)
  {
    FunctionSpace_ptr coefficientspace;
    if (systemname ==  "Poisson")
    {
      if (solvername ==  "SNES")
      {
        if (uflsymbol ==  "f")
        {
          coefficientspace.reset(new PoissonSNES::CoefficientSpace_f(mesh));
        }
        else if (uflsymbol ==  "g")
        {
          coefficientspace.reset(new PoissonSNES::CoefficientSpace_g(mesh));
        }
        else
        {
          dolfin::error("Unknown uflsymbol in ufc_fetch_coefficientspace_from_solver");
        }
      }
      else
      {
        dolfin::error("Unknown solvername in ufc_fetch_coefficientspace_from_solver");
      }
    }
    else if (systemname ==  "Error")
    {
      if (solvername ==  "solver")
      {
        dolfin::error("Unknown uflsymbol in ufc_fetch_coefficientspace_from_solver");
      }
      else
      {
        dolfin::error("Unknown solvername in ufc_fetch_coefficientspace_from_solver");
      }
    }
    else
    {
      dolfin::error("Unknown systemname in ufc_fetch_coefficientspace_from_solver");
    }
    return coefficientspace;
  }

  // A function to return a form for a solver from a system given a functionspace, a solvername, a solvertype and a formname.
  Form_ptr ufc_fetch_form(const std::string &systemname, const std::string &solvername, const std::string &solvertype, const std::string &formname, const FunctionSpace_ptr functionspace)
  {
    Form_ptr form;
    if (systemname ==  "Poisson")
    {
      if (solvername ==  "SNES")
      {
        if (solvertype == "SNES")
        {
          if (formname == "Residual")
          {
            form.reset(new PoissonSNES::Form_F(functionspace));
          }
          else if (formname == "Jacobian")
          {
            form.reset(new PoissonSNES::Form_J(functionspace, functionspace));
          }
          else
          {
            dolfin::error("Unknown formname in ufc_fetch_form");
          }
        }
        else
        {
          dolfin::error("Unknown solvertype in ufc_fetch_form");
        }
      }
      else
      {
        dolfin::error("Unknown systemname in ufc_fetch_form");
      }
    }
    else if (systemname ==  "Error")
    {
      if (solvername ==  "solver")
      {
        if (solvertype == "SNES")
        {
          if (formname == "Residual")
          {
            form.reset(new Errorsolver::Form_F(functionspace));
          }
          else if (formname == "Jacobian")
          {
            form.reset(new Errorsolver::Form_J(functionspace, functionspace));
          }
          else
          {
            dolfin::error("Unknown formname in ufc_fetch_form");
          }
        }
        else
        {
          dolfin::error("Unknown solvertype in ufc_fetch_form");
        }
      }
      else
      {
        dolfin::error("Unknown systemname in ufc_fetch_form");
      }
    }
    else
    {
      dolfin::error("Unknown systemname in ufc_fetch_form");
    }
    return form;
  }

}

