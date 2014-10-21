
#include "SystemSolversWrapper.h"
#include "BoostTypes.h"
#include <dolfin.h>
#include "poissonNewton.h"

namespace buckettools
{
  // A function to return a functionspace from a system given a mesh (defaults to first solver in system as they should all be the same).
  FunctionSpace_ptr ufc_fetch_functionspace(const std::string &systemname,                                             Mesh_ptr mesh, PythonPeriodicMap_ptr periodicmap, MeshFunction_size_t_ptr facetdomains, 
                                            const std::vector<std::size_t> &masterids, const std::vector<std::size_t> &slaveids)
  {
    FunctionSpace_ptr functionspace;
    if (systemname ==  "poisson")
    {
      // All solvers within a system should return the same functionspace so just take the first one
      if (periodicmap && facetdomains)
      {
        functionspace.reset( new poissonNewton::FunctionSpace(mesh, periodicmap, facetdomains, masterids, slaveids) );
      }
      else
      {
        functionspace.reset( new poissonNewton::FunctionSpace(mesh) );
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
    if (systemname ==  "poisson")
    {
      // All solvers within a system should return the same functionspace
      if (solvername ==  "Newton")
      {
        if (periodicmap && facetdomains)
        {
          functionspace.reset(new poissonNewton::FunctionSpace(mesh, periodicmap, facetdomains, masterids, slaveids) );
        }
        else
        {
          functionspace.reset(new poissonNewton::FunctionSpace(mesh));
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
    if (systemname ==  "poisson")
    {
      if (solvername ==  "Newton")
      {
        if (uflsymbol ==  "f")
        {
          coefficientspace.reset(new poissonNewton::CoefficientSpace_f(mesh));
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
    if (systemname ==  "poisson")
    {
      if (solvername ==  "Newton")
      {
        if (solvertype == "SNES")
        {
          if (formname == "Residual")
          {
            form.reset(new poissonNewton::Form_F(functionspace));
          }
          else if (formname == "Jacobian")
          {
            form.reset(new poissonNewton::Form_J(functionspace, functionspace));
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

