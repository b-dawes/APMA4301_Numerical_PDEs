
#include "SystemFunctionalsWrapper.h"
#include "BoostTypes.h"
#include <dolfin.h>

namespace buckettools
{
  // A function to return a functionspace (for a coefficient) from a system given a mesh, a functionname and a uflsymbol.
  FunctionSpace_ptr ufc_fetch_coefficientspace_from_functional(const std::string &systemname, const std::string &functionname, const std::string &functionalname, const std::string &uflsymbol, 
                                                               Mesh_ptr mesh, PythonPeriodicMap_ptr periodicmap, MeshFunction_size_t_ptr facetdomains, 
                                                               const std::vector<std::size_t> &masterids, const std::vector<std::size_t> &slaveids)
  {
    FunctionSpace_ptr coefficientspace;
    if (systemname ==  "poisson")
    {
      if (functionname ==  "u")
      {
        dolfin::error("Unknown functionalname in ufc_fetch_coefficientspace_from_functional");
      }
      else if (functionname ==  "RHS")
      {
        dolfin::error("Unknown functionalname in ufc_fetch_coefficientspace_from_functional");
      }
      else
      {
        dolfin::error("Unknown functionname in ufc_fetch_coefficientspace_from_functional");
      }
    }
    else
    {
      dolfin::error("Unknown systemname in ufc_fetch_coefficientspace_from_functional");
    }
    return coefficientspace;
  }

  // A function to return a functionspace (for a coefficient) from a system given a mesh, a coefficientname and a uflsymbol.
  FunctionSpace_ptr ufc_fetch_coefficientspace_from_functional(const std::string &systemname, const std::string &coefficientname, const std::string &uflsymbol, 
                                                               Mesh_ptr mesh, PythonPeriodicMap_ptr periodicmap, MeshFunction_size_t_ptr facetdomains, 
                                                               const std::vector<std::size_t> &masterids, const std::vector<std::size_t> &slaveids)
  {
    FunctionSpace_ptr coefficientspace;
    if (systemname ==  "poisson")
    {
      if (coefficientname ==  "RHS")
      {
        dolfin::error("Unknown functional in ufc_fetch_coefficientspace_from_functional");
      }
      else
      {
        dolfin::error("Unknown coefficientname in ufc_fetch_coefficientspace_from_functional");
      }
    }
    else
    {
      dolfin::error("Unknown systemname in ufc_fetch_coefficientspace_from_functional");
    }
    return coefficientspace;
  }

  // A function to return a functional from a system-function set given a mesh and a functionalname.
  Form_ptr ufc_fetch_functional(const std::string &systemname, const std::string &functionname, const std::string &functionalname, Mesh_ptr mesh)
  {
    Form_ptr functional;
    if (systemname ==  "poisson")
    {
      if (functionname ==  "u")
      {
        dolfin::error("Unknown functionalname in ufc_fetch_functional");
      }
      else if (functionname ==  "RHS")
      {
        dolfin::error("Unknown functionalname in ufc_fetch_functional");
      }
      else
      {
        dolfin::error("Unknown functionname in ufc_fetch_functional");
      }
    }
    else
    {
      dolfin::error("Unknown systemname in ufc_fetch_functional");
    }
    return functional;
  }

  // A function to return a functional for a constant from a system-function set given a mesh.
  Form_ptr ufc_fetch_functional(const std::string &systemname, const std::string &coefficientname, Mesh_ptr mesh)
  {
    Form_ptr functional;
    if (systemname ==  "poisson")
    {
      dolfin::error("Unknown coefficientname in ufc_fetch_functional");
    }
    else
    {
      dolfin::error("Unknown systemname in ufc_fetch_functional");
    }
    return functional;
  }

}

