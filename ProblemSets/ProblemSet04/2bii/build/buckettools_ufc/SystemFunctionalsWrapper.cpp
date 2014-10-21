
#include "SystemFunctionalsWrapper.h"
#include "BoostTypes.h"
#include <dolfin.h>
#include "PoissonuL2NormErrorSquared.h"

namespace buckettools
{
  // A function to return a functionspace (for a coefficient) from a system given a mesh, a functionname and a uflsymbol.
  FunctionSpace_ptr ufc_fetch_coefficientspace_from_functional(const std::string &systemname, const std::string &functionname, const std::string &functionalname, const std::string &uflsymbol, 
                                                               Mesh_ptr mesh, PythonPeriodicMap_ptr periodicmap, MeshFunction_size_t_ptr facetdomains, 
                                                               const std::vector<std::size_t> &masterids, const std::vector<std::size_t> &slaveids)
  {
    FunctionSpace_ptr coefficientspace;
    if (systemname ==  "Poisson")
    {
      if (functionname ==  "u")
      {
        if (functionalname ==  "L2NormErrorSquared")
        {
          if (uflsymbol ==  "ue")
          {
            coefficientspace.reset(new PoissonuL2NormErrorSquared::CoefficientSpace_ue(mesh));
          }
          else
          {
            dolfin::error("Unknown uflsymbol in ufc_fetch_coefficientspace_from_functional");
          }
        }
        else
        {
          dolfin::error("Unknown functionalname in ufc_fetch_coefficientspace_from_functional");
        }
      }
      else if (functionname ==  "f")
      {
        dolfin::error("Unknown functionalname in ufc_fetch_coefficientspace_from_functional");
      }
      else if (functionname ==  "g")
      {
        dolfin::error("Unknown functionalname in ufc_fetch_coefficientspace_from_functional");
      }
      else if (functionname ==  "AnalyticSolution")
      {
        dolfin::error("Unknown functionalname in ufc_fetch_coefficientspace_from_functional");
      }
      else
      {
        dolfin::error("Unknown functionname in ufc_fetch_coefficientspace_from_functional");
      }
    }
    else if (systemname ==  "Error")
    {
      if (functionname ==  "error")
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
    if (systemname ==  "Poisson")
    {
      if (coefficientname ==  "f")
      {
        dolfin::error("Unknown functional in ufc_fetch_coefficientspace_from_functional");
      }
      else if (coefficientname ==  "g")
      {
        dolfin::error("Unknown functional in ufc_fetch_coefficientspace_from_functional");
      }
      else if (coefficientname ==  "AnalyticSolution")
      {
        dolfin::error("Unknown functional in ufc_fetch_coefficientspace_from_functional");
      }
      else
      {
        dolfin::error("Unknown coefficientname in ufc_fetch_coefficientspace_from_functional");
      }
    }
    else if (systemname ==  "Error")
    {
      dolfin::error("Unknown coefficentname in ufc_fetch_coefficientspace_from_functional");
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
    if (systemname ==  "Poisson")
    {
      if (functionname ==  "u")
      {
        if (functionalname ==  "L2NormErrorSquared")
        {
          functional.reset(new PoissonuL2NormErrorSquared::Form_e2(mesh));
        }
        else
        {
          dolfin::error("Unknown functionalname in ufc_fetch_functional");
        }
      }
      else if (functionname ==  "f")
      {
        dolfin::error("Unknown functionalname in ufc_fetch_functional");
      }
      else if (functionname ==  "g")
      {
        dolfin::error("Unknown functionalname in ufc_fetch_functional");
      }
      else if (functionname ==  "AnalyticSolution")
      {
        dolfin::error("Unknown functionalname in ufc_fetch_functional");
      }
      else
      {
        dolfin::error("Unknown functionname in ufc_fetch_functional");
      }
    }
    else if (systemname ==  "Error")
    {
      if (functionname ==  "error")
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
    if (systemname ==  "Poisson")
    {
      dolfin::error("Unknown coefficientname in ufc_fetch_functional");
    }
    else if (systemname ==  "Error")
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

