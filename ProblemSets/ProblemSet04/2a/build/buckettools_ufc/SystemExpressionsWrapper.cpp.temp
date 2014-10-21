
#include "SystemExpressionsWrapper.h"
#include "BoostTypes.h"
#include <dolfin.h>

namespace buckettools
{
  // A function to return an expression for a coefficient from a system given a systemname and a functionname (and its size, shape and private members bucket, system and time.
  Expression_ptr cpp_fetch_expression(const std::string &systemname, const std::string &functionname, const std::string &expressiontype, const std::string &expressionname, const std::size_t &size, const std::vector<std::size_t> &shape, const Bucket *bucket, const SystemBucket *system, const double_ptr time)
  {
    Expression_ptr expression;
    if (systemname ==  "poisson")
    {
      if (functionname == "u")
      {
        if (expressiontype == "initial_condition")
        {
          dolfin::error("Unknown expressionname in cpp_fetch_expression.");
        }
        else if (expressiontype == "boundary_condition")
        {
          dolfin::error("Unknown expressionname in cpp_fetch_expression.");
        }
        else if (expressiontype == "value")
        {
          dolfin::error("Unknown expressionname in cpp_fetch_expression.");
        }
        else
        {
          dolfin::error("Unknown expressiontype in cpp_fetch_expression.");
        }
      }
      else if (functionname ==  "RHS")
      {
        if (expressiontype == "initial_condition")
        {
          dolfin::error("Unknown expressionname in cpp_fetch_expression.");
        }
        else if (expressiontype == "boundary_condition")
        {
          dolfin::error("Unknown expressionname in cpp_fetch_expression.");
        }
        else if (expressiontype == "value")
        {
          dolfin::error("Unknown expressionname in cpp_fetch_expression.");
        }
        else
        {
          dolfin::error("Unknown expressiontype in cpp_fetch_expression.");
        }
      }
      else
      {
        dolfin::error("Unknown functionname in cpp_fetch_expression.");
      }
    }
    else
    {
      dolfin::error("Unknown systemname in cpp_fetch_expression");
    }
    return expression;
  }

  // A function to initialize an expression for a cpp expression given a systemname and a functionname (and a boost shared pointer to the expression to initialize.
  void cpp_init_expression(Expression_ptr expression, const std::string &systemname, const std::string &functionname, const std::string &expressiontype, const std::string &expressionname)
  {
    if (systemname ==  "poisson")
    {
      if (functionname == "u")
      {
        if (expressiontype == "initial_condition")
        {
          dolfin::error("Unknown expressionname in cpp_fetch_expression.");
        }
        else if (expressiontype == "boundary_condition")
        {
          dolfin::error("Unknown expressionname in cpp_fetch_expression.");
        }
        else if (expressiontype == "value")
        {
          dolfin::error("Unknown expressionname in cpp_fetch_expression.");
        }
        else
        {
          dolfin::error("Unknown expressiontype in cpp_init_expression.");
        }
      }
      else if (functionname ==  "RHS")
      {
        if (expressiontype == "initial_condition")
        {
          dolfin::error("Unknown expressionname in cpp_fetch_expression.");
        }
        else if (expressiontype == "boundary_condition")
        {
          dolfin::error("Unknown expressionname in cpp_fetch_expression.");
        }
        else if (expressiontype == "value")
        {
          dolfin::error("Unknown expressionname in cpp_fetch_expression.");
        }
        else
        {
          dolfin::error("Unknown expressiontype in cpp_init_expression.");
        }
      }
      else
      {
        dolfin::error("Unknown functionname in cpp_init_expression.");
      }
    }
    else
    {
      dolfin::error("Unknown systemname in cpp_init_expression");
    }
  }

}

