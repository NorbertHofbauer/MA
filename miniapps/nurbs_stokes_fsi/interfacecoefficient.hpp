#ifndef INTERFACECOEFFICIENT_HPP
#define INTERFACECOEFFICIENT_HPP

#include "mfem.hpp" // include mfem project


class InterfaceCoefficient : public mfem::Coefficient
{
protected:
   mfem::GridFunction *gf_source; // gridfunction source
   mfem::GridFunction *gf_target; // gridfunction target

public:
   const double k; // heat transfer coefficient

   InterfaceCoefficient(const double K)
      : gf_source(NULL),gf_target(NULL), k(K){ }

   void SetGridFunctionSource(mfem::GridFunction &gf_source_) { gf_source = &gf_source_; }
   void SetGridFunctionTarget(mfem::GridFunction &gf_target_) { gf_target = &gf_target_; }
   
   virtual double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip);
};

#endif // INTERFACECOEFFICIENT_HPP