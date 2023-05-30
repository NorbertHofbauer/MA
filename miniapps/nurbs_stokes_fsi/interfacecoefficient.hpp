#ifndef INTERFACECOEFFICIENT_HPP
#define INTERFACECOEFFICIENT_HPP

#include "mfem.hpp" // include mfem project


class InterfaceDirichletCoefficient : public mfem::Coefficient
{
protected:
   mfem::GridFunction *gf_source; // gridfunction source
   mfem::GridFunction *gf_target; // gridfunction target

public:
   const double beta_t; // relaxation
   InterfaceDirichletCoefficient(const double BETA_T)
      : gf_source(NULL),gf_target(NULL),beta_t(BETA_T){ }

   void SetGridFunctionSource(mfem::GridFunction &gf_source_) { gf_source = &gf_source_; }
   void SetGridFunctionTarget(mfem::GridFunction &gf_target_) { gf_target = &gf_target_; }
   
   virtual double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip);
}
;

class InterfaceFluxCoefficient : public mfem::Coefficient
{
protected:
   mfem::GridFunction *gf_source; // gridfunction source
   mfem::GridFunction *gf_target; // gridfunction target

public:
   const double k_source; // heat transfer coefficient
   const double k_target; // heat transfer coefficient
   const double beta_q; // relaxation

   InterfaceFluxCoefficient(const double K_SOURCE,const double K_TARGET,const double BETA_Q)
      : gf_source(NULL),gf_target(NULL), k_source(K_SOURCE), k_target(K_TARGET), beta_q(BETA_Q){ }

   void SetGridFunctionSource(mfem::GridFunction &gf_source_) { gf_source = &gf_source_; }
   void SetGridFunctionTarget(mfem::GridFunction &gf_target_) { gf_target = &gf_target_; }
   
   virtual double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip);
}
;

#endif // INTERFACECOEFFICIENT_HPP