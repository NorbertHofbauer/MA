#ifndef VISCOSITYMODELS_HPP
#define VISCOSITYMODELS_HPP

#include "mfem.hpp" // include mfem project

class ShearRate
{

public:
   ShearRate(){ }
   
   double Calc(mfem::DenseMatrix &grad);
};

class ShearRateCoefficient : public mfem::Coefficient
{
protected:
   mfem::GridFunction *u; // displacement
   
   mfem::DenseMatrix grad; // auxiliary matrix, used in Eval

public:
   ShearRateCoefficient()
      : u(NULL) { }

   void SetVelocity(mfem::GridFunction &u_) { u = &u_; }
   
   virtual double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip);
};

// A Coefficient for computing the generalized shear rate and the belonging viscosity model
class CarreauModelCoefficient : public mfem::Coefficient
{
protected:
   mfem::GridFunction *u; // displacement
      
   mfem::DenseMatrix grad; // auxiliary matrix, used in Eval

public:
   const double a;
   const double b;
   const double c;

   CarreauModelCoefficient(const double A = 6500, const double B = 0.13, const double C = 0.725)
      : u(NULL), a(A), b(B), c(C) { }

   void SetVelocity(mfem::GridFunction &u_) { u = &u_; }
   
   virtual double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip);
};

#endif // VISCOSITYMODELS_HPP