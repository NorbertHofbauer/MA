#ifndef VISCOSITYMODELS_HPP
#define VISCOSITYMODELS_HPP

#include "mfem.hpp" // include mfem project

// generalized shear rate
// book 
// Wilczyński
// Rheology in Polymer Processing, page 48

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

// Carreau Constitutive Equation
// book 
// Hopmann 
// Extrusion Dies
// for Plastics and Rubber, page 14
class CarreauModelCoefficient : public mfem::Coefficient
{
protected:
   mfem::GridFunction *u; // displacement
      
   mfem::DenseMatrix grad; // auxiliary matrix, used in Eval

public:
   const double a;
   const double b;
   const double c;
   const double density;

   CarreauModelCoefficient(const double A, const double B, const double C, const double D)
      : u(NULL), a(A), b(B), c(C), density(D) { }

   void SetVelocity(mfem::GridFunction &u_) { u = &u_; }
   
   virtual double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip);
};

#endif // VISCOSITYMODELS_HPP