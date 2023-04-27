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
// Oswald
// Polymer Rheology,page 69
class CarreauModelCoefficient : public mfem::Coefficient
{
protected:
   mfem::GridFunction *u; // displacement
      
   mfem::DenseMatrix grad; // auxiliary matrix, used in Eval

public:
   const double k1;
   const double k2;
   const double k3;
   const double density;

   CarreauModelCoefficient(const double K1, const double K2, const double K3, const double D)
      : u(NULL), k1(K1), k2(K2), k3(K3), density(D) { }

   void SetVelocity(mfem::GridFunction &u_) { u = &u_; }
   
   virtual double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip);
};

// Carreau Constitutive Equation
// book 
// Oswald
// Polymer Rheology,page 69
class CarreauWLFModelCoefficient : public mfem::Coefficient
{
protected:
   mfem::GridFunction *u; // displacement
   mfem::GridFunction *t; // temperature
      
   mfem::DenseMatrix grad; // auxiliary matrix, used in Eval

public:
   const double k1;
   const double k2;
   const double k3;
   const double k4;
   const double k5;
   const double density;

   CarreauWLFModelCoefficient(const double K1, const double K2, const double K3, const double K4, const double K5, const double D)
      : u(NULL), k1(K1), k2(K2), k3(K3), k4(K4), k5(K5), density(D) { }

   void SetVelocity(mfem::GridFunction &u_) { u = &u_; }
   void SetTemperature(mfem::GridFunction &t_) { t = &t_; }
   
   virtual double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip);
};

class PowerLawModelCoefficient : public mfem::Coefficient
{
protected:
   mfem::GridFunction *u; // displacement
   mfem::GridFunction *t; // temperature
      
   mfem::DenseMatrix grad; // auxiliary matrix, used in Eval

public:
   const double m0;
   const double n;
   const double a;
   const double T0;
   const double shearrate0;
   const double density;

   PowerLawModelCoefficient(const double M0, const double N, const double A, const double TC,const double SHEARRATE0, const double D)
      : u(NULL), m0(M0), n(N), a(A), T0(TC), shearrate0(SHEARRATE0), density(D) { }

   void SetVelocity(mfem::GridFunction &u_) { u = &u_; }
   void SetTemperature(mfem::GridFunction &t_) { t = &t_; }
   
   virtual double Eval(mfem::ElementTransformation &T, const mfem::IntegrationPoint &ip);
};

#endif // VISCOSITYMODELS_HPP