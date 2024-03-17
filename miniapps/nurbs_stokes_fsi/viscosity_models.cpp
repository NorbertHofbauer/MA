#include "viscosity_models.hpp"
#include <fstream>  // fstream for input and output in textfiles
#include <iostream> // for output to the console

#include <algorithm> // to get access to math functions
#include <cmath>     // to get access to math functions


double ShearRate::Calc(mfem::DenseMatrix &grad)
{
   double shearrate = 0;
   for (size_t n = 0; n < grad.Size(); n++)
   {
      //for (size_t m = 0; m < grad.Size(); m++)5
      for (size_t m = n; m < grad.Size(); m++)
      {
         if (n==m)
         {
            shearrate += 2*grad(n,m)*grad(n,m);
         }else
         {
            //shearrate += 0.5*(grad(n,m)+grad(m,n))*(grad(n,m)+grad(m,n));
            shearrate += (grad(n,m)+grad(m,n))*(grad(n,m)+grad(m,n));
         }
         //std::cout << "shearrate " << shearrate << " n " << n << " m " << m << "\n";
      }
   }
   shearrate = std::sqrt(shearrate);
   
   return shearrate;
}

double ShearRateCoefficient::Eval(mfem::ElementTransformation &T,
                               const mfem::IntegrationPoint &ip)
{
   MFEM_ASSERT(u != NULL, "velocity field is not set");

   //double L = lambda.Eval(T, ip);
   //double M = mu.Eval(T, ip);
   ShearRate shearrate;
   u->GetVectorGradient(T, grad);
   
   /*for (size_t n = 0; n < grad.Size(); n++)
   {
      for (size_t m = 0; m < grad.Size(); m++)
      {
         if (n==m)
         {
            shearrate += 2*grad(n,m)*grad(n,m);
         }else
         {
            shearrate += 0.5*(grad(n,m)+grad(n,m))*(grad(n,m)+grad(n,m));
         }
         //std::cout << "shearrate " << shearrate << " n " << n << " m " << m << "\n";
      }
   }
   shearrate = std::sqrt(shearrate);
   */
   return shearrate.Calc(grad);
}

double CarreauModelCoefficient::Eval(mfem::ElementTransformation &T,
                               const mfem::IntegrationPoint &ip)
{
   MFEM_ASSERT(u != NULL, "velocity field is not set");

   //double L = lambda.Eval(T, ip);
   //double M = mu.Eval(T, ip);
   ShearRate shearrate;
   double dynamic_viscosity = 0;
   double kinematic_viscosity = 0;
   
   u->GetVectorGradient(T, grad);
   
   dynamic_viscosity = k1/std::pow((1+k2*shearrate.Calc(grad)),k3);
   kinematic_viscosity = dynamic_viscosity/density;
   //std::cout << "kinematic_viscosity " << kinematic_viscosity << "\n";

   return kinematic_viscosity;
}

double CarreauWLFModelCoefficient::Eval(mfem::ElementTransformation &T,
                               const mfem::IntegrationPoint &ip)
{
   MFEM_ASSERT(u != NULL, "velocity field is not set");
   MFEM_ASSERT(t != NULL, "temperature field is not set");

   //double L = lambda.Eval(T, ip);
   //double M = mu.Eval(T, ip);
   mfem::GridFunctionCoefficient temperature(t); 
   ShearRate shearrate;
   double temp = temperature.Eval(T, ip);
   double at = 0;
   double logat = 0;
   double dynamic_viscosity = 0;
   double kinematic_viscosity = 0;
   
   u->GetVectorGradient(T, grad);
   
   logat = (8.86*(k4-k5))/(101.6 + k4 -k5) - (8.86*(temp - k5))/(101.6 + temp - k5);
   at = std::pow(10,logat);
   dynamic_viscosity = k1*at/std::pow((1+k2*shearrate.Calc(grad)*at),k3);
   //dynamic_viscosity = k1/std::pow((1+k2*shearrate.Calc(grad)),k3);
   kinematic_viscosity = dynamic_viscosity/density;
   //std::cout << "logat " << logat << " at " << at << " temp " << temp << "\n";
   //std::cout << "kinematic_viscosity " << kinematic_viscosity << "\n";

   return kinematic_viscosity;
}

double PowerLawModelCoefficient::Eval(mfem::ElementTransformation &T,
                               const mfem::IntegrationPoint &ip)
{
   MFEM_ASSERT(u != NULL, "velocity field is not set");
   MFEM_ASSERT(t != NULL, "temperature field is not set");

   //double L = lambda.Eval(T, ip);
   //double M = mu.Eval(T, ip);
   mfem::GridFunctionCoefficient temperature(t); 
   ShearRate shearrate;
   double sr;
   double temp = temperature.Eval(T, ip);
   double m = 0;
   double dynamic_viscosity = 0;
   double kinematic_viscosity = 0;
   
   u->GetVectorGradient(T, grad);
   sr = shearrate.Calc(grad);

   m = m0 * std::exp(-a * (temp-T0));
   if (sr>shearrate0)
   {
      dynamic_viscosity = m * std::pow(sr,n-1);
   }else{
      dynamic_viscosity = m ;
   }
   kinematic_viscosity = dynamic_viscosity/density;
   //std::cout << "m " << m << " sr " << sr << "\n";
   //std::cout << "kinematic_viscosity " << kinematic_viscosity << "\n";

   return kinematic_viscosity;
}
;