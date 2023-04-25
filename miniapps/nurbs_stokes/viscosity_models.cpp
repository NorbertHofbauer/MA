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
   
   u->GetVectorGradient(T, grad);
   
   dynamic_viscosity = a/std::pow((1+b*shearrate.Calc(grad)),c);

   return dynamic_viscosity;
}
;