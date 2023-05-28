#include "interfacecoefficient.hpp"
#include <fstream>  // fstream for input and output in textfiles
#include <iostream> // for output to the console

#include <algorithm> // to get access to math functions
#include <cmath>     // to get access to math functions


double InterfaceCoefficient::Eval(mfem::ElementTransformation &T,
                               const mfem::IntegrationPoint &ip)
{
   MFEM_ASSERT(gf_source != NULL, "gridfunction source is not set");
   MFEM_ASSERT(gf_target != NULL, "gridfunction target is not set");

   //double L = lambda.Eval(T, ip);
   //double M = mu.Eval(T, ip);

   mfem::GridFunctionCoefficient gfc_source(gf_source); 
   mfem::GridFunctionCoefficient gfc_target(gf_target); 
   int sdim = gf_target->FESpace()->GetMesh()->SpaceDimension();
   //double source = gfc_source.Eval(T, ip);
   mfem::Vector phys_point;
   T.Transform(ip, phys_point);

   int ret;
   mfem::IntegrationPoint ip_source;
   int elem_idx;
   mfem::ElementTransformation* T_source;
   for (int i=0; i<gf_source->FESpace()->GetMesh()->GetNBE(); ++i)
   {
      T_source = gf_source->FESpace()->GetMesh()->GetBdrElementTransformation(i);
      mfem::InverseElementTransformation invtran(T_source);
      ret = invtran.Transform(phys_point, ip_source);
      if (ret == 0)
      {
         elem_idx = i;
         break;
      }
   }

   if(ret != 0){
      std::cout << "interface ret !=0 \n";
      return 0;
   }

   double target = 0;
   double source = 0;
   target = gfc_target.Eval(T, ip);
   source = gfc_source.Eval(*T_source, ip_source);
   
   /*
   std::cout << "k " << k << " source " << source << " target " << target;
   for (size_t i = 0; i < sdim; i++)
   {
      std::cout << " phys_point[ "<< i << "] " << phys_point[i];
   }
   std::cout << "\n";
   */
   // CONJUGATE BOUNDARY NEEDS TO BE MADE

   return k*(target-source);
}
