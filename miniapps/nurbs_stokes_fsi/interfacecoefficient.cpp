#include "interfacecoefficient.hpp"
#include <fstream>  // fstream for input and output in textfiles
#include <iostream> // for output to the console

#include <algorithm> // to get access to math functions
#include <cmath>     // to get access to math functions


double InterfaceDirichletCoefficient::Eval(mfem::ElementTransformation &T,
                               const mfem::IntegrationPoint &ip)
{
   MFEM_ASSERT(gf_source != NULL, "gridfunction source is not set");
   MFEM_ASSERT(gf_target != NULL, "gridfunction target is not set");

   //double L = lambda.Eval(T, ip);
   //double M = mu.Eval(T, ip);

   mfem::GridFunctionCoefficient gfc_source(gf_source); 
   mfem::GridFunctionCoefficient gfc_target(gf_target); 
   int sdim = gf_target->FESpace()->GetMesh()->SpaceDimension();
   mfem::Vector phys_point;
   T.Transform(ip, phys_point);

   int ret;
   mfem::IntegrationPoint ip_source;
   int elem_idx;
   mfem::ElementTransformation* T_source;
   for (int i=0; i<gf_source->FESpace()->GetMesh()->GetNE(); ++i)
   {
      T_source = gf_source->FESpace()->GetMesh()->GetElementTransformation(i);
      mfem::InverseElementTransformation invtran(T_source);
      ret = invtran.Transform(phys_point, ip_source);
      if (ret == 0)
      {
         elem_idx = i;
         //std::cout << " source " << gf_source->GetValue(elem_idx, ip_source,1) << " element " << i; 
         //std::cout << " ElementNo " << T_source->ElementNo << " ElementType " << T_source->ElementType << " \n";
         //std::cout << " ElementNo " << T.ElementNo << " ElementType " << T.ElementType << " \n";
         
         break;
      }
   }

   if(ret != 0){
      std::cout << "interface ret !=0 \n";
      return 0;
   }

   double target = 0;
   double source = 0;
   double dirichlet = 0;
   target = gfc_target.Eval(T, ip);
   source = gfc_source.Eval(*T_source, ip_source);
   //source = gf_source->GetValue(elem_idx, ip_source,1);
   /*
   std::cout << " source " << source << " target " << target;
   for (size_t i = 0; i < sdim; i++)
   {
      std::cout << " phys_point[ "<< i << "] " << phys_point[i];
   }
   //std::cout << "\n";
   */
   if (sdim==2)
   {
      //std::cout << " ip.x        " << ip.x <<        " ip.y        " << ip.y << "\n";
      //std::cout << " ip_source.x " << ip_source.x << " ip_source.y " << ip_source.y << "\n";
   }
   
   dirichlet = beta_t*target + (1-beta_t)*source;
   //std::cout << " dirichlet " << dirichlet << " \n";

   return dirichlet;
}

double InterfaceFluxCoefficient::Eval(mfem::ElementTransformation &T,
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
   for (int i=0; i<gf_source->FESpace()->GetMesh()->GetNE(); ++i)
   {
      T_source = gf_source->FESpace()->GetMesh()->GetElementTransformation(i);
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
   double flux = 0;
   double flux_source = 0;
   double flux_target = 0;
   mfem::Vector grad_source(sdim);
   mfem::Vector grad_target(sdim);
   mfem::Vector nhat(sdim);
   T_source->Reset();
   T_source->SetIntPoint(&ip_source);
   T.Reset();
   T.SetIntPoint(&ip);
   gf_source->GetGradient(*T_source, grad_source);
   gf_target->GetGradient(T, grad_target);

   //std::cout << " height " << T_source->Jacobian().Height() << " width " << T_source->Jacobian().Width() << " \n " ;
   //std::cout << " ElementNo " << T_source->ElementNo << " ElementType " << T_source->ElementType << " \n";
   
   //mfem::CalcOrtho(T_source->Jacobian(), nhat);
   mfem::CalcOrtho(T.Jacobian(), nhat);
   nhat *= 1.0 / nhat.Norml1();

   target = gfc_target.Eval(T, ip);
   source = gfc_source.Eval(*T_source, ip_source);
   /*
   std::cout << " source " << source << " target " << target ;
   //std::cout << " source " << source << " target " << target << "k_source " << k_source << "k_target " << k_target;
   for (size_t i = 0; i < sdim; i++)
   {
      std::cout << " phys_point[ "<< i << "] " << phys_point[i];
   }
   //std::cout << "\n";
   */
 /*
   for (size_t i = 0; i < grad_source.Size(); i++)
   {
      std::cout << " grad_source[ "<< i << "] " << grad_source[i];
   }
   std::cout << "\n";
   
   
   for (size_t i = 0; i < nhat.Size(); i++)
   {
      std::cout << " nhat[ "<< i << "] " << nhat[i];
   }
   std::cout << "\n";
   */
   for (size_t i = 0; i < sdim; i++)
   {
      flux_source += grad_source[i] * nhat[i];
      flux_target += grad_target[i] * nhat[i];
   }
     
   flux_source *= k_source;
   flux_target *= k_target;
   flux = beta_q*flux_target + (1-beta_q)*flux_source;
   //std::cout << " delta_flux " << (flux-flux_source);
   dflux->push_back(flux-flux_source);
   //std::cout << " flux " << flux << "\n";
   /*
   if (sdim==2)
   {
      std::cout << " ip.x        " << ip.x <<        " ip.y        " << ip.y << "\n";
      std::cout << " ip_source.x " << ip_source.x << " ip_source.y " << ip_source.y << "\n";
   }*/

   return flux;
}