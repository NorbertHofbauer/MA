//                          MFEM - NURBS Version
//
// equations can be viewed in overleaf
// here we use strongly imposed dirichlet boundary conditions
// https://www.overleaf.com/project/62a3a05e216fe69ec1a7b6ac
//
// Compile with: make nurbs_stokes
//
// Description:  nurbs_stokes with the testcase poiseuille flow
//
// TODO: - for building an MPI application we need to integrate the "Par" everywhere!


#include "mfem.hpp" // include mfem project
#include "nurbs_stokes_fsi_solver.hpp"
#include "viscosity_models.hpp"
#include <fstream>
#include <iostream>

// we don`t want to use namespaces here, otherwise code is harder to unterstand, e.g. mfem::Vector != std::vector
//using namespace std; // bad idea!!!!
//using namespace mfem;

// we build one main in this file, self build classes have to be included later if necessary
int main(int argc, char *argv[])
{
   // Setup
   // we define the standard values for our boundaries, equation constants and the meshfile
   double density = 1;    // value for density
   double temp_diffusion_const_fluid = 0.1; //fluid
   double temp_diffusion_const_solid = 0.5; //fluid
   const char *jobname = "jobname";  //jobname will be a prefix for results
   const char *mesh_file_fluid = "../../../MA/data/nurbs_fluid_domain.mesh";  //our standard test mesh
   const char *mesh_file_solid = "../../../MA/data/nurbs_solid_domain.mesh";  //our standard test mesh
   int vis_model = 1;         // type of the viscosity model
   // 1 Carreau
   mfem::Vector model_parameters; // array to get the model parameters
   double relaxation = 0.3;
   double beta_q = 0.3;
   double beta_t = 0.3;
   int ref_levels = 0;              // standard number of refinements for the mesh
   int maxIter = 5000;
   int maxIter2 = 10;
   double rtol = 1.e-6; // convergence criteria
   double atol = 1.e-5; // convergence criteria
   bool visualization = 1;
   int orderelevationvelocity = 0;
   int orderelevationpressure = 0;
   int orderelevationtemperature_fluid = 0;
   int orderelevationtemperature_solid = 0;
   //double v_error_norm = 1;
   //double p_error_norm = 1;
   //double tf_error_norm = 1;
   //double ts_error_norm = 1;
   //double cht_tf_error_norm = 1;
   //double cht_ts_error_norm = 1;
   //double cht_dflux_error_norm = 1;
   double v_error_norm_abs = 1;
   double p_error_norm_abs = 1;
   double tf_error_norm_abs = 1;
   double ts_error_norm_abs = 1;
   double cht_tf_error_norm_abs = 1;
   double cht_ts_error_norm_abs = 1;
   double cht_dflux_error_norm_abs = 1;
   double v_error_norm_rel = 1;
   double p_error_norm_rel = 1;
   double tf_error_norm_rel = 1;
   double ts_error_norm_rel = 1;
   double cht_tf_error_norm_rel = 1;
   double cht_ts_error_norm_rel = 1;
   double cht_dflux_error_norm_rel = 1;
   std::vector<double> flux0;
   std::vector<double> flux;
   double max_error = 1e-3;
   mfem::Array<int> vdbc_bdr_noslip;
   mfem::Array<int> vdbc_bdr;
   mfem::Vector vdbc_bdr_values;
   mfem::Array<int> pdbc_bdr;
   mfem::Vector pdbc_bdr_values;
   mfem::Array<int> tfdbc_bdr;
   mfem::Vector tfdbc_bdr_values;
   mfem::Array<int> tsdbc_bdr;
   mfem::Vector tsdbc_bdr_values;
   mfem::Array<int> tfiface_bdr;
   mfem::Array<int> tsiface_bdr;
   double init_temp;
   
   //./nurbs_stokes_fsi -mf ../../../MA/data/nurbs_fluid_domain.mesh -ms ../../../MA/data/nurbs_solid_domain.mesh -r 2 -vm 3 -mp '2e+4 0.28 -0.025 170 10'  -vnos '1 2 9 6 4 7 10 12' -vdbc '3 8' -vdbc_values '5.5 0 5.5 0' -pdbc '5 11' -pdbc_values '10 10' -tfdbc '3 8' -tfdbc_values '220 220' -tsdbc '1' -tsdbc_values '250' -d 1 -tfdc 0.1 -tsdc 0.5 -tfiface '2 6 9' -tsiface '2 3 4' -mi 10000 -mi2 20 -oev 2 -oep 1 -oetf 0 -oets 0 -rel 0.5 -betaq 0.3 -betat 0.3 -me 1e-8 -at 1e-7 -rt 1e-8  -job test

   // casacade
   //./nurbs_stokes_fsi -mf ../../../MA/data/fluid.mesh -ms ../../../MA/data/solid.mesh -r 2 -vm 3 -mp '2e+4 0.28 -0.025 170 10'  -vnos '1 3 4 5 12 15 14 16 18 7 9 11 10 13 24 23 26 28 17 20 22 25' -vdbc '2 6 8' -vdbc_values '5.5 0 5.5 0 5.5 0' -pdbc '27 21 19' -pdbc_values '0 0 0' -tfdbc '2 6 8' -tfdbc_values '220 220 220' -tsdbc '1 7' -tsdbc_values '250 250' -d 1 -tfdc 0.1 -tsdc 0.5 -tfiface '12 13' -tsiface '3 5' -mi 10000 -mi2 2 -oev 2 -oep 1 -oetf 0 -oets 0 -rel 0.5 -betaq 0.3 -betat 0.3 -me 1e-8 -at 1e-7 -rt 1e-8 -job test

   //cascade2
   //./nurbs_stokes_fsi -mf ../../../MA/data/fluid.mesh -ms ../../../MA/data/solid.mesh -r 2 -vm 3 -mp '2e+4 0.28 -0.025 170 10'  -vnos '1 3 4 5 12 15 14 16 18 7 9 11 10 13 24 23 26 28 17 20 22 25' -vdbc '2 6 8' -vdbc_values '5.5 0 5.5 0 5.5 0' -pdbc '27 21 19' -pdbc_values '0 0 0' -tfdbc '2 6 8' -tfdbc_values '220 220 220' -tsdbc '1 7' -tsdbc_values '250 250' -d 1 -tfdc 0.1 -tsdc 0.5 -tfiface '12 13 5 15 10 24' -tsiface '3 5 6 8 4 2' -mi 10000 -mi2 2 -oev 2 -oep 1 -oetf 0 -oets 0 -rel 1 -betaq 0.3 -betat 0.3 -me 1e-8 -at 1e-7 -rt 1e-8 -job test

   // validation
   // flow and visco models
   //./nurbs_stokes_fsi -mf ../../../MA/data/quad_nurbs.mesh -ms ../../../MA/data/quad_nurbs.mesh -r 5 -vm 0 -mp '200'  -vnos '1' -vdbc '3' -vdbc_values '10 0' -pdbc '2' -pdbc_values '999' -tfdbc '1 3' -tfdbc_values '220 230' -tsdbc '1' -tsdbc_values '250' -d 1 -tfdc 0.1 -tsdc 0.5 -tfiface '' -tsiface '' -mi 10000 -mi2 20 -oev 2 -oep 1 -oetf 0 -oets 0 -rel 1 -betaq 0.3 -betat 0.3 -me 1e-8 -at 1e-8 -rt 1e-8  -job test
   //./nurbs_stokes_fsi -mf ../../../MA/data/quad_nurbs.mesh -ms ../../../MA/data/quad_nurbs.mesh -r 3 -vm 1 -mp '305 0.00046 0.48'  -vnos '1' -vdbc '3' -vdbc_values '10 0' -pdbc '2' -pdbc_values '10' -tfdbc '1 3' -tfdbc_values '210 250' -tsdbc '1' -tsdbc_values '250' -d 1 -tfdc 0.1 -tsdc 0.5 -tfiface '' -tsiface '' -mi 10000 -mi2 10 -oev 2 -oep 1 -oetf 0 -oets 0 -rel 1 -betaq 0.3 -betat 0.3 -me 1e-8 -at 1e-7 -rt 1e-8  -job test
   //./nurbs_stokes_fsi -mf ../../../MA/data/quad_nurbs.mesh -ms ../../../MA/data/quad_nurbs.mesh -r 3 -vm 2 -mp '305 0.00046 0.48 320 153'  -vnos '1' -vdbc '3' -vdbc_values '10 0' -pdbc '2' -pdbc_values '10' -tfdbc '1 3' -tfdbc_values '210 250' -tsdbc '1' -tsdbc_values '250' -d 1 -tfdc 0.1 -tsdc 0.5 -tfiface '' -tsiface '' -mi 10000 -mi2 10 -oev 2 -oep 1 -oetf 0 -oets 0 -rel 1 -betaq 0.3 -betat 0.3 -me 1e-8 -at 1e-7 -rt 1e-8  -job test
   //./nurbs_stokes_fsi -mf ../../../MA/data/quad_nurbs.mesh -ms ../../../MA/data/quad_nurbs.mesh -r 3 -vm 3 -mp '2.8e+4 0.28 -0.025 170 10'  -vnos '1' -vdbc '3' -vdbc_values '10 0' -pdbc '2' -pdbc_values '10' -tfdbc '1 3' -tfdbc_values '210 250' -tsdbc '1' -tsdbc_values '250' -d 1 -tfdc 0.1 -tsdc 0.5 -tfiface '' -tsiface '' -mi 10000 -mi2 10 -oev 2 -oep 1 -oetf 0 -oets 0 -rel 1 -betaq 0.3 -betat 0.3 -me 1e-8 -at 1e-7 -rt 1e-8  -job test

   // temperature field and coupling

   //./nurbs_stokes_fsi -mf ../../../MA/data/validate_2_fluid.mesh -ms ../../../MA/data/validate_2_solid.mesh -r 3 -vm 0 -mp '200' -vnos '' -vdbc '' -vdbc_values '' -pdbc '' -pdbc_values '' -tfdbc '' -tfdbc_values '' -tsdbc '4 6' -tsdbc_values '200 220' -d 1 -tfdc 0.1 -tsdc 0.5 -tfiface '2 4' -tsiface '2 8' -mi 1000 -mi2 10 -oev 2 -oep 1 -oetf 0 -oets 0 -rel 0.5 -betaq 0.3 -betat 0.3 -me 1e-8 -at 1e-7 -rt 1e-8  -job test

   // Parse command-line options.
   // input options for our executable
   mfem::OptionsParser args(argc, argv);
   args.AddOption(&jobname, "-job", "--jobname",
                  "Jobname prefix for results",true);
   args.AddOption(&mesh_file_fluid, "-mf", "--mesh",
                  "Fluid mesh file to use.",true);
   args.AddOption(&mesh_file_solid, "-ms", "--mesh",
                  "Solid mesh file to use.",true);
   args.AddOption(&ref_levels, "-r", "--refine",
                  "Number of times to refine the mesh uniformly, -1 for auto.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&vis_model, "-vm", "--viscosityModel", 
                  "Choose viscosity Model e.g. -vm 1 \n "
                  "Available Models: "
                  "0 No Model, Parameters dynamic viscosity"
                  "1 Carreau, Parameters k1,k2,k3"
                  "2 CarreauWLF, Parameters k1,k2,k3,k4,k5"
                  "3 PowerLaw, Parameters m0,n,a,T0,shearrate0",true);
   args.AddOption(&model_parameters, "-mp", "--modelParameters", 
                  " e.g. -mp 1 2 6",true);
   args.AddOption(&orderelevationvelocity, "-oev", "--orderelevationvelocity", 
                  "the order of the elevation for velocity. e.g. -oev 1");
   args.AddOption(&orderelevationpressure, "-oep", "--orderelevationpressure", 
                  "the order of the elevation for pressure. e.g. -oep 1");
   args.AddOption(&orderelevationtemperature_fluid, "-oetf", "--orderelevationtemperature_fluid", 
                  "the order of the elevation for fluid temperature. e.g. -oetf 1");
   args.AddOption(&orderelevationtemperature_solid, "-oets", "--orderelevationtemperature_solid", 
                  "the order of the elevation for solid temperature. e.g. -oets 1");
   args.AddOption(&density, "-d", "--density", 
                  "fluid density. e.g. -d 1000",true);
   args.AddOption(&relaxation, "-rel", "--relaxation", 
                  "Relaxation for outer loop. e.g. -rel 0.5");
   args.AddOption(&temp_diffusion_const_fluid, "-tfdc", "--temp_diffusion_constant_fluid", 
                  "temperature diffusion constant fluid. e.g. -tfdc 1000",true);
   args.AddOption(&temp_diffusion_const_solid, "-tsdc", "--temp_diffusion_constant_solid", 
                  "temperature diffusion constant solid. e.g. -tsdc 5",true);
   args.AddOption(&maxIter, "-mi", "--maxIteration", 
                  "maximum Iterations for the inner solver. e.g. -mi 500");
   args.AddOption(&maxIter2, "-mi2", "--maxIteration2", 
                  "maximum Iterations for the outer loop. e.g. -mi2 10");
   //deprecated
   args.AddOption(&max_error, "-me", "--maxError", 
                  "convergence criteria: maximum error between iterations in the outer loop. e.g. -me 1.e-3");
   // ^^^
   args.AddOption(&rtol, "-rt", "--rtol", 
                  "convergence criteria relative tolerance. e.g. -rt 1.e-6");
   args.AddOption(&atol, "-at", "--atol", 
                  "convergence criteria absolute tolerance. e.g. -at 1.e-5");
   args.AddOption(&vdbc_bdr_noslip, "-vnos", "--velocity_bdr_noslip", 
                  "Velocity Boundaries with no slip condition. e.g. -vnos 1 2 6 5 ...",true);
   args.AddOption(&vdbc_bdr, "-vdbc", "--velocity_dirchlet_bdr", 
                  "Velocity Dirichlet Boundaries. e.g. -vdbc 1 2 6 5 ...",true);   
   args.AddOption(&vdbc_bdr_values, "-vdbc_values", "--velocity_dirchlet_bdr_values", 
                  "Vector Values for the Velocity Dirichlet Boundaries. e.g. for 3D you need 3 Values for 1 boundary -vdbc_values 1.1 2.9 6.7 ...",true);
   args.AddOption(&pdbc_bdr, "-pdbc", "--pressure_dirchlet_bdr", 
                  "Pressure Dirichlet Boundaries. e.g. -pdbc 1 2 6 5 ...");   
   args.AddOption(&pdbc_bdr_values, "-pdbc_values", "--pressure_dirchlet_bdr_values", 
                  "Values for the Pressure Dirichlet Boundaries. e.g. -pdbc_values 1.1 2.9 6.7 5.22 ...");
   args.AddOption(&tfdbc_bdr, "-tfdbc", "--temperature_dirchlet_bdr_fluid", 
                  "Temperature Dirichlet Boundaries. e.g. -tfdbc 1 2 6 5 ...",true);   
   args.AddOption(&tfdbc_bdr_values, "-tfdbc_values", "--temperature_dirchlet_bdr_values_fluid", 
                  "Values for the Fluid Temperature Dirichlet Boundaries. e.g. -tfdbc_values 1.1 2.9 6.7 5.22 ...",true);
   args.AddOption(&tsdbc_bdr, "-tsdbc", "--temperature_dirchlet_bdr_solid", 
                  "Temperature Dirichlet Boundaries. e.g. -tsdbc 1 2 6 5 ...",true);   
   args.AddOption(&tsdbc_bdr_values, "-tsdbc_values", "--temperature_dirchlet_bdr_values_solid", 
                  "Values for the Solid Temperature Dirichlet Boundaries. e.g. -tsdbc_values 1.1 2.9 6.7 5.22 ...",true);
   args.AddOption(&tfiface_bdr, "-tfiface", "--temperature_iface_fluid", 
                  "Temperature interface fluid. e.g. -tfiface 2 9 6 ...",true);
   args.AddOption(&tsiface_bdr, "-tsiface", "--temperature_iface_solid", 
                  "Temperature interface solid. e.g. -tsiface 2 3 4 ...",true);   
   args.AddOption(&beta_q, "-betaq", "--relaxation_betaq", 
                  "Relaxation for heat transfer flux. e.g. -betaq 0.5");
   args.AddOption(&beta_t, "-betat", "--relaxation_betat", 
                  "Relaxation for heat transfer dirichlet. e.g. -betat 0.5");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(std::cout);
      return 1;
   }
   args.PrintOptions(std::cout);

   NurbsStokesSolver nssolver;
   nssolver.jobname = jobname;
   nssolver.meshfile_fluid = mesh_file_fluid;
   nssolver.meshfile_solid = mesh_file_solid;
   nssolver.density = density;     // density of the fluid
   nssolver.temp_diffusion_const_fluid = temp_diffusion_const_fluid; // value for temperature diffusion constant coefficient
   nssolver.temp_diffusion_const_solid = temp_diffusion_const_solid; // value for temperature diffusion constant coefficient
   nssolver.ref_levels = ref_levels;
   nssolver.set_order_elevation_velocity(orderelevationvelocity);
   nssolver.set_order_elevation_pressure(orderelevationpressure);
   nssolver.set_order_elevation_temperature_fluid(orderelevationtemperature_fluid);
   nssolver.set_order_elevation_temperature_solid(orderelevationtemperature_solid);
   nssolver.maxIter = maxIter;

   nssolver.rtol = rtol; // convergence criteria
   nssolver.atol = atol; // convergence criteria
   //nssolver.max_error = max_error;
   nssolver.beta_q = 0.3;
   
   nssolver.set_dirichletbc_velocity_noslip(vdbc_bdr_noslip);
   nssolver.set_dirichletbc_velocity(vdbc_bdr, vdbc_bdr_values);
   nssolver.set_dirichletbc_pressure(pdbc_bdr, pdbc_bdr_values);
   nssolver.set_dirichletbc_temperature_fluid(tfdbc_bdr, tfdbc_bdr_values);
   nssolver.set_dirichletbc_temperature_solid(tsdbc_bdr, tsdbc_bdr_values);
   nssolver.set_iface_fluid(tfiface_bdr);
   nssolver.set_iface_solid(tsiface_bdr);

   nssolver.init();

   mfem::GridFunction vr(nssolver.vfes),pr(nssolver.pfes),tfr(nssolver.tffes),tsr(nssolver.tsfes); // gridfunctions needed for relaxation 
   mfem::GridFunction v_last(nssolver.vfes),p_last(nssolver.pfes),tf_last(nssolver.tffes),ts_last(nssolver.tsfes); // gridfunctions needed for convergence of outer loop
   mfem::GridFunction v_current(nssolver.vfes),p_current(nssolver.pfes),tf_current(nssolver.tffes),ts_current(nssolver.tsfes); // gridfunctions needed for convergence of outer loop
   mfem::GridFunction cht_tf0_last(nssolver.tffes),cht_tf0_current(nssolver.tffes),cht_ts0_last(nssolver.tsfes),cht_ts0_current(nssolver.tsfes); //gridfunctions needed for convergence of inner loop
   mfem::GridFunction v0(nssolver.vfes),p0(nssolver.pfes),tf0(nssolver.tffes),v(nssolver.vfes),p(nssolver.pfes),tf(nssolver.tffes),ts0(nssolver.tsfes),ts(nssolver.tsfes);
   mfem::GridFunction cht_tf0(nssolver.tffes),cht_ts0(nssolver.tsfes);
   nssolver.visualization = 0;
   nssolver.calc_dirichletbc_fluid(v0,p0,tf0);
   nssolver.visualization = 0;
   nssolver.calc_dirichletbc_solid(ts0);
   nssolver.visualization = 0;

   vr = v0;
   pr = p0;
   tfr = tf0;
   tsr = ts0;

   //for (size_t i = 0; i < 28; i++)
   int iter=0;
   //while ((v_error_norm>max_error)||(p_error_norm>max_error)||(tf_error_norm>max_error)||(ts_error_norm>max_error))
   while (((v_error_norm_abs>atol)&&(v_error_norm_rel>rtol))||((p_error_norm_abs>atol)&&(p_error_norm_rel>rtol))||((tf_error_norm_abs>atol)&&(tf_error_norm_rel>rtol))||((ts_error_norm_abs>atol)&&(ts_error_norm_rel>rtol)))
   //while (false)
   {  
      iter+=1;
      nssolver.iter = iter;
      if ((iter==-1)|(iter==maxIter2))
      {
         nssolver.visualization = 1;
      }

      // viscosity model - must be a mfem::coefficient or a child class from coefficient
      if (vis_model==0)
      {
         MFEM_ASSERT(model_parameters.Size() == 1, "1 Parameter needed dynamic viscosity!");
         mfem::ConstantCoefficient kin_vis(model_parameters[0]/density);
         nssolver.solve_flow(v0,p0,tf0,v,p,tf,kin_vis);
      }else if (vis_model==1)
      {
         MFEM_ASSERT(model_parameters.Size() == 3, "Carreau Model needs 3 Parameters k1,k2,k3!");
         CarreauModelCoefficient kin_vis(model_parameters[0],model_parameters[1],model_parameters[2], density);
         nssolver.solve_flow(v0,p0,tf0,v,p,tf,kin_vis);
      }else if (vis_model==2)
      {
         MFEM_ASSERT(model_parameters.Size() == 5, "CarreauWLF Model needs 5 Parameters k1,k2,k3,k4,k5!");
         CarreauWLFModelCoefficient kin_vis(model_parameters[0],model_parameters[1],model_parameters[2],model_parameters[3],model_parameters[4], density);
         //set reference temperature on first iteration, otherwise model crashes
         if (iter==1)
         {
            for (size_t i = 0; i < tf0.Size(); i++)
            {
               if (tf0[i] == 0)
               {
                  tf0[i] = model_parameters[3];
               }
               //std::cout <<  t0[i] << " t0 \n";
            }
         }
         nssolver.solve_flow(v0,p0,tf0,v,p,tf,kin_vis);
      }else if (vis_model==3)
      {
         MFEM_ASSERT(model_parameters.Size() == 5, "PowerLaw Model needs 5 Parameters m0,n,a,T0,shearrate0!");
         PowerLawModelCoefficient kin_vis(model_parameters[0],model_parameters[1],model_parameters[2],model_parameters[3],model_parameters[4], density);
         if (iter==1)
         {
            for (size_t i = 0; i < tf0.Size(); i++)
            {
               if (tf0[i] == 0)
               {
                  tf0[i] = model_parameters[3];
               }
               //std::cout <<  t0[i] << " t0 \n";
            }
         }
         nssolver.solve_flow(v0,p0,tf0,v,p,tf,kin_vis);
      }else{
         MFEM_ASSERT(false, "Viscosity Model not available!");
      }
      
      v_current = v;
      p_current = p;

      //v0 = (1-relaxation)*vr + relaxation*v;
      //p0 = (1-relaxation)*tr + relaxation*p;
      //v0 = v;
      //p0 = p;
      
      for (size_t i = 0; i < vr.Size(); i++)
      {
         bool skip=false;
         for (size_t ii = 0; ii < nssolver.vel_ess_tdof_list.Size(); ii++)
         {
            if (i==nssolver.vel_ess_tdof_list[ii])
            {
               skip=true;
            }
         }
         if (!skip)
         {
            v0[i] = (1-relaxation)*vr[i] + relaxation*v[i];
         }
      }
      for (size_t i = 0; i < pr.Size(); i++)
      {  
         bool skip=false;
         for (size_t ii = 0; ii < nssolver.pres_ess_tdof_list.Size(); ii++)
         {
            if (i==nssolver.pres_ess_tdof_list[ii])
            {
               skip=true;
            }
         }
         if (!skip)
         {
            p0[i] = (1-relaxation)*pr[i] + relaxation*p[i];
         }
      }

      vr = v0;
      pr = p0;
      

      nssolver.visualization = 0;
      if ((iter==-1)|(iter==maxIter2))
      {
         nssolver.visualization = 1;
      }
      
      //cht_tf_error_norm = 1;
      //cht_ts_error_norm = 1;
      //cht_dflux_error_norm = 1;

      cht_tf_error_norm_abs = 1;
      cht_ts_error_norm_abs = 1;
      cht_dflux_error_norm_abs = 1;
      cht_tf_error_norm_rel = 1;
      cht_ts_error_norm_rel = 1;
      cht_dflux_error_norm_rel = 1;

      flux.clear();

      cht_tf0 = tf0;
      cht_ts0 = ts0;
      cht_tf0_last = tf0;
      cht_ts0_last = ts0;
      
      int iter2 = 0;
      while (((cht_tf_error_norm_abs>atol)&&(cht_tf_error_norm_rel>rtol))||((cht_ts_error_norm_abs>atol)&&(cht_ts_error_norm_rel>rtol))||((cht_dflux_error_norm_abs>atol)&&(cht_dflux_error_norm_rel>rtol)))
      //while ((cht_tf_error_norm>max_error)||(cht_ts_error_norm>max_error))
      {  
         iter2+=1;
                  
         //cht_tf_error_norm = cht_tf0.Norml2();
         //cht_ts_error_norm = cht_ts0.Norml2();

         if (iter2!=1)
         {  
            //cht_dflux_error_norm = 0;
            flux0 = flux;
            flux.clear();
         }
         
         if ((iter==1)&&(iter2==1))
         {
            nssolver.calc_temperaturesystem_strongbc_solid_init(ts0, tf0, ts, tf);
            //cht_tf0 = tf;
            cht_ts0 = ts;
         }

         nssolver.solve_temperature(v0,cht_tf0,cht_ts0,v,tf,ts,flux);
      //   for (size_t i = 0; i < v0.Size(); i++)
      //{
         //std::cout << " flux ";
         //std::cout <<  v_last[i] << " --- " << v0[i]  << " --- " << v[i] << " \n";
      //}

         if (iter2==1)
         {  
            //cht_dflux_error_norm = 0;
            flux0 = flux;
            for (size_t i = 0; i < flux0.size(); i++)
            {
               flux0[i] = 0;
            }
         }
         
         //cht_tf_error_norm_l2 = std::abs((cht_tf_error_norm_l2 - tf.Norml2())/tf.Norml2());
         //cht_ts_error_norm_l2 = std::abs((cht_ts_error_norm_l2 - ts.Norml2())/ts.Norml2());
         //cht_tf_error_norm = std::abs((cht_tf_error_norm - tf.Norml2()));
         //cht_ts_error_norm = std::abs((cht_ts_error_norm - ts.Norml2()));
         /*
         if (iter2!=1)
         {
            for (size_t i = 0; i < flux.size(); i++)
            {
               //std::cout << " flux ";
               //std::cout <<  flux[i] << " \n";
               if (cht_dflux_error_norm < std::abs(flux0[i] - flux[i]))
               {
                  cht_dflux_error_norm = std::abs(flux0[i] - flux[i]);
               }
            }
         }
         */
         

         cht_tf0_last = cht_tf0;
         cht_ts0_last = cht_ts0; 
         cht_tf0_current = tf;
         cht_ts0_current = ts;

         cht_tf0 = tf;
         cht_ts0 = ts;
         tf_current = tf;
         ts_current = ts;

         cht_tf_error_norm_abs = nssolver.error_norm_abs(cht_tf0_last,cht_tf0_current);
         cht_ts_error_norm_abs = nssolver.error_norm_abs(cht_ts0_last,cht_ts0_current);
         cht_dflux_error_norm_abs = nssolver.error_norm_abs_vector(flux0,flux);
         cht_tf_error_norm_rel = nssolver.error_norm_rel(cht_tf0_last,cht_tf0_current);
         cht_ts_error_norm_rel = nssolver.error_norm_rel(cht_ts0_last,cht_ts0_current);
         cht_dflux_error_norm_rel = nssolver.error_norm_rel_vector(flux0,flux);
         

         std::cout << "CHT tf_error_norm_abs ";
         std::cout <<  cht_tf_error_norm_abs << " \n";
         std::cout << "CHT ts_error_norm_abs ";
         std::cout <<  cht_ts_error_norm_abs << " \n";
         std::cout << "CHT dflux_error_norm_abs ";
         std::cout <<  cht_dflux_error_norm_abs << " \n";
         std::cout << "CHT tf_error_norm_rel ";
         std::cout <<  cht_tf_error_norm_rel << " \n";
         std::cout << "CHT ts_error_norm_rel ";
         std::cout <<  cht_ts_error_norm_rel << " \n";
         std::cout << "CHT dflux_error_norm_rel ";
         std::cout <<  cht_dflux_error_norm_rel << " \n";
         
         std::cout << "CHT ITERATION " + std::to_string(iter2) + " END\n";
      
         if (((cht_tf_error_norm_abs>atol)&&(cht_tf_error_norm_rel>rtol))||((cht_ts_error_norm_abs>atol)&&(cht_ts_error_norm_rel>rtol))||((cht_dflux_error_norm_abs>atol)&&(cht_dflux_error_norm_rel>rtol)))
         {
            std::cout << "TRUE all\n";
         }else{
            std::cout << "FALSE all\n";
         }
         if (((cht_tf_error_norm_abs>atol))||((cht_ts_error_norm_abs>atol))||((cht_dflux_error_norm_abs>atol)))
         {
            std::cout << "TRUE absolut\n";
         }else{
            std::cout << "FALSE absolut\n";
         }
         if (((cht_tf_error_norm_rel>rtol))||((cht_ts_error_norm_rel>rtol))||((cht_dflux_error_norm_rel>rtol)))
         {
            std::cout << "TRUE relativ\n";
         }else{
            std::cout << "FALSE relativ\n";
         }
         if (iter2==100)
         {
            std::cout << " too many Interface Iterations! \n";
            break;
         }
      }

      if (nssolver.visualization)
      {
         char vishost[] = "localhost";
         int  visport   = 19916;
         mfem::socketstream sol_sock(vishost, visport);
         sol_sock.precision(8);
         sol_sock << "solution\n" << *nssolver.mesh_fluid << tf << std::flush;
      }
      
      if (nssolver.visualization)
      {
         char vishost[] = "localhost";
         int  visport   = 19916;
         mfem::socketstream sol_sock(vishost, visport);
         sol_sock.precision(8);
         sol_sock << "solution\n" << *nssolver.mesh_solid << ts << std::flush;
      }

      nssolver.visualization = 0;
      //v = v0;
      
      for (size_t i = 0; i < tfr.Size(); i++)
      {  
         bool skip=false;
         for (size_t ii = 0; ii < nssolver.tempf_ess_tdof_list.Size(); ii++)
         {
            if (i==nssolver.tempf_ess_tdof_list[ii])
            {
               skip=true;
            }
         }
         if (!skip)
         {
            tf0[i] = (1-relaxation)*tfr[i] + relaxation*tf[i];
         }
      }
      for (size_t i = 0; i < tsr.Size(); i++)
      {  
         bool skip=false;
         for (size_t ii = 0; ii < nssolver.temps_ess_tdof_list.Size(); ii++)
         {
            if (i==nssolver.temps_ess_tdof_list[ii])
            {
               skip=true;
            }
         }
         if (!skip)
         {
            ts0[i] = (1-relaxation)*tsr[i] + relaxation*ts[i];
         }
      }
      
      //t0 = t;
      tfr = tf0;
      tsr = ts0;
      
      
      /*
      v_error_norm = std::abs((v_error_norm - v0.Norml2())/v0.Norml2());
      //p_error_norm = std::abs((p_error_norm - p0.Norml2())/p0.Norml2());
      p_error_norm = std::abs((p_error_norm - p0.Norml2()));
      tf_error_norm = std::abs((tf_error_norm - tf0.Norml2())/tf0.Norml2());
      ts_error_norm = std::abs((ts_error_norm - ts0.Norml2())/ts0.Norml2());
      */
      
      //v_error_norm = std::abs((v_error_norm - v0.Norml2()));
      //p_error_norm = std::abs((p_error_norm - p0.Norml2()));
      //tf_error_norm = std::abs((tf_error_norm - tf0.Norml2()));
      //ts_error_norm = std::abs((ts_error_norm - ts0.Norml2()));

      
      //for (size_t i = 0; i < v0.Size(); i++)
      //{
         //std::cout << " flux ";
      //   std::cout <<  v_last[i]  << " --- " << v_current[i] << " \n";
      //}

      v_error_norm_abs = nssolver.error_norm_abs(v_last,v_current);
      p_error_norm_abs = nssolver.error_norm_abs(p_last,p_current);
      v_error_norm_rel = nssolver.error_norm_rel(v_last,v_current);
      p_error_norm_rel = nssolver.error_norm_rel(p_last,p_current);

      tf_error_norm_abs = nssolver.error_norm_abs(tf_last,tf);
      ts_error_norm_abs = nssolver.error_norm_abs(ts_last,ts);
      tf_error_norm_rel = nssolver.error_norm_rel(tf_last,tf);
      ts_error_norm_rel = nssolver.error_norm_rel(ts_last,ts);

      std::cout << " v_error_norm_abs ";
      std::cout << v_error_norm_abs << " \n";
      std::cout << " p_error_norm_abs ";
      std::cout <<  p_error_norm_abs << " \n";
      std::cout << " tf_error_norm_abs ";
      std::cout <<  tf_error_norm_abs << " \n";
      std::cout << " ts_error_norm_abs ";
      std::cout <<  ts_error_norm_abs << " \n";
      
      std::cout << " v_error_norm_rel ";
      std::cout << v_error_norm_rel << " \n";
      std::cout << " p_error_norm_rel ";
      std::cout <<  p_error_norm_rel << " \n";
      std::cout << " tf_error_norm_rel ";
      std::cout <<  tf_error_norm_rel << " \n";
      std::cout << " ts_error_norm_rel ";
      std::cout <<  ts_error_norm_rel << " \n";

      v_last = v_current;
      p_last = p_current;
      tf_last = tf_current;
      ts_last = ts_current;


      std::cout << "NURBS STOKES ITERATION " + std::to_string(iter) + " END ######################################################################\n";

      if (((v_error_norm_abs>atol)&&(v_error_norm_rel>rtol))||((p_error_norm_abs>atol)&&(p_error_norm_rel>rtol))||((tf_error_norm_abs>atol)&&(tf_error_norm_rel>rtol))||((ts_error_norm_abs>atol)&&(ts_error_norm_rel>rtol)))
      {
         std::cout << "TRUE all\n";
      }else{
         std::cout << "FALSE all\n";
      }
      if (((v_error_norm_abs>atol))||((p_error_norm_abs>atol))||((tf_error_norm_abs>atol))||((ts_error_norm_abs>atol)))
      {
         std::cout << "TRUE absolut\n";
      }else{
         std::cout << "FALSE absolut\n";
      }
      if (((v_error_norm_rel>rtol))||((p_error_norm_rel>rtol))||((tf_error_norm_rel>rtol))||((ts_error_norm_rel>rtol)))
      {
         std::cout << "TRUE relativ\n";
      }else{
         std::cout << "FALSE relativ\n";
      }

      if (iter==maxIter2)
      //if (iter==1)
      {
         break;
      }     
   }

   std::cout << "NURBS STOKES END\n";

   // POSTPROCESSING
   int nop = 20; // number of points - 1
   double x_coor = 0.5; // xcoord for extracting results
   double y_coor = 0.5; // xcoord for extracting results
   std::vector<std::vector<double>> post_vector; // vector results postprocessing 
   // SCALARFIELD tf0
   mfem::GridFunction* gf_source = new mfem::GridFunction(tf0);
   mfem::GridFunctionCoefficient gfc_source(gf_source); 
   int sdim = gf_source->FESpace()->GetMesh()->SpaceDimension();

   std::vector<mfem::Vector> phys_points;
   for (size_t i = 0; i < nop + 1; i++)
   {
      mfem::Vector phys_point(sdim);
      phys_point[0] = x_coor;
      phys_point[1] = double(i)/double(nop);
      phys_points.push_back(phys_point);
   }
   
   for (size_t i = 0; i < phys_points.size(); i++)
   {
      int ret;
      mfem::IntegrationPoint ip_source;
      int elem_idx;
      mfem::ElementTransformation* T_source;
      for (int ii=0; ii<gf_source->FESpace()->GetMesh()->GetNE(); ++ii)
      {
         T_source = gf_source->FESpace()->GetMesh()->GetElementTransformation(ii);
         mfem::InverseElementTransformation invtran(T_source);
         ret = invtran.Transform(phys_points[i], ip_source);
         if (ret == 0)
         {
            elem_idx = ii;
            //std::cout << " source " << gf_source->GetValue(elem_idx, ip_source,1) << " element " << i; 
            //std::cout << " ElementNo " << T_source->ElementNo << " ElementType " << T_source->ElementType << " \n";
            //std::cout << " ElementNo " << T.ElementNo << " ElementType " << T.ElementType << " \n";         
            break;
         }
      }
      if(ret != 0){
         std::cout << "phys point not found ret !=0 \n";
         break;
      }
      double source = 0;
      //std::cout << " ref ip_source.x " << ip_source.x << " ip_source.y " << ip_source.y << "\n";
      T_source->TransformBack(phys_points[i], ip_source);
      //std::cout << " ip_source.x " << ip_source.x << " ip_source.y " << ip_source.y << "\n";
      //std::cout << " phys_points[i][0] " << phys_points[i][0] << " phys_points[i][1] " << phys_points[i][1] << "\n";
      T_source->Reset();
      T_source->SetIntPoint(&ip_source);
      source = gfc_source.Eval(*T_source, ip_source);
      std::cout << "SOURCE " + std::to_string(source) + "\n";
      post_vector.push_back({phys_points[i][0],phys_points[i][1],source});
   }
   // VECTORFIELD v0
   //mfem::GridFunction* gf_source = new mfem::GridFunction(v0);
   gf_source = new mfem::GridFunction(v0);
   mfem::VectorGridFunctionCoefficient vgfc_source(gf_source); 
   sdim = gf_source->FESpace()->GetMesh()->SpaceDimension();

   /*std::vector<mfem::Vector> phys_points;
   for (size_t i = 0; i < 11; i++)
   {
      mfem::Vector phys_point(sdim);
      phys_point[0] = 0.5;
      phys_point[1] = double(i)/10;
      phys_points.push_back(phys_point);
   }
   */
   for (size_t i = 0; i < phys_points.size(); i++)
   {
      int ret;
      mfem::IntegrationPoint ip_source;
      int elem_idx;
      mfem::ElementTransformation* T_source;
      for (int ii=0; ii<gf_source->FESpace()->GetMesh()->GetNE(); ++ii)
      {
         T_source = gf_source->FESpace()->GetMesh()->GetElementTransformation(ii);
         mfem::InverseElementTransformation invtran(T_source);
         ret = invtran.Transform(phys_points[i], ip_source);
         if (ret == 0)
         {
            elem_idx = ii;
            //std::cout << " source " << gf_source->GetValue(elem_idx, ip_source,1) << " element " << i; 
            //std::cout << " ElementNo " << T_source->ElementNo << " ElementType " << T_source->ElementType << " \n";
            //std::cout << " ElementNo " << T.ElementNo << " ElementType " << T.ElementType << " \n";         
            break;
         }
      }
      if(ret != 0){
         std::cout << "phys point not found ret !=0 \n";
         break;
      }
      mfem::Vector vsource(sdim);
      vsource[0] = 0.;
      vsource[1] = 0.;
      T_source->TransformBack(phys_points[i], ip_source);
      T_source->Reset();
      T_source->SetIntPoint(&ip_source);
      vgfc_source.Eval(vsource,*T_source, ip_source);
      std::cout << "SOURCE[0] " + std::to_string(vsource[0]) + " SOURCE[1] " + std::to_string(vsource[1]) + "\n";
      post_vector.push_back({phys_points[i][0],phys_points[i][1],vsource[0],vsource[1]});
   }

   std::string filename = "tf0.res";
   std::ofstream output_file;
   output_file.open(filename.c_str(), std::ofstream::out | std::ofstream::trunc);
   output_file << "fluid temperature\n";
   for (size_t i = 0; i < nop + 1; i++)
   {
      for (size_t ii = 0; ii < post_vector[i].size(); ii++)
      {
         std::cout << "post_vector[" + std::to_string(ii)+ "] " + std::to_string(post_vector[i][ii])+ " ";
         output_file << std::to_string(post_vector[i][ii]) << " ";
      }
      std::cout << "\n";
      output_file << "\n";
   }
   output_file.close();

   //std::string filename = "vectorfield.res";
   filename = "v0.res";
   int ic = 0;
   //std::ofstream output_file;
   output_file.open(filename.c_str(), std::ofstream::out | std::ofstream::trunc);
   output_file << "velocity\n";
   for (size_t i = nop + 1; i < 2*nop + 2; i++)
   {
      for (size_t ii = 0; ii < post_vector[i].size(); ii++)
      {
         std::cout << "post_vector[" + std::to_string(ii)+ "] " + std::to_string(post_vector[i][ii])+ " ";
         output_file << std::to_string(post_vector[i][ii]) << " ";
      }
      std::cout << "\n";
      output_file << "\n";
      ic = i;
   }
   output_file.close();

   // FSI TEST
   // POSTPROCESSING
   nop = 20; // number of points - 1
   y_coor = 0.5; // ycoord for extracting results
   
   // SCALARFIELD ts0 part1
   gf_source = new mfem::GridFunction(ts0);
   gfc_source = mfem::GridFunctionCoefficient(gf_source); 
   sdim = gf_source->FESpace()->GetMesh()->SpaceDimension();

   phys_points.clear();
   for (size_t i = 0; i < nop + 1; i++)
   {
      mfem::Vector phys_point(sdim);
      phys_point[0] = -1 + double(i)/double(nop);
      phys_point[1] = y_coor;
      phys_points.push_back(phys_point);
   }
   
   for (size_t i = 0; i < phys_points.size(); i++)
   {
      int ret;
      mfem::IntegrationPoint ip_source;
      int elem_idx;
      mfem::ElementTransformation* T_source;
      for (int ii=0; ii<gf_source->FESpace()->GetMesh()->GetNE(); ++ii)
      {
         T_source = gf_source->FESpace()->GetMesh()->GetElementTransformation(ii);
         mfem::InverseElementTransformation invtran(T_source);
         ret = invtran.Transform(phys_points[i], ip_source);
         if (ret == 0)
         {
            elem_idx = ii;
            break;
         }
      }
      if(ret != 0){
         std::cout << "phys point not found ret !=0 \n";
         break;
      }
      double source = 0;
      T_source->TransformBack(phys_points[i], ip_source);
      T_source->Reset();
      T_source->SetIntPoint(&ip_source);
      source = gfc_source.Eval(*T_source, ip_source);
      std::cout << "SOURCE " + std::to_string(source) + "\n";
      post_vector.push_back({phys_points[i][0],phys_points[i][1],source});
   }

   // SCALARFIELD tf0
   gf_source = new mfem::GridFunction(tf0);
   gfc_source = mfem::GridFunctionCoefficient(gf_source); 
   sdim = gf_source->FESpace()->GetMesh()->SpaceDimension();

   phys_points.clear();
   for (size_t i = 0; i < nop + 1; i++)
   {
      mfem::Vector phys_point(sdim);
      phys_point[0] = double(i)/double(nop);
      phys_point[1] = y_coor;
      phys_points.push_back(phys_point);
   }
   
   for (size_t i = 0; i < phys_points.size(); i++)
   {
      int ret;
      mfem::IntegrationPoint ip_source;
      int elem_idx;
      mfem::ElementTransformation* T_source;
      for (int ii=0; ii<gf_source->FESpace()->GetMesh()->GetNE(); ++ii)
      {
         T_source = gf_source->FESpace()->GetMesh()->GetElementTransformation(ii);
         mfem::InverseElementTransformation invtran(T_source);
         ret = invtran.Transform(phys_points[i], ip_source);
         if (ret == 0)
         {
            elem_idx = ii;
            break;
         }
      }
      if(ret != 0){
         std::cout << "phys point not found ret !=0 \n";
         break;
      }
      double source = 0;
      T_source->TransformBack(phys_points[i], ip_source);
      T_source->Reset();
      T_source->SetIntPoint(&ip_source);
      source = gfc_source.Eval(*T_source, ip_source);
      std::cout << "SOURCE " + std::to_string(source) + "\n";
      post_vector.push_back({phys_points[i][0],phys_points[i][1],source});
   }

   // SCALARFIELD ts0 part2
   gf_source = new mfem::GridFunction(ts0);
   gfc_source = mfem::GridFunctionCoefficient(gf_source); 
   sdim = gf_source->FESpace()->GetMesh()->SpaceDimension();

   phys_points.clear();
   for (size_t i = 0; i < nop + 1; i++)
   {
      mfem::Vector phys_point(sdim);
      phys_point[0] = 1 + double(i)/double(nop);
      phys_point[1] = y_coor;
      phys_points.push_back(phys_point);
   }
   
   for (size_t i = 0; i < phys_points.size(); i++)
   {
      int ret;
      mfem::IntegrationPoint ip_source;
      int elem_idx;
      mfem::ElementTransformation* T_source;
      for (int ii=0; ii<gf_source->FESpace()->GetMesh()->GetNE(); ++ii)
      {
         T_source = gf_source->FESpace()->GetMesh()->GetElementTransformation(ii);
         mfem::InverseElementTransformation invtran(T_source);
         ret = invtran.Transform(phys_points[i], ip_source);
         if (ret == 0)
         {
            elem_idx = ii;
            break;
         }
      }
      if(ret != 0){
         std::cout << "phys point not found ret !=0 \n";
         break;
      }
      double source = 0;
      T_source->TransformBack(phys_points[i], ip_source);
      T_source->Reset();
      T_source->SetIntPoint(&ip_source);
      source = gfc_source.Eval(*T_source, ip_source);
      std::cout << "SOURCE " + std::to_string(source) + "\n";
      post_vector.push_back({phys_points[i][0],phys_points[i][1],source});
   }

   filename = "fsi.res";
   //std::ofstream output_file;
   output_file.open(filename.c_str(), std::ofstream::out | std::ofstream::trunc);
   output_file << "fsi temperature\n";
   for (size_t i = ic + 1; i < post_vector.size(); i++)
   {
      for (size_t ii = 0; ii < post_vector[i].size(); ii++)
      {
         std::cout << "post_vector[" + std::to_string(ii)+ "] " + std::to_string(post_vector[i][ii])+ " ";
         output_file << std::to_string(post_vector[i][ii]) << " ";
      }
      std::cout << "\n";
      output_file << "\n";
   }
   output_file.close();

   return 0;
}