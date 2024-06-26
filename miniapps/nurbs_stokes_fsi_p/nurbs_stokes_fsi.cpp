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

// we don`t want to use namespaces here, otherwise code is harder to unterstand, e.g. mfem::Vector != std::vector
//using namespace std; // bad idea!!!!
//using namespace mfem;

// we build one main in this file, self build classes have to be included later if necessary
int main(int argc, char *argv[])
{
   // Initialize MPI and HYPRE.
   mfem::Mpi::Init(argc, argv);
   int num_procs = mfem::Mpi::WorldSize();
   int myid = mfem::Mpi::WorldRank();
   mfem::Hypre::Init();

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
   int par_ref_levels = 0;
   int maxIter = 5000;
   int maxIter2 = 10;
   double rtol = 1.e-6; // convergence criteria
   double atol = 1.e-5; // convergence criteria
   bool visualization = 1;
   int orderelevationvelocity = 0;
   int orderelevationpressure = 0;
   int orderelevationtemperature_fluid = 0;
   int orderelevationtemperature_solid = 0;
   double v_error_norm = 1;
   double p_error_norm = 1;
   double tf_error_norm = 1;
   double ts_error_norm = 1;
   double cht_tf_error_norm = 1;
   double cht_ts_error_norm = 1;
   double cht_dflux_error_norm = 1;
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
   
   //mpirun -np 4 nurbs_stokes_fsi_p -mf ../../../MA/data/nurbs_fluid_domain.mesh -ms ../../../MA/data/nurbs_solid_domain.mesh -r 2 -vm 3 -mp '2e+4 0.28 -0.025 170 10'  -vnos '1 2 9 6 4 7 10 12' -vdbc '3 8' -vdbc_values '5.5 0 5.5 0' -pdbc '5 11' -pdbc_values '10 10' -tfdbc '3 8' -tfdbc_values '220 220' -tsdbc '1' -tsdbc_values '250' -d 1 -tfdc 0.1 -tsdc 0.5 -tfiface '2 6 9' -tsiface '2 3 4' -mi 10000 -mi2 20 -oev 2 -oep 1 -oetf 0 -oets 0 -rel 0.5 -betaq 0.3 -betat 0.3 -me 1e-8 -at 1e-7 -rt 1e-8 -job test

   // casacade
   //mpirun -np 4 nurbs_stokes_fsi_p -mf ../../../MA/data/fluid.mesh -ms ../../../MA/data/solid.mesh -r 2 -vm 3 -mp '2e+4 0.28 -0.025 170 10'  -vnos '1 3 4 5 12 15 14 16 18 7 9 11 10 13 24 23 26 28 17 20 22 25' -vdbc '2 6 8' -vdbc_values '5.5 0 5.5 0 5.5 0' -pdbc '27 21 19' -pdbc_values '0 0 0' -tfdbc '2 6 8' -tfdbc_values '220 220 220' -tsdbc '1 7' -tsdbc_values '250 250' -d 1 -tfdc 0.1 -tsdc 0.5 -tfiface '12 13' -tsiface '3 5' -mi 10000 -mi2 2 -oev 2 -oep 1 -oetf 0 -oets 0 -rel 0.5 -betaq 0.3 -betat 0.3 -me 1e-8 -at 1e-7 -rt 1e-8  -job test
      
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
   args.AddOption(&par_ref_levels, "-rp", "--refine-parallel",
                  "Number of times to refine the mesh uniformly in parallel.");
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
   args.AddOption(&max_error, "-me", "--maxError", 
                  "convergence criteria: maximum error between iterations in the outer loop. e.g. -me 1.e-3");                  
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
   nssolver.myid = myid;
   nssolver.num_procs = num_procs;
   nssolver.jobname = jobname;
   nssolver.meshfile_fluid = mesh_file_fluid;
   nssolver.meshfile_solid = mesh_file_solid;
   nssolver.density = density;     // density of the fluid
   nssolver.temp_diffusion_const_fluid = temp_diffusion_const_fluid; // value for temperature diffusion constant coefficient
   nssolver.temp_diffusion_const_solid = temp_diffusion_const_solid; // value for temperature diffusion constant coefficient
   nssolver.ref_levels = ref_levels;
   nssolver.par_ref_levels = par_ref_levels;
   nssolver.set_order_elevation_velocity(orderelevationvelocity);
   nssolver.set_order_elevation_pressure(orderelevationpressure);
   nssolver.set_order_elevation_temperature_fluid(orderelevationtemperature_fluid);
   nssolver.set_order_elevation_temperature_solid(orderelevationtemperature_solid);
   nssolver.maxIter = maxIter;

   nssolver.rtol = rtol; // convergence criteria
   nssolver.atol = atol; // convergence criteria
   nssolver.max_error = max_error;
   nssolver.beta_q = 0.3;
   
   nssolver.set_dirichletbc_velocity_noslip(vdbc_bdr_noslip);
   nssolver.set_dirichletbc_velocity(vdbc_bdr, vdbc_bdr_values);
   nssolver.set_dirichletbc_pressure(pdbc_bdr, pdbc_bdr_values);
   nssolver.set_dirichletbc_temperature_fluid(tfdbc_bdr, tfdbc_bdr_values);
   nssolver.set_dirichletbc_temperature_solid(tsdbc_bdr, tsdbc_bdr_values);
   nssolver.set_iface_fluid(tfiface_bdr);
   nssolver.set_iface_solid(tsiface_bdr);

   nssolver.init();

   mfem::ParGridFunction vr(nssolver.vfes),pr(nssolver.pfes),tfr(nssolver.tffes),tsr(nssolver.tsfes);
   mfem::ParGridFunction v0(nssolver.vfes),p0(nssolver.pfes),tf0(nssolver.tffes),v(nssolver.vfes),p(nssolver.pfes),tf(nssolver.tffes),ts0(nssolver.tsfes),ts(nssolver.tsfes);
   mfem::ParGridFunction cht_tf0(nssolver.tffes),cht_ts0(nssolver.tsfes);
   
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
   while ((v_error_norm>max_error)||(p_error_norm>max_error)||(tf_error_norm>max_error)||(ts_error_norm>max_error))
   //while (false)
   {  
      iter+=1;
      nssolver.iter = iter;
      if ((iter==-1)|(iter==maxIter2))
      {
         nssolver.visualization = 1;
      }

      v_error_norm = v0.Norml2();
      p_error_norm = p0.Norml2();
      tf_error_norm = tf0.Norml2();
      ts_error_norm = ts0.Norml2();

      // viscosity model - must be a mfem::coefficient or a child class from coefficient
      if (vis_model==0)
      {
         MFEM_ASSERT(model_parameters.Size() == 1, "1 Parameter needed dynamic viscosity!");
         mfem::ConstantCoefficient kin_vis(model_parameters[0]/density);
         //nssolver.solve_flow(v0,p0,tf0,v,p,tf,kin_vis);
      }else if (vis_model==1)
      {
         MFEM_ASSERT(model_parameters.Size() == 3, "Carreau Model needs 3 Parameters k1,k2,k3!");
         CarreauModelCoefficient kin_vis(model_parameters[0],model_parameters[1],model_parameters[2], density);
         //nssolver.solve_flow(v0,p0,tf0,v,p,tf,kin_vis);
      }else if (vis_model==2)
      {
         MFEM_ASSERT(model_parameters.Size() == 5, "CarreauWLF Model needs 3 Parameters k1,k2,k3,k4,k5!");
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
         //nssolver.solve_flow(v0,p0,tf0,v,p,tf,kin_vis);
      }else if (vis_model==3)
      {
         MFEM_ASSERT(model_parameters.Size() == 5, "PowerLaw Model needs 4 Parameters m0,n,a,T0,shearrate0!");
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
         //nssolver.solve_flow(v0,p0,tf0,v,p,tf,kin_vis);
      }else{
         MFEM_ASSERT(false, "Viscosity Model not available!");
      }
            
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
      
      cht_tf_error_norm = 1;
      cht_ts_error_norm = 1;
      cht_dflux_error_norm = 1;
      flux.clear();

      cht_tf0 = tf0;
      cht_ts0 = ts0;

      int iter2 = 0;
      while ((cht_tf_error_norm>max_error)||(cht_ts_error_norm>max_error)||(cht_dflux_error_norm>max_error))
      //while ((cht_tf_error_norm>max_error)||(cht_ts_error_norm>max_error))
      {  
         iter2+=1;
                  
         cht_tf_error_norm = cht_tf0.Norml2();
         cht_ts_error_norm = cht_ts0.Norml2();

         if (iter2!=1)
         {  
            cht_dflux_error_norm = 0;
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

         cht_tf0 = tf;
         cht_ts0 = ts;

         //cht_tf_error_norm_l2 = std::abs((cht_tf_error_norm_l2 - tf.Norml2())/tf.Norml2());
         //cht_ts_error_norm_l2 = std::abs((cht_ts_error_norm_l2 - ts.Norml2())/ts.Norml2());
         cht_tf_error_norm = std::abs((cht_tf_error_norm - tf.Norml2()));
         cht_ts_error_norm = std::abs((cht_ts_error_norm - ts.Norml2()));

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

         std::cout << "CHT tf_error_norm ";
         std::cout <<  cht_tf_error_norm << " \n";
         std::cout << "CHT ts_error_norm ";
         std::cout <<  cht_ts_error_norm << " \n";
         std::cout << "CHT dflux_error_norm ";
         std::cout <<  cht_dflux_error_norm << " \n";

         if (iter2==2)
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
         sol_sock << "parallel " << num_procs << " " << myid << "\n";
         sol_sock.precision(8);
         sol_sock << "solution\n" << *nssolver.pmesh_fluid << tf << std::flush;
      }
      
      if (nssolver.visualization)
      {
         char vishost[] = "localhost";
         int  visport   = 19916;
         mfem::socketstream sol_sock(vishost, visport);
         sol_sock << "parallel " << num_procs << " " << myid << "\n";
         sol_sock.precision(8);
         sol_sock << "solution\n" << *nssolver.pmesh_solid << ts << std::flush;
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
      
      v_error_norm = std::abs((v_error_norm - v0.Norml2()));
      p_error_norm = std::abs((p_error_norm - p0.Norml2()));
      tf_error_norm = std::abs((tf_error_norm - tf0.Norml2()));
      ts_error_norm = std::abs((ts_error_norm - ts0.Norml2()));
      

      std::cout << " v_error_norm ";
      std::cout << v_error_norm << " \n";
      std::cout << " p_error_norm ";
      std::cout <<  p_error_norm << " \n";
      std::cout << " tf_error_norm ";
      std::cout <<  tf_error_norm << " \n";
      std::cout << " ts_error_norm ";
      std::cout <<  ts_error_norm << " \n";
      
      std::cout << "NURBS STOKES ITERATION " + std::to_string(iter) + " END\n";
      if (iter==maxIter2)
      {
         break;
      }
      
   }

   std::cout << "NURBS STOKES END\n";
         
   return 0;
}
