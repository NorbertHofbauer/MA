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
#include "nurbs_stokes_solverp.hpp"
#include "viscosity_modelsp.hpp"

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
   double temp_diffusion_const = 0.1;
   const char *mesh_file = "../../../MA/data/quad_lin_nurbs.mesh";  //our standard test mesh
   int vis_model = 1;         // type of the viscosity model
   // 1 Carreau
   mfem::Vector model_parameters; // array to get the model parameters
   double relaxation = 1;
   int ref_levels = 0;              // standard number of refinements for the mesh
   int par_ref_levels = 0;
   int maxIter = 5000;
   int maxIter2 = 10;
   double rtol = 1.e-6; // convergence criteria
   double atol = 1.e-5; // convergence criteria
   bool visualization = 1;
   int orderelevationvelocity = 0;
   int orderelevationpressure = 0;
   int orderelevationtemperature = 0;
   double v_error_norm_l2 = 1;
   double p_error_norm_l2 = 1;
   double t_error_norm_l2 = 1;
   double max_error = 1e-3;
   mfem::Array<int> vdbc_bdr_noslip;
   mfem::Array<int> vdbc_bdr;
   mfem::Vector vdbc_bdr_values;
   mfem::Array<int> pdbc_bdr;
   mfem::Vector pdbc_bdr_values;
   mfem::Array<int> tdbc_bdr;
   mfem::Vector tdbc_bdr_values;
   double init_temp;
   
   //mpirun -np 4 nurbs_stokesp -m ../../../MA/data/quad_lin_nurbs.mesh -r 3 -vm 0 -mp '205.9'  -vnos '1 3' -vdbc '4' -vdbc_values '28.5 0' -pdbc '2' -pdbc_values '10' -tdbc '1 3 4' -tdbc_values '28 5.5 -10' -d 1 -tdc -0.1 -mi 10000 -mi2 10 -oev 2 -oep 1 -oet 2 -rel 1
   //mpirun -np 4 nurbs_stokesp -m ../../../MA/data/quad_lin_nurbs.mesh -r 3 -vm 1 -mp '1777 0.064 0.73'  -vnos '1 3' -vdbc '4' -vdbc_values '28.5 0' -pdbc '2' -pdbc_values '10' -tdbc '1 3 4' -tdbc_values '28 5.5 -10' -d 1 -tdc -0.1 -mi 10000 -mi2 10 -oev 2 -oep 1 -oet 2 -rel 1
   //mpirun -np 4 nurbs_stokesp -m ../../../MA/data/quad_lin_nurbs.mesh -r 3 -vm 2 -mp '1777 0.064 0.73 200 123'  -vnos '1 3' -vdbc '4' -vdbc_values '28.5 0' -pdbc '2' -pdbc_values '10' -tdbc '1 3 4' -tdbc_values '220 220 220' -d 1 -tdc -0.1 -mi 10000 -mi2 10 -oev 2 -oep 1 -oet 2 -rel 1
   //mpirun -np 4 nurbs_stokesp -m ../../../MA/data/quad_lin_nurbs.mesh -r 3 -vm 3 -mp '2e+4 0.28 -0.025 170 10'  -vnos '1 3' -vdbc '4' -vdbc_values '28.5 0' -pdbc '2' -pdbc_values '10' -tdbc '1 3 4' -tdbc_values '220 220 220' -d 1 -tdc -0.1 -mi 10000 -mi2 10 -oev 2 -oep 1 -oet 2 -rel 1

   // Parse command-line options.
   // input options for our executable
   mfem::OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.",true);
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
   args.AddOption(&orderelevationtemperature, "-oet", "--orderelevationtemperature", 
                  "the order of the elevation for temperature. e.g. -oet 1");
   args.AddOption(&density, "-d", "--density", 
                  "fluid density. e.g. -d 1000",true);
   args.AddOption(&relaxation, "-rel", "--relaxation", 
                  "Relaxation for outer loop. e.g. -rel 0.5");
   args.AddOption(&temp_diffusion_const, "-tdc", "--temp_diffusion_constant", 
                  "temperature diffusion constant. e.g. -tdc 1000",true);
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
   args.AddOption(&tdbc_bdr, "-tdbc", "--temperature_dirchlet_bdr", 
                  "Temperature Dirichlet Boundaries. e.g. -tdbc 1 2 6 5 ...",true);   
   args.AddOption(&tdbc_bdr_values, "-tdbc_values", "--temperature_dirchlet_bdr_values", 
                  "Values for the Temperature Dirichlet Boundaries. e.g. -tdbc_values 1.1 2.9 6.7 5.22 ...",true);
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(std::cout);
      return 1;
   }
   if (myid == 0)
   {
      args.PrintOptions(std::cout);
   }

   NurbsStokesSolver nssolver;
   nssolver.myid = myid;
   nssolver.num_procs = num_procs;
   nssolver.density = density;     // density of the fluid
   nssolver.temp_diffusion_const = temp_diffusion_const; // value for temperature diffusion constant coefficient
   nssolver.ref_levels = ref_levels;
   nssolver.par_ref_levels = par_ref_levels;
   nssolver.set_order_elevation_velocity(orderelevationvelocity);
   nssolver.set_order_elevation_pressure(orderelevationpressure);
   nssolver.set_order_elevation_temperature(orderelevationtemperature);
   nssolver.maxIter = maxIter;

   nssolver.rtol = rtol; // convergence criteria
   nssolver.atol = atol; // convergence criteria
   
   nssolver.set_dirichletbc_velocity_noslip(vdbc_bdr_noslip);
   nssolver.set_dirichletbc_velocity(vdbc_bdr, vdbc_bdr_values);
   nssolver.set_dirichletbc_pressure(pdbc_bdr, pdbc_bdr_values);
   nssolver.set_dirichletbc_temperature(tdbc_bdr, tdbc_bdr_values);

   nssolver.init();
   mfem::ParGridFunction vr(nssolver.vfes),pr(nssolver.pfes),tr(nssolver.tfes);
   mfem::ParGridFunction v0(nssolver.vfes),p0(nssolver.pfes),t0(nssolver.tfes),v(nssolver.vfes),p(nssolver.pfes),t(nssolver.tfes);
   
   nssolver.visualization = 1;
   nssolver.calc_dirichletbc(v0,p0,t0);
   nssolver.visualization = 0;
   //MFEM_ASSERT(false, "TEST STOP");

   vr = v0;
   pr = p0;
   tr = t0;
   
   //for (size_t i = 0; i < 28; i++)
   int iter=0;
   while ((v_error_norm_l2>max_error)||(p_error_norm_l2>max_error)||(t_error_norm_l2>max_error))
   {  
      iter+=1;
      if ((iter==-1)|(iter==maxIter2))
      {
         nssolver.visualization = 1;
      }

      v_error_norm_l2 = v0.Norml2();
      p_error_norm_l2 = p0.Norml2();
      t_error_norm_l2 = t0.Norml2();

      // viscosity model - must be a mfem::coefficient or a child class from coefficient
      if (vis_model==0)
      {
         MFEM_ASSERT(model_parameters.Size() == 1, "1 Parameter needed dynamic viscosity!");
         mfem::ConstantCoefficient kin_vis(model_parameters[0]/density);
         nssolver.solve_flow(v0,p0,t0,v,p,t,kin_vis);
      }else if (vis_model==1)
      {
         MFEM_ASSERT(model_parameters.Size() == 3, "Carreau Model needs 3 Parameters k1,k2,k3!");
         CarreauModelCoefficient kin_vis(model_parameters[0],model_parameters[1],model_parameters[2], density);
         nssolver.solve_flow(v0,p0,t0,v,p,t,kin_vis);
      }else if (vis_model==2)
      {
         MFEM_ASSERT(model_parameters.Size() == 5, "CarreauWLF Model needs 3 Parameters k1,k2,k3,k4,k5!");
         CarreauWLFModelCoefficient kin_vis(model_parameters[0],model_parameters[1],model_parameters[2],model_parameters[3],model_parameters[4], density);
         //set reference temperature on first iteration, otherwise model crashes
         if (iter==1)
         {
            for (size_t i = 0; i < t0.Size(); i++)
            {
               if (t0[i] == 0)
               {
                  t0[i] = model_parameters[3];
               }
               //std::cout <<  t0[i] << " t0 \n";
            }
         }
         nssolver.solve_flow(v0,p0,t0,v,p,t,kin_vis);
      }else if (vis_model==3)
      {
         MFEM_ASSERT(model_parameters.Size() == 5, "PowerLaw Model needs 4 Parameters m0,n,a,T0,shearrate0!");
         PowerLawModelCoefficient kin_vis(model_parameters[0],model_parameters[1],model_parameters[2],model_parameters[3],model_parameters[4], density);
         if (iter==1)
         {
            for (size_t i = 0; i < t0.Size(); i++)
            {
               if (t0[i] == 0)
               {
                  t0[i] = model_parameters[3];
               }
               //std::cout <<  t0[i] << " t0 \n";
            }
         }
         nssolver.solve_flow(v0,p0,t0,v,p,t,kin_vis);
      }else{
         MFEM_ASSERT(false, "Viscosity Model not available!");
      }

      //v0 = (1-relaxation)*vr + relaxation*v;
      //p0 = (1-relaxation)*tr + relaxation*p;
      //v0 = v;
      //p0 = p;

      for (size_t i = 0; i < vr.Size(); i++)
      {
         v0[i] = (1-relaxation)*vr[i] + relaxation*v[i];
      }
      for (size_t i = 0; i < pr.Size(); i++)
      {
         p0[i] = (1-relaxation)*pr[i] + relaxation*p[i];
      }

      vr = v0;
      pr = p0;

      nssolver.visualization = 0;
      if ((iter==-1)|(iter==maxIter2))
      {
         nssolver.visualization = 1;
      }
      nssolver.solve_temperature(v0,t0,v,t);
      nssolver.visualization = 0;
      //v = v0;
      
      for (size_t i = 0; i < tr.Size(); i++)
      {
         t0[i] = (1-relaxation)*tr[i] + relaxation*t[i];
      }
      //t0 = t;
      tr = t0;
      

      v_error_norm_l2 = std::abs(v_error_norm_l2 - v0.Norml2());
      p_error_norm_l2 = std::abs(p_error_norm_l2 - p0.Norml2());
      t_error_norm_l2 = std::abs(t_error_norm_l2 - t0.Norml2());
      
      std::cout << " v_error_norm l2 ";
      std::cout << v_error_norm_l2 << " \n";
      std::cout << " p_error_norm l2 ";
      std::cout <<  p_error_norm_l2 << " \n";
      std::cout << " t_error_norm l2 ";
      std::cout <<  t_error_norm_l2 << " \n";
      
      //std::cout << "v0\n";
      //std::cout << v0 << "\n";
      //std::cout << "p0\n";
      //std::cout << p0 << "\n";
      //std::cout << "t0\n";
      //std::cout << t0 << "\n";
      //std::cout << "v\n";
      //std::cout << v << "\n";
      //std::cout << "p\n";
      //std::cout << p << "\n";
      //std::cout << "t\n";
      //std::cout << t << "\n";
      
      std::cout << "NURBS STOKES ITERATION " + std::to_string(iter) + " END\n";
      if (iter==maxIter2)
      {
         break;
      }
      
   }

   std::cout << "NURBS STOKES END\n";
         
   return 0;
}
