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
#include "nurbs_stokes_solver.hpp"
#include "viscosity_models.hpp"

// we don`t want to use namespaces here, otherwise code is harder to unterstand, e.g. mfem::Vector != std::vector
//using namespace std; // bad idea!!!!
//using namespace mfem;

// we build one main in this file, self build classes have to be included later if necessary
int main(int argc, char *argv[])
{
   // Setup
   // we define the standard values for our boundaries, equation constants and the meshfile
   double v_max = 28;               // max velocity for our boundary on the inlet
   double p_val = 100;           // value for pressure boundary
   double kin_viscosity = 20000;    // value for kinematic visosity
   double density = 2;    // value for density
   const char *mesh_file = "../../../MA/data/quad_lin_nurbs.mesh";  //our standard test mesh
   int ref_levels = 0;              // standard number of refinements for the mesh
   bool visualization = 1;
   mfem::Array<int> order(4);       // to store order from mesh and order elevation. the order for the finite element collections and spaces will be, velocity=order[0] + 1 + order[1], pressure=order[0] + order[1]
   order[0] = 0;                    // mesh order 
   order[1] = 0;
   order[2] = 0;
   order[3] = 0;
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

   //./nurbs_stokes -m ../../../MA/data/quad_lin_nurbs.mesh -r 3 -vnos '1 3' -vdbc '4' -vdbc_values '28.5 0' -pdbc '2' -pdbc_values '10' -tdbc '1 3 4' -tdbc_values '28 5.5 -10'

   // Parse command-line options.
   // input options for our executable
   mfem::OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&ref_levels, "-r", "--refine",
                  "Number of times to refine the mesh uniformly, -1 for auto.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.AddOption(&vdbc_bdr_noslip, "-vnos", "--velocity_bdr_noslip", 
                  "Velocity Boundaries with no slip condition. e.g. -vnos 1 2 6 5 ...");
   args.AddOption(&vdbc_bdr, "-vdbc", "--velocity_dirchlet_bdr", 
                  "Velocity Dirichlet Boundaries. e.g. -vdbc 1 2 6 5 ...");   
   args.AddOption(&vdbc_bdr_values, "-vdbc_values", "--velocity_dirchlet_bdr_values", 
                  "Vector Values for the Velocity Dirichlet Boundaries. e.g. for 3D you need 3 Values for 1 boundary -vdbc_values 1.1 2.9 6.7 ...");
   args.AddOption(&pdbc_bdr, "-pdbc", "--pressure_dirchlet_bdr", 
                  "Pressure Dirichlet Boundaries. e.g. -pdbc 1 2 6 5 ...");   
   args.AddOption(&pdbc_bdr_values, "-pdbc_values", "--pressure_dirchlet_bdr_values", 
                  "Values for the Pressure Dirichlet Boundaries. e.g. -pdbc_values 1.1 2.9 6.7 5.22 ...");
   args.AddOption(&tdbc_bdr, "-tdbc", "--temperature_dirchlet_bdr", 
                  "Temperature Dirichlet Boundaries. e.g. -tdbc 1 2 6 5 ...");   
   args.AddOption(&tdbc_bdr_values, "-tdbc_values", "--temperature_dirchlet_bdr_values", 
                  "Values for the Temperature Dirichlet Boundaries. e.g. -tdbc_values 1.1 2.9 6.7 5.22 ...");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(std::cout);
      return 1;
   }
   args.PrintOptions(std::cout);

   // viscosity model - must be a mfem::coefficient or a child class from coefficient
   CarreauModelCoefficient kin_vis(6500, 0.13, 0.725, density);
   //mfem::ConstantCoefficient kin_vis(kin_viscosity);
   //kin_vis = new CarreauModelCoefficient();
   //CarreauModelCoefficient kin_vis(1,2,3);

   NurbsStokesSolver nssolver;
   nssolver.density = density;     // density of the fluid
   nssolver.temp_diffusion_const = 0.1; // value for temperature diffusion constant coefficient
   nssolver.ref_levels = ref_levels;
   nssolver.set_order_elevation_velocity(1);
   nssolver.set_order_elevation_pressure(1);
   nssolver.set_order_elevation_temperature(1);
   nssolver.maxIter = 5000;

   nssolver.rtol = 1.e-6; // convergence criteria
   nssolver.atol = 1.e-5; // convergence criteria
   
   nssolver.set_dirichletbc_velocity_noslip(vdbc_bdr_noslip);
   nssolver.set_dirichletbc_velocity(vdbc_bdr, vdbc_bdr_values);
   nssolver.set_dirichletbc_pressure(pdbc_bdr, pdbc_bdr_values);
   nssolver.set_dirichletbc_temperature(tdbc_bdr, tdbc_bdr_values);

   nssolver.init();

   mfem::GridFunction v0(nssolver.vfes),p0(nssolver.pfes),t0(nssolver.tfes),v(nssolver.vfes),p(nssolver.pfes),t(nssolver.tfes);
   nssolver.visualization = 0;
   nssolver.calc_dirichletbc(v0,p0,t0);
   nssolver.visualization = 0;
   
   //for (size_t i = 0; i < 28; i++)
   int iter=0;
   int max_iter=10;
   while ((v_error_norm_l2>max_error)||(p_error_norm_l2>max_error)||(t_error_norm_l2>max_error))
   {  
      iter+=1;
      if ((iter==-1)|(iter==max_iter))
      {
         nssolver.visualization = 0;
      }

      v_error_norm_l2 = v0.Norml2();
      p_error_norm_l2 = p0.Norml2();
      t_error_norm_l2 = t0.Norml2();

      nssolver.solve_flow(v0,p0,t0,v,p,t,kin_vis);
      v0 = v;
      p0 = p;

      nssolver.visualization = 0;
      nssolver.solve_temperature(v0,t0,v,t);
      nssolver.visualization = 0;
      //v = v0;
      t0 = t;

      v_error_norm_l2 = std::abs(v_error_norm_l2 - v0.Norml2());
      p_error_norm_l2 = std::abs(p_error_norm_l2 - p0.Norml2());
      t_error_norm_l2 = std::abs(t_error_norm_l2 - t0.Norml2());
      
      std::cout << "v_error_norm l2 ## p_error_norm l2 ## t_error_norm l2 \n";
      std::cout << v_error_norm_l2 << " ## "<< p_error_norm_l2 << " ## "<< t_error_norm_l2 << " \n";
      
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
      if (iter==max_iter)
      {
         break;
      }
      
   }

   std::cout << "NURBS STOKES END\n";
         
   return 0;
}
