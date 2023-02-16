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
   const char *mesh_file = "../../../MA/data/quad_nurbs.mesh";  //our standard test mesh
   int ref_levels = 0;              // standard number of refinements for the mesh
   bool visualization = 1;
   mfem::Array<int> order(4);       // to store order from mesh and order elevation. the order for the finite element collections and spaces will be, velocity=order[0] + 1 + order[1], pressure=order[0] + order[1]
   order[0] = 0;                    // mesh order 
   order[1] = 0;
   order[2] = 0;
   order[3] = 0;

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
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(std::cout);
      return 1;
   }
   args.PrintOptions(std::cout);

   NurbsStokesSolver nssolver;

   nssolver.init();

   mfem::GridFunction v0(nssolver.vfes),p0(nssolver.pfes),t0(nssolver.tfes),v(nssolver.vfes),p(nssolver.pfes),t(nssolver.tfes);

   nssolver.calc_dirichletbc(v0,p0,t0);

   nssolver.solve_flow(v0,p0,t0,v,p,t);
   nssolver.solve_temperature(v0,t0,v,t);
   
   return 0;
}