#ifndef NURBSSTOKESSOLVER_HPP
#define NURBSSTOKESSOLVER_HPP
//                          MFEM - NURBS Version
//
// equations can be viewed in overleaf
// here we use strongly imposed dirichlet boundary conditions
// https://www.overleaf.com/project/62a3a05e216fe69ec1a7b6ac
//
//
// TODO: - for building an MPI application we need to integrate the "Par" everywhere!

// we don`t want to use namespaces here, otherwise code is harder to unterstand, e.g. mfem::Vector != std::vector
//using namespace std; // bad idea!!!!
//using namespace mfem;
#include "mfem.hpp" // include mfem project

class NurbsStokesSolver
{
public:
   NurbsStokesSolver();
   ~NurbsStokesSolver();
   
   //visualization
   bool visualization = 0;

   // parameters
   double v_max;               // max velocity for our boundary on the inlet
   double p_val;           // value for pressure boundary
   double kin_viscosity;    // value for kinematic visosity
   double temp_1;               // value for temperature
   double temp_2;               // value for temperature
   double temp_diffusion_const; // value for temperature diffusion constant coefficient

   const char *meshfile; // path to meshfile
   mfem::Mesh *mesh;  // mfem mesh object
   int sdim;                        // dimension in space from the mesh 1D, 2D, 3D
   int ref_levels;                  // mesh refinement levels
   mfem::Array<int> order;       // to store order from mesh and order elevation
   // collection of finite elements from the same family (H1,RT,L2, nurbs,...) in multiple dimensions,
   // this is used to match the dof's of a finite element space between elements
   // https://docs.mfem.org/html/fe__coll_8hpp_source.html
   // we need an NURBSFECollection for every different component in our PDE, e.g.: velocity, pressure, ...
   mfem::FiniteElementCollection *vfec; // velocity
   mfem::FiniteElementCollection *pfec; // pressure
   mfem::FiniteElementCollection *tfec; // temperature
   // we want to use a nurbs mesh, so we need the nurbs extension from mfem
   // for every physical field in our equations we need an extension, here we also declare the order for our nurbs
   mfem::NURBSExtension *vNURBSext = NULL; // velocity
   mfem::NURBSExtension *pNURBSext = NULL; // pressure
   mfem::NURBSExtension *tNURBSext = NULL; // temperature
   // declaration for our finite element spaces
   // here the mesh, NURBSextension, finite element collection and the space dimensions are used to define our desired fe space
   mfem::FiniteElementSpace *vfes; // velocity finite element space, with dimension sdim
   mfem::FiniteElementSpace *pfes;  // pressure finite element space, with dimension 1 (scalar field)
   mfem::FiniteElementSpace *tfes;  // temperature finite element space, with dimension 1 (scalar field)
   // dirichlet boundaries
   mfem::Array<int> vdbc_bdr; // contains the whole boundary markers
   mfem::Array<int> vdbc_bdr_noslip; // to select only the boundary markers for the no slip walls
   mfem::Array<int> vdbc_bdr_inlet;   // to select only the boundary markers for the inlet
   mfem::Array<int> pdbc_bdr; // contains the whole boundary markers
   mfem::Array<int> tdbc_bdr; // contains the whole boundary markers
   mfem::Array<int> tdbc_bdr_inlet;   // to select only the boundary markers for the inlet
   mfem::Array<int> tdbc_bdr_walls;   // to select only the boundary markers for the walls
   mfem::Array<int> vdummy_bdr; // contains the whole boundary markers
   mfem::Array<int> pdummy_bdr; // contains the whole boundary markers
   mfem::Array<int> tdummy_bdr; // contains the whole boundary markers

   mfem::Array<int> vel_ess_tdof_list;
   mfem::Array<int> pres_ess_tdof_list;
   mfem::Array<int> temp_ess_tdof_list;
   mfem::Array<int> vel_ess_tdof_list_dummy;
   mfem::Array<int> pres_ess_tdof_list_dummy;
   mfem::Array<int> temp_ess_tdof_list_dummy;

   bool is_initialized = false;
   bool bcstrong = true;
   bool bcweak = false;
   // SOLVER
   int maxIter=10; // maximal number of iterations
   double rtol = 1.e-10; // convergence criteria
   double atol = 1.e-10; // convergence criteria
   int solver_type = 1;    // type 1  MINRES

   bool init(); // initialize, set mesh properties and elevations before init
   bool update(); // update solver
   bool reset(); // delete all data and initialize afterwards
   bool check_initialized(); // check if object is initialized
   bool use_bcstrong(bool use); // use strong bc
   bool use_bcweak(bool use); // use weak bc
   bool init_dirichletbc(); // inits the boundary markers for the dirichlet bc's
   bool set_meshfile(std::string meshfilepath); //set the path to the meshfile
   bool set_mesh_refinement_level(int refinement_level); //set the refinment level for the mesh refinement
   bool set_order_elevation_velocity(int elevationorder); //set the order elevation for the velocity field
   bool set_order_elevation_pressure(int elevationorder); //set the order elevation for the pressure field
   bool set_order_elevation_temperature(int elevationorder); //set the order elevation for the temperature field
   bool set_parameters(std::vector<double> parameters); // set the parameter values
   bool set_dirichletbc_velocity_noslip(std::vector<int> boundary_marker); //sets the dirichlet boundary in the velocity field for the noslip conditions
   bool set_dirichletbc_velocity_inlet(std::vector<int> boundary_marker); //sets the dirichlet boundary in the velocity field for the inlet
   bool set_dirichletbc_pressure(std::vector<int> boundary_marker); //sets the dirichlet boundary in the pressure field
   bool set_dirichletbc_temperature_inlet(std::vector<int> boundary_marker); //sets the dirichlet boundary in the temperature field for the inlet
   bool set_dirichletbc_temperature_walls(std::vector<int> boundary_marker); //sets the dirichlet boundary in the temperature field for the walls

   bool calc_dirichletbc(mfem::GridFunction &v0, mfem::GridFunction &p0, mfem::GridFunction &t0); // calculate our gridfunction on the dirchlet bc
   
   bool calc_flowsystem_strongbc(mfem::GridFunction &v0, mfem::GridFunction &p0, mfem::GridFunction &t0, mfem::GridFunction &v, mfem::GridFunction &p, mfem::GridFunction &t); // assemble and compute our system matrices with strong boundary conditions
   bool calc_temperaturesystem_strongbc(mfem::GridFunction &v0, mfem::GridFunction &t0, mfem::GridFunction &v, mfem::GridFunction &t); // assemble and compute our system matrices with strong boundary conditions

   bool solve_flow(mfem::GridFunction &v0, mfem::GridFunction &p0, mfem::GridFunction &t0, mfem::GridFunction &v, mfem::GridFunction &p, mfem::GridFunction &t);
   bool solve_temperature(mfem::GridFunction &v0, mfem::GridFunction &t0, mfem::GridFunction &v, mfem::GridFunction &t);
};

#endif // NURBSSTOKESSOLVER_HPP
