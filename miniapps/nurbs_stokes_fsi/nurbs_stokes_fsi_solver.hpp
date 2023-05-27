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
#include "viscosity_models.hpp" // include our shear rate models

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
   double density;    // value for density
   double temp_1;               // value for temperature
   double temp_2;               // value for temperature
   double temp_diffusion_const_fluid; // value for temperature diffusion constant coefficient
   double temp_diffusion_const_solid; // value for temperature diffusion constant coefficient
   // boundaries user data
   mfem::Array<int> user_vdbc_bdr_noslip;
   mfem::Array<int> user_vdbc_bdr;
   mfem::Vector user_vdbc_bdr_values;
   mfem::Array<int> user_pdbc_bdr;
   mfem::Vector user_pdbc_bdr_values;
   mfem::Array<int> user_tfdbc_bdr;
   mfem::Vector user_tfdbc_bdr_values;
   mfem::Array<int> user_tsdbc_bdr;
   mfem::Vector user_tsdbc_bdr_values;

   const char *meshfile_fluid; // path to meshfile fluid
   const char *meshfile_solid; // path to meshfile solid
   mfem::Mesh *mesh_fluid;  // mfem mesh object
   mfem::Mesh *mesh_solid;  // mfem mesh object
   int sdim;                        // dimension in space from the mesh 1D, 2D, 3D
   int ref_levels;                  // mesh refinement levels
   mfem::Array<int> order;       // to store order from mesh and order elevation
   // collection of finite elements from the same family (H1,RT,L2, nurbs,...) in multiple dimensions,
   // this is used to match the dof's of a finite element space between elements
   // https://docs.mfem.org/html/fe__coll_8hpp_source.html
   // we need an NURBSFECollection for every different component in our PDE, e.g.: velocity, pressure, ...
   mfem::FiniteElementCollection *vfec; // velocity
   mfem::FiniteElementCollection *pfec; // pressure
   mfem::FiniteElementCollection *tffec; // temperature fluid
   mfem::FiniteElementCollection *tsfec; // temperature solid
   // we want to use a nurbs mesh, so we need the nurbs extension from mfem
   // for every physical field in our equations we need an extension, here we also declare the order for our nurbs
   mfem::NURBSExtension *vNURBSext = NULL; // velocity
   mfem::NURBSExtension *pNURBSext = NULL; // pressure
   mfem::NURBSExtension *tfNURBSext = NULL; // temperature fluid
   mfem::NURBSExtension *tsNURBSext = NULL; // temperature solid
   // declaration for our finite element spaces
   // here the mesh, NURBSextension, finite element collection and the space dimensions are used to define our desired fe space
   mfem::FiniteElementSpace *vfes; // velocity finite element space, with dimension sdim
   mfem::FiniteElementSpace *pfes;  // pressure finite element space, with dimension 1 (scalar field)
   mfem::FiniteElementSpace *tffes;  // temperature fluid finite element space, with dimension 1 (scalar field)
   mfem::FiniteElementSpace *tsfes;  // temperature fluid finite element space, with dimension 1 (scalar field)
   // dirichlet boundaries
   mfem::Array<int> vdbc_bdr; // contains the set boundary markers
   mfem::Array<int> vdbc_bdr_all; // contains the whole boundary markers
   mfem::Array<int> vdbc_bdr_noslip; // to select only the boundary markers for the no slip walls
   mfem::Array<int> pdbc_bdr; // contains the whole boundary markers
   mfem::Array<int> tfdbc_bdr; // contains the whole boundary markers
   mfem::Array<int> vdummy_bdr; // contains the whole boundary markers
   mfem::Array<int> pdummy_bdr; // contains the whole boundary markers
   mfem::Array<int> tfdummy_bdr; // contains the whole boundary markers
   mfem::Array<int> tsdbc_bdr; // contains the whole boundary markers
   mfem::Array<int> tsdummy_bdr; // contains the whole boundary markers

   mfem::Array<int> vel_ess_tdof_list;
   mfem::Array<int> pres_ess_tdof_list;
   mfem::Array<int> tempf_ess_tdof_list;
   mfem::Array<int> vel_ess_tdof_list_dummy;
   mfem::Array<int> pres_ess_tdof_list_dummy;
   mfem::Array<int> tempf_ess_tdof_list_dummy;
   mfem::Array<int> temps_ess_tdof_list;
   mfem::Array<int> temps_ess_tdof_list_dummy;

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
   bool set_meshfile_fluid(std::string meshfilepath); //set the path to the meshfile
   bool set_meshfile_solid(std::string meshfilepath); //set the path to the meshfile
   bool set_mesh_refinement_level(int refinement_level); //set the refinment level for the mesh refinement
   bool set_order_elevation_velocity(int elevationorder); //set the order elevation for the velocity field
   bool set_order_elevation_pressure(int elevationorder); //set the order elevation for the pressure field
   bool set_order_elevation_temperature_fluid(int elevationorder); //set the order elevation for the temperature field
   bool set_order_elevation_temperature_solid(int elevationorder); //set the order elevation for the temperature field
   bool set_dirichletbc_velocity_noslip(mfem::Array<int> boundaries); //sets the dirichlet boundary in the velocity field for the noslip conditions
   bool set_dirichletbc_velocity(mfem::Array<int> boundaries,mfem::Vector boundaries_values); //sets the dirichlet boundary in the velocity field
   bool set_dirichletbc_pressure(mfem::Array<int> boundaries,mfem::Vector boundaries_values); //sets the dirichlet boundary in the pressure field
   bool set_dirichletbc_temperature_fluid(mfem::Array<int> boundaries,mfem::Vector boundaries_values); //sets the dirichlet boundary in the temperature field
   bool set_dirichletbc_temperature_solid(mfem::Array<int> boundaries,mfem::Vector boundaries_values); //sets the dirichlet boundary in the temperature field

   bool calc_dirichletbc_fluid(mfem::GridFunction &v0, mfem::GridFunction &p0, mfem::GridFunction &tf0); // calculate our gridfunction on the dirchlet bc
   bool calc_dirichletbc_solid(mfem::GridFunction &ts0); // calculate our gridfunction on the dirchlet bc
   
   bool calc_flowsystem_strongbc(mfem::GridFunction &v0, mfem::GridFunction &p0, mfem::GridFunction &tf0, mfem::GridFunction &v, mfem::GridFunction &p, mfem::GridFunction &tf, mfem::Coefficient &kin_vis); // assemble and compute our system matrices with strong boundary conditions
   
   bool calc_temperaturesystem_strongbc_fluid(mfem::GridFunction &v0, mfem::GridFunction &tf0, mfem::GridFunction &v, mfem::GridFunction &tf); // assemble and compute our system matrices with strong boundary conditions
   bool calc_temperaturesystem_strongbc_solid(mfem::GridFunction &ts0, mfem::GridFunction &ts); // assemble and compute our system matrices with strong boundary conditions

   bool solve_flow(mfem::GridFunction &v0, mfem::GridFunction &p0, mfem::GridFunction &tf0, mfem::GridFunction &v, mfem::GridFunction &p, mfem::GridFunction &tf,mfem::Coefficient &kin_vis);
   bool solve_temperature(mfem::GridFunction &v0, mfem::GridFunction &tf0, mfem::GridFunction &ts0, mfem::GridFunction &v, mfem::GridFunction &tf, mfem::GridFunction &ts);
};


#endif // NURBSSTOKESSOLVER_HPP
