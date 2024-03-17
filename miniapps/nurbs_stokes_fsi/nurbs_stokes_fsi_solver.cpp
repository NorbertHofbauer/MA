#include "nurbs_stokes_fsi_solver.hpp"
#include <fstream>  // fstream for input and output in textfiles
#include <iostream> // for output to the console

#include <algorithm> // to get access to math functions
#include <cmath>     // to get access to math functions

NurbsStokesSolver::NurbsStokesSolver()
{
   v_max = 28;               // max velocity for our boundary on the inlet
   p_val = 100;           // value for pressure boundary
   kin_viscosity = 200;    // value for kinematic visosity
   density = 1;
   temp_1 = 0;               // value for temperature
   temp_2 = 50;               // value for temperature
   temp_diffusion_const_fluid = 0.1; // value for temperature diffusion constant coefficient
   temp_diffusion_const_solid = 0.5; // value for temperature diffusion constant coefficient
   meshfile_fluid = "../../../MA/data/nurbs_fluid_domain.mesh";  //our standard test mesh
   meshfile_solid = "../../../MA/data/nurbs_solid_domain.mesh";  //our standard test mesh
   ref_levels = 0;
   order = mfem::Array<int>(6);
   order[0] = 0;                    // mesh order fluid
   order[1] = 0;                    // order elevation velocity
   order[2] = 0;                    // order elevation pressure
   order[3] = 0;                    // order elevation temperature fluid
   order[4] = 0;                    // mesh order solid
   order[5] = 0;                    // order elevation temperature solid

}

NurbsStokesSolver::~NurbsStokesSolver()
{}

bool NurbsStokesSolver::init()
{
  if (is_initialized)
  {
      return false; // already initialized
  }else{
      //FLUID
      // Read the mesh from the given mesh file. 
      mesh_fluid = new mfem::Mesh(meshfile_fluid, 1, 1);
      sdim = mesh_fluid->SpaceDimension(); // get dimension in space from the mesh 1D, 2D, 3D
      order[0] = mesh_fluid->NURBSext->GetOrder();
      // mesh refinement
      if (ref_levels < 0) // if ref_level parameter is smaller then 0, set a standard refinement
      {
         ref_levels =
            (int)floor(log(500./mesh_fluid->GetNE())/log(2.)/sdim);
      }
      for (int l = 0; l < ref_levels; l++) // loop for refining of the mesh
      {
         mesh_fluid->UniformRefinement();
      }
      //mesh_fluid->DegreeElevate(1,3);
      mesh_fluid->PrintInfo();
      std::cout << "Using Nurbs mesh for Fluid." << std::endl;
      std::cout << "mesh Space Dimension " << sdim << std::endl;
      std::cout << "given Mesh Order " << order[0] << std::endl;
      std::cout << "given Order elevation velocity " << order[1] << std::endl;
      std::cout << "given Order elevation pressure " << order[2] << std::endl;
      std::cout << "given Order elevation temperature " << order[3] << std::endl;
      std::cout << "number of bdr attributes " << mesh_fluid->bdr_attributes.Max() << std::endl; // output the number of boundary attributes
      // fe collections
      vfec = new mfem::NURBSFECollection(order[0] + order[1]); // Pm+n
      pfec = new mfem::NURBSFECollection(order[0] + order[2]); // Pm+n
      tffec = new mfem::NURBSFECollection(order[0] + order[3]); // Pm+n
      std::cout << "velocity finite element collection Order " << vfec->GetOrder() << std::endl;
      std::cout << "pressure finite element collection Order " << pfec->GetOrder() << std::endl;
      std::cout << "temperature fluid finite element collection Order " << tffec->GetOrder() << std::endl;
      // nurbsextensions
      vNURBSext = new mfem::NURBSExtension(mesh_fluid->NURBSext, order[0]+order[1]);
      pNURBSext = new mfem::NURBSExtension(mesh_fluid->NURBSext, order[0]+order[2]);
      tfNURBSext = new mfem::NURBSExtension(mesh_fluid->NURBSext, order[0]+order[3]);
      std::cout << "velocity nurbs ext Order " << vNURBSext->GetOrder() << std::endl;
      std::cout << "pressure nurbs ext Order " << pNURBSext->GetOrder() << std::endl;
      std::cout << "temperature nurbs ext Order " << tfNURBSext->GetOrder() << std::endl;

      // declaration for our finite element spaces
      vfes = new mfem::FiniteElementSpace(mesh_fluid, vNURBSext, vfec, sdim); // velocity finite element space, with dimension sdim
      pfes = new mfem::FiniteElementSpace(mesh_fluid, pNURBSext, pfec);  // pressure finite element space, with dimension 1 (scalar field)
      tffes = new mfem::FiniteElementSpace(mesh_fluid, tfNURBSext, tffec);  // temperature finite element space, with dimension 1 (scalar field)
      std::cout << "vfes velocity finite element space Order " << vfes->GetMaxElementOrder() << std::endl;
      std::cout << "pfes pressure finite element space Order " << pfes->GetMaxElementOrder() << std::endl;
      std::cout << "tffes temperature finite element space Order " << tffes->GetMaxElementOrder() << std::endl;
      std::cout << "Number of finite element unknowns vfes: "
      << vfes->GetTrueVSize() << std::endl;
      std::cout << "Number of finite element unknowns pfes: "
      << pfes->GetTrueVSize() << std::endl;
      std::cout << "Number of finite element unknowns tffes: "
      << tffes->GetTrueVSize() << std::endl;

      // SOLID
      // Read the mesh from the given mesh file. 
      mesh_solid = new mfem::Mesh(meshfile_solid, 1, 1);
      sdim = mesh_solid->SpaceDimension(); // get dimension in space from the mesh 1D, 2D, 3D
      order[4] = mesh_solid->NURBSext->GetOrder();
      // mesh refinement
      if (ref_levels < 0) // if ref_level parameter is smaller then 0, set a standard refinement
      {
         ref_levels =
            (int)floor(log(500./mesh_solid->GetNE())/log(2.)/sdim);
      }
      for (int l = 0; l < ref_levels; l++) // loop for refining of the mesh
      {
         mesh_solid->UniformRefinement();
      }
      mesh_solid->PrintInfo();
      std::cout << "Using Nurbs mesh for Solid." << std::endl;
      std::cout << "mesh Space Dimension " << sdim << std::endl;
      std::cout << "given Mesh Order " << order[4] << std::endl;
      std::cout << "given Order elevation temperature " << order[5] << std::endl;
      std::cout << "number of bdr attributes " << mesh_solid->bdr_attributes.Max() << std::endl; // output the number of boundary attributes
      // fe collections
      tsfec = new mfem::NURBSFECollection(order[4] + order[5]); // Pm+n
      std::cout << "temperature solid finite element collection Order " << tsfec->GetOrder() << std::endl;
      // nurbsextensions
      tsNURBSext = new mfem::NURBSExtension(mesh_solid->NURBSext, order[4]+order[5]);

      // declaration for our finite element spaces
      tsfes = new mfem::FiniteElementSpace(mesh_solid, tsNURBSext, tsfec);  // temperature finite element space, with dimension 1 (scalar field)
      std::cout << "tsfes temperature finite element space Order " << tsfes->GetMaxElementOrder() << std::endl;
      std::cout << "Number of finite element unknowns tsfes: "
      << tsfes->GetTrueVSize() << std::endl;

      // init the dirichlet bc
      this->init_dirichletbc();

      is_initialized = true;  
      return true;
  }
}

bool NurbsStokesSolver::update()
{ 
  return true;
}

bool NurbsStokesSolver::reset()
{
  //data.clear();
  init();
  return true;
}

bool NurbsStokesSolver::check_initialized()
{
  return is_initialized;
}

bool NurbsStokesSolver::use_bcstrong(bool use)
{
   bcstrong = use;
   bcweak = false;
   return true;
}

bool NurbsStokesSolver::use_bcweak(bool use)
{
   bcweak = use;
   bcstrong = false;
   return true;
}

bool NurbsStokesSolver::init_dirichletbc()
{   
   MFEM_ASSERT(user_vdbc_bdr.Size()*sdim == user_vdbc_bdr_values.Size(), "number of given velocity boundaries and values does'nt match!");
   MFEM_ASSERT(user_pdbc_bdr.Size() == user_pdbc_bdr_values.Size(), "number of given pressure boundaries and values does'nt match!");
   MFEM_ASSERT(user_tfdbc_bdr.Size() == user_tfdbc_bdr_values.Size(), "number of given fluid temperature boundaries and values does'nt match!");
   MFEM_ASSERT(user_tsdbc_bdr.Size() == user_tsdbc_bdr_values.Size(), "number of given solid temperature boundaries and values does'nt match!");
   
   // init our boundary markers
   vdbc_bdr_all = mfem::Array<int>(mesh_fluid->bdr_attributes.Max()); // contains the whole boundary markers
   vdbc_bdr = mfem::Array<int>(mesh_fluid->bdr_attributes.Max());   // to select only the set boundary markers
   vdbc_bdr_noslip = mfem::Array<int>(mesh_fluid->bdr_attributes.Max()); // to select only the boundary markers for the no slip walls
   pdbc_bdr = mfem::Array<int>(mesh_fluid->bdr_attributes.Max()); // contains the whole boundary markers
   tfdbc_bdr = mfem::Array<int>(mesh_fluid->bdr_attributes.Max()); // contains the whole boundary markers
   tfdbc_bdr_all = mfem::Array<int>(mesh_fluid->bdr_attributes.Max()); // contains the whole boundary markers
   tfiface_bdr = mfem::Array<int>(mesh_fluid->bdr_attributes.Max()); // contains the whole boundary markers
   tsdbc_bdr = mfem::Array<int>(mesh_solid->bdr_attributes.Max()); // contains the whole boundary markers
   tsiface_bdr = mfem::Array<int>(mesh_solid->bdr_attributes.Max()); // contains the whole boundary markers
   tsdbc_bdr_all = mfem::Array<int>(mesh_solid->bdr_attributes.Max()); // contains the whole boundary markers
   vdummy_bdr = mfem::Array<int>(mesh_fluid->bdr_attributes.Max()); // contains the whole boundary markers
   pdummy_bdr = mfem::Array<int>(mesh_fluid->bdr_attributes.Max()); // contains the whole boundary markers
   tfdummy_bdr = mfem::Array<int>(mesh_fluid->bdr_attributes.Max()); // contains the whole boundary markers
   tsdummy_bdr = mfem::Array<int>(mesh_solid->bdr_attributes.Max()); // contains the whole boundary markers

   //set marker
   // velocity
   vdbc_bdr = 0;
   vdbc_bdr_all = 0;
   vdbc_bdr_noslip = 0;

   // pressure
   pdbc_bdr = 0;

   // temperature
   tfdbc_bdr = 0;
   tfiface_bdr = 0;
   tfdbc_bdr_all = 0;
   tsdbc_bdr = 0;
   tsiface_bdr = 0;
   tsdbc_bdr_all = 0;

   // mark dummy
   vdummy_bdr = 0;
   pdummy_bdr = 0;
   tfdummy_bdr = 0;
   tsdummy_bdr = 0;

   // get the true dofs of the boundaries
   vfes->GetEssentialTrueDofs(vdbc_bdr, vel_ess_tdof_list);
   pfes->GetEssentialTrueDofs(pdbc_bdr, pres_ess_tdof_list);   
   tffes->GetEssentialTrueDofs(tfdbc_bdr, tempf_ess_tdof_list);
   tsfes->GetEssentialTrueDofs(tsdbc_bdr, temps_ess_tdof_list);
   vfes->GetEssentialTrueDofs(vdummy_bdr, vel_ess_tdof_list_dummy);
   pfes->GetEssentialTrueDofs(pdummy_bdr, pres_ess_tdof_list_dummy);   
   tffes->GetEssentialTrueDofs(tfdummy_bdr, tempf_ess_tdof_list_dummy);
   tsfes->GetEssentialTrueDofs(tsdummy_bdr, temps_ess_tdof_list_dummy);

   for (size_t i = 0; i < user_vdbc_bdr_noslip.Size(); i++)
   {
      vdbc_bdr_all[user_vdbc_bdr_noslip[i]-1] = 1;
      vdbc_bdr_noslip[user_vdbc_bdr_noslip[i]-1] = 1;
      //std::cout << user_vdbc_bdr_noslip[i]-1 << " set velocity no slip \n";
   }
   
   for (size_t i = 0; i < user_vdbc_bdr.Size(); i++)
   {
      vdbc_bdr_all[user_vdbc_bdr[i]-1] = 1;
      vdbc_bdr[user_vdbc_bdr[i]-1] = 1;
   }
   for (size_t i = 0; i < user_pdbc_bdr.Size(); i++)
   {
      pdbc_bdr[user_pdbc_bdr[i]-1] = 1;
   }
   for (size_t i = 0; i < user_tfdbc_bdr.Size(); i++)
   {
      tfdbc_bdr[user_tfdbc_bdr[i]-1] = 1;
      tfdbc_bdr_all[user_tfdbc_bdr[i]-1] = 1;
   }
   for (size_t i = 0; i < user_tfiface_bdr.Size(); i++)
   {
      tfiface_bdr[user_tfiface_bdr[i]-1] = 1;
      tfdbc_bdr_all[user_tfiface_bdr[i]-1] = 1;
   }

   for (size_t i = 0; i < user_tsdbc_bdr.Size(); i++)
   {
      tsdbc_bdr[user_tsdbc_bdr[i]-1] = 1;
      tsdbc_bdr_all[user_tsdbc_bdr[i]-1] = 1;
   }
   for (size_t i = 0; i < user_tsiface_bdr.Size(); i++)
   {
      tsiface_bdr[user_tsiface_bdr[i]-1] = 1;
      tsdbc_bdr_all[user_tsiface_bdr[i]-1] = 1;
   }

   // get the true dofs of the boundaries
   vfes->GetEssentialTrueDofs(vdbc_bdr_all, vel_ess_tdof_list);
   pfes->GetEssentialTrueDofs(pdbc_bdr, pres_ess_tdof_list);   
   tffes->GetEssentialTrueDofs(tfdbc_bdr_all, tempf_ess_tdof_list);
   tsfes->GetEssentialTrueDofs(tsdbc_bdr, temps_ess_tdof_list);

   return true;
}

bool NurbsStokesSolver::set_meshfile_fluid(std::string meshfilepath)
{
   meshfile_fluid = meshfilepath.c_str();
   return true;
}

bool NurbsStokesSolver::set_meshfile_solid(std::string meshfilepath)
{
   meshfile_solid = meshfilepath.c_str();
   return true;
}

bool NurbsStokesSolver::set_mesh_refinement_level(int refinement_level)
{
   ref_levels = refinement_level;
   return true;
}

bool NurbsStokesSolver::set_order_elevation_velocity(int elevationorder)
{
   order[1] = elevationorder;
   return true;
}

bool NurbsStokesSolver::set_order_elevation_pressure(int elevationorder)
{
   order[2] = elevationorder;
   return true;
}

bool NurbsStokesSolver::set_order_elevation_temperature_fluid(int elevationorder)
{
   order[3] = elevationorder;
   return true;
}

bool NurbsStokesSolver::set_order_elevation_temperature_solid(int elevationorder)
{
   order[5] = elevationorder;
   return true;
}

bool NurbsStokesSolver::set_dirichletbc_velocity_noslip(mfem::Array<int> boundaries)
{
   user_vdbc_bdr_noslip = boundaries;
   
   return true;
}

bool NurbsStokesSolver::set_dirichletbc_velocity(mfem::Array<int> boundaries,mfem::Vector boundaries_values)
{
   user_vdbc_bdr = boundaries;
   user_vdbc_bdr_values = boundaries_values;
   
   return true;
}

bool NurbsStokesSolver::set_dirichletbc_pressure(mfem::Array<int> boundaries,mfem::Vector boundaries_values)
{
   user_pdbc_bdr = boundaries;
   user_pdbc_bdr_values = boundaries_values;
   
   return true;
}
   
bool NurbsStokesSolver::set_dirichletbc_temperature_fluid(mfem::Array<int> boundaries,mfem::Vector boundaries_values)
{
   user_tfdbc_bdr = boundaries;
   user_tfdbc_bdr_values = boundaries_values;

   return true;
}

bool NurbsStokesSolver::set_dirichletbc_temperature_solid(mfem::Array<int> boundaries,mfem::Vector boundaries_values)
{
   user_tsdbc_bdr = boundaries;
   user_tsdbc_bdr_values = boundaries_values;

   return true;
}

bool NurbsStokesSolver::set_iface_fluid(mfem::Array<int> boundaries)
{
   user_tfiface_bdr = boundaries;
   
   return true;
}

bool NurbsStokesSolver::set_iface_solid(mfem::Array<int> boundaries)
{
   user_tsiface_bdr = boundaries;
   
   return true;
}

bool NurbsStokesSolver::calc_dirichletbc_fluid(mfem::GridFunction &v0, mfem::GridFunction &p0, mfem::GridFunction &tf0)
{
   auto lambda_inlet = [this](const mfem::Vector &QuadraturPointPosition, mfem::Vector &VelocityValue) -> void
   {
      double h=1;
      VelocityValue[0] = 4*(h-QuadraturPointPosition(1))*QuadraturPointPosition(1)/(h*h)*28;
      //VelocityValue[0] = v_max;
      //VelocityValue[1] = 4*(h-QuadraturPointPosition(1))*QuadraturPointPosition(1)/(h*h)*v_max;
      VelocityValue[1] = 0;
      //std::cout << " qp(0) = " << QuadraturPointPosition(0) << " qp(1) = " << QuadraturPointPosition(1) << std::endl;
      //std::cout << " v(0) = " << VelocityValue(0) << " v(1) = " << VelocityValue(1) << std::endl;      
      return;
   };

   auto lambda_inlet2 = [this](const mfem::Vector &QuadraturPointPosition, double &VelocityValue) -> void
   {
      double h=1;
      VelocityValue = 4*(h-QuadraturPointPosition(1))*QuadraturPointPosition(1)/(h*h)*v_max;
      //VelocityValue[0] = v_max;
      //VelocityValue[1] = 4*(h-QuadraturPointPosition(1))*QuadraturPointPosition(1)/(h*h)*v_max;
      //VelocityValue[1] = 0;
      //std::cout << " qp(0) = " << QuadraturPointPosition(0) << " qp(1) = " << QuadraturPointPosition(1) << std::endl;
      //std::cout << " v(0) = " << VelocityValue(0) << " v(1) = " << VelocityValue(1) << std::endl;      
      return;
   };

   mfem::GridFunction v_bc(vfes),p_bc(pfes), tf_bc(tffes); // to calculate our gridfunction on the dirichlet boundary

   // we need grid functions to first compute the controlpoint values on the boundary, so we can project them on to our system
   // means we will build a system that needed to be solved for the desired boundary values

   // create vectors for coefficients and boundary markers

   std::vector<mfem::Array<int>> vdbc_bdr_marker;
   std::vector<mfem::Vector> vdbc_bdr_vector;
   std::vector<mfem::VectorConstantCoefficient> vdbc_bdr_vectorcoefficient;

   for (size_t i = 0; i < user_vdbc_bdr.Size(); i++)
   {
      vdbc_bdr_marker.push_back(mfem::Array<int>(mesh_fluid->bdr_attributes.Max()));
      vdbc_bdr_marker[i] = 0;
      vdbc_bdr_marker[i][user_vdbc_bdr[i]-1] = 1;

      mfem::Vector vvector(sdim);
      vdbc_bdr_vector.push_back(vvector);
      vdbc_bdr_vector[i] = 0.;
      for (size_t ii = 0; ii < sdim; ii++)
      {
         vdbc_bdr_vector[i][ii] = user_vdbc_bdr_values[i*sdim+ii];
         //std::cout << user_vdbc_bdr_values[i*sdim+ii] << " velocity values \n";
      }
      mfem::VectorConstantCoefficient vcc(vdbc_bdr_vector[i]);
      vdbc_bdr_vectorcoefficient.push_back(vcc);
   }
   
   std::vector<mfem::Array<int>> pdbc_bdr_marker;
   std::vector<mfem::ConstantCoefficient> pdbc_bdr_coefficient;

   for (size_t i = 0; i < user_pdbc_bdr.Size(); i++)
   {
      pdbc_bdr_marker.push_back(mfem::Array<int>(mesh_fluid->bdr_attributes.Max()));
      pdbc_bdr_marker[i] = 0;
      pdbc_bdr_marker[i][user_pdbc_bdr[i]-1] = 1;

      mfem::ConstantCoefficient cc(user_pdbc_bdr_values[i]);
      pdbc_bdr_coefficient.push_back(cc);
      //std::cout << user_pdbc_bdr_values[i] << " pressure values \n";
   }

   std::vector<mfem::Array<int>> tfdbc_bdr_marker;
   std::vector<mfem::ConstantCoefficient> tfdbc_bdr_coefficient;

   for (size_t i = 0; i < user_tfdbc_bdr.Size(); i++)
   {
      tfdbc_bdr_marker.push_back(mfem::Array<int>(mesh_fluid->bdr_attributes.Max()));
      tfdbc_bdr_marker[i] = 0;
      tfdbc_bdr_marker[i][user_tfdbc_bdr[i]-1] = 1;

      mfem::ConstantCoefficient cc(user_tfdbc_bdr_values[i]);
      tfdbc_bdr_coefficient.push_back(cc);
      //std::cout << user_tdbc_bdr_values[i] << " temperature values \n";
   }

   // VELOCITY
   // define rhs with the desired boundary condition values
   mfem::VectorFunctionCoefficient vfc_inlet(sdim, lambda_inlet); // function for our desired boundary condition
   //mfem::FunctionCoefficient vfc_inlet(lambda_inlet2); // function for our desired boundary condition
   mfem::LinearForm *f_bc(new mfem::LinearForm(vfes)); // define linear form for rhs
   // define bilinear form add the boundary, means the nurbs add the boundary
   mfem::BilinearForm a_bc(vfes); // define the bilinear form results in n x n matrix, we use the velocity finite element space
   mfem::ConstantCoefficient One_bc(1); // coefficient for the kinematic viscosity
   //f_bc->AddBoundaryIntegrator(new mfem::VectorBoundaryLFIntegrator(vfc_inlet),vdbc_bdr); // define integrator on desired boundary
   for (size_t i = 0; i < vdbc_bdr_marker.size(); i++)
   {
      f_bc->AddBoundaryIntegrator(new mfem::VectorBoundaryLFIntegrator(vdbc_bdr_vectorcoefficient[i]),vdbc_bdr_marker[i]); // define integrator on desired boundary
      //f_bc->AddBoundaryIntegrator(new mfem::VectorBoundaryLFIntegrator(vfc_inlet),vdbc_bdr_marker[i]); // define integrator on desired boundary
      //std::cout << user_vdbc_bdr[i] << " velocity boundary \n";
      a_bc.AddBoundaryIntegrator(new mfem::VectorMassIntegrator(One_bc),vdbc_bdr_marker[i]); // bilinear form (lambda*u_vector),(v_vector))
   }
   //f_bc->AddBoundaryIntegrator(new mfem::BoundaryNormalLFIntegrator(vfc_inlet),vdbc_bdr_inlet); // define integrator on desired boundary
   f_bc->Assemble(); // assemble the linear form (vector)
   a_bc.Assemble(); // assemble the bilinear form (matrix)
   
   // PRESSURE
   // define rhs with the desired boundary condition values
   //mfem::ConstantCoefficient pfc_outlet(p_val); // function for our desired boundary condition
   mfem::LinearForm *g_bc(new mfem::LinearForm(pfes)); // define linear form for rhs
   // define bilinear form add the boundary, means the nurbs add the boundary
   mfem::BilinearForm b_bc(pfes); // define the bilinear form results in n x n matrix, we use the pressure finite element space
   //g_bc->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(pfc_outlet),pdbc_bdr); // define integrator on desired boundary
   for (size_t i = 0; i < pdbc_bdr_marker.size(); i++)
   {
      g_bc->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(pdbc_bdr_coefficient[i]),pdbc_bdr_marker[i]); // define integrator on desired boundary
      b_bc.AddBoundaryIntegrator(new mfem::MassIntegrator(One_bc),pdbc_bdr_marker[i]); // bilinear form (lambda*u_vector),(v_vector))
   }
   g_bc->Assemble(); // assemble the linear form (vector)
   b_bc.Assemble(); // assemble the bilinear form (matrix)

   // TEMPERATURE
   // define rhs with the desired boundary condition values
   //mfem::ConstantCoefficient tfc_inlet(temp_1); // function for our desired boundary condition
   //mfem::ConstantCoefficient tfc_walls(temp_2); // function for our desired boundary condition
   mfem::LinearForm *h_bc(new mfem::LinearForm(tffes)); // define linear form for rhs
   // define bilinear form add the boundary, means the nurbs add the boundary
   mfem::BilinearForm d_bc(tffes); // define the bilinear form results in n x n matrix, we use the velocity finite element space
   //h_bc->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(tfc_inlet),tdbc_bdr); // define integrator on desired boundary
   //h_bc->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(tfc_walls),tdbc_bdr); // define integrator on desired boundary
   for (size_t i = 0; i < tfdbc_bdr_marker.size(); i++)
   {
      h_bc->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(tfdbc_bdr_coefficient[i]),tfdbc_bdr_marker[i]); // define integrator on desired boundary
      d_bc.AddBoundaryIntegrator(new mfem::MassIntegrator(One_bc),tfdbc_bdr_marker[i]); // bilinear form (lambda*u_vector),(v_vector))
   }
   h_bc->Assemble(); // assemble the linear form (vector)   
   d_bc.Assemble(); // assemble the bilinear form (matrix)

   mfem::SparseMatrix A_BC, B_BC, D_BC;
   mfem::Vector V_BC, F_BC;
   mfem::Vector P_BC, G_BC;
   mfem::Vector T_BC, H_BC;
   a_bc.SetDiagonalPolicy(mfem::Matrix::DiagonalPolicy::DIAG_ZERO); // important, otherwise a different policy will be used, which results in false building of our matrix
   a_bc.FormLinearSystem(vel_ess_tdof_list_dummy, v_bc, *f_bc, A_BC, V_BC, F_BC); // form A_BC
   b_bc.SetDiagonalPolicy(mfem::Matrix::DiagonalPolicy::DIAG_ZERO); // important, otherwise a different policy will be used, which results in false building of our matrix
   b_bc.FormLinearSystem(pres_ess_tdof_list_dummy, p_bc, *g_bc, B_BC, P_BC, G_BC); // form B_BC
   d_bc.SetDiagonalPolicy(mfem::Matrix::DiagonalPolicy::DIAG_ZERO); // important, otherwise a different policy will be used, which results in false building of our matrix
   d_bc.FormLinearSystem(tempf_ess_tdof_list_dummy, tf_bc, *h_bc, D_BC, T_BC, H_BC); // form D_BC
    

   /*
   for (size_t i = 0; i < F_BC.Size(); i++)
   {
      std::cout << F_BC[i] << " F_BC \n";
   }

   for (size_t i = 0; i < G_BC.Size(); i++)
   {
      std::cout << G_BC[i] << " G_BC \n";
   }

   for (size_t i = 0; i < H_BC.Size(); i++)
   {
      std::cout << H_BC[i] << " H_BC \n";
   }
   *//*
   for (size_t i = 0; i < vdbc_bdr_marker.size(); i++)
   {
      std::cout << user_vdbc_bdr[i] << " user_vdbc_bdr \n";
      for (size_t ii = 0; ii < vdbc_bdr_marker[i].Size(); ii++)
      {
         std::cout << vdbc_bdr_marker[i][ii] << " boundary " << ii << " vdbc_bdr_marker \n";
      }
   }
   for (size_t i = 0; i < pdbc_bdr_marker.size(); i++)
   {
      std::cout << user_pdbc_bdr[i] << " user_pdbc_bdr \n";
      for (size_t ii = 0; ii < pdbc_bdr_marker[i].Size(); ii++)
      {
         std::cout << pdbc_bdr_marker[i][ii] << " boundary " << ii << " pdbc_bdr_marker \n";
      }
   }
   for (size_t i = 0; i < tdbc_bdr_marker.size(); i++)
   {
      std::cout << user_tdbc_bdr[i] << " user_tdbc_bdr \n";
      for (size_t ii = 0; ii < tdbc_bdr_marker[i].Size(); ii++)
      {
         std::cout << tdbc_bdr_marker[i][ii] << " boundary " << ii << " tdbc_bdr_marker \n";
      }
   }
   */

   //A_BC.PrintMatlab();
   //B_BC.PrintMatlab();
   //D_BC.PrintMatlab();

   // SOLVER
   // setup solver
   //int maxIter(100); // maximal number of iterations
   //double rtol(1.e-12); // convergence criteria
   //double atol(1.e-12); // convergence criteria

   // setup solver
   //1
   //mfem::MINRESSolver solver;
   //2
   //mfem::CGSolver bc_solver;
   //mfem::Solver *bc_prec;
   //bc_prec = new mfem::GSSmoother();
   //bc_solver.SetPreconditioner(*bc_prec);
   mfem::GMRESSolver bc_solver;
   
   bc_solver.SetAbsTol(atol);
   bc_solver.SetRelTol(rtol);
   bc_solver.SetMaxIter(maxIter);
   bc_solver.SetOperator(A_BC);
   bc_solver.SetKDim((int)maxIter/5);
   bc_solver.SetPrintLevel(3);

   // solve the system
   bc_solver.Mult(F_BC, V_BC);

   // check if solver converged
   if (bc_solver.GetConverged())
   {
      std::cout << "Solver converged in " << bc_solver.GetNumIterations()
                << " iterations with a residual norm of "
                << bc_solver.GetFinalNorm() << ".\n";
   }
   else
   {
      std::cout << "Solver did not converge in " << bc_solver.GetNumIterations()
                << " iterations. Residual norm is " << bc_solver.GetFinalNorm()
                << ".\n";
   }
   //std::cout << "MINRES solver took " << chrono.RealTime() << "s.\n";

   // Save the mesh and the solution
   {
      //std::ofstream mesh_ofs(jobname + "_" + std::string(meshfile_fluid));
      std::ofstream mesh_ofs(jobname + "_mesh_v");
      mesh_ofs.precision(8);
      mesh_fluid->Print(mesh_ofs);
      //vfes->GetMesh()->DegreeElevate(1,3);
      //vfes->GetMesh()->Print(mesh_ofs);

      std::ofstream v_bc_ofs(jobname + "_" +"sol_v_bc.gf");
      v_bc_ofs.precision(8);
      v_bc.Save(v_bc_ofs);
      v0.Save(v_bc_ofs);

      //std::ofstream p_ofs("sol_p.gf");
      //p_ofs.precision(8);
      //p.Save(p_ofs);
   }

   //std::cout << "v_bc " << v_bc << "\n";

   bc_solver.SetOperator(B_BC);

   // solve the system
   bc_solver.Mult(G_BC, P_BC);

   // check if solver converged
   if (bc_solver.GetConverged())
   {
      std::cout << "Solver converged in " << bc_solver.GetNumIterations()
                << " iterations with a residual norm of "
                << bc_solver.GetFinalNorm() << ".\n";
   }
   else
   {
      std::cout << "Solver did not converge in " << bc_solver.GetNumIterations()
                << " iterations. Residual norm is " << bc_solver.GetFinalNorm()
                << ".\n";
   }
   //std::cout << "MINRES solver took " << chrono.RealTime() << "s.\n";

   // Save the mesh and the solution
   {
      //std::ofstream mesh_ofs(jobname + "_" + std::string(meshfile_fluid));
      std::ofstream mesh_ofs(jobname + "_mesh_p");
      mesh_ofs.precision(8);
      mesh_fluid->Print(mesh_ofs);

      std::ofstream p_bc_ofs(jobname + "_" +"sol_p_bc.gf");
      p_bc_ofs.precision(8);
      p_bc.Save(p_bc_ofs);
      p0.Save(p_bc_ofs);
   }

   bc_solver.SetOperator(D_BC);

   // solve the system
   bc_solver.Mult(H_BC, T_BC);

   // check if solver converged
   if (bc_solver.GetConverged())
   {
      std::cout << "Solver converged in " << bc_solver.GetNumIterations()
                << " iterations with a residual norm of "
                << bc_solver.GetFinalNorm() << ".\n";
   }
   else
   {
      std::cout << "Solver did not converge in " << bc_solver.GetNumIterations()
                << " iterations. Residual norm is " << bc_solver.GetFinalNorm()
                << ".\n";
   }
   //std::cout << "MINRES solver took " << chrono.RealTime() << "s.\n";

   // Save the mesh and the solution
   {
      //std::ofstream mesh_ofs(jobname + "_" + std::string(meshfile_fluid));
      std::ofstream mesh_ofs(jobname + "_mesh_tf");
      mesh_ofs.precision(8);
      mesh_fluid->Print(mesh_ofs);

      std::ofstream tf_bc_ofs(jobname + "_" + "sol_tf_bc.gf");
      tf_bc_ofs.precision(8);
      tf_bc.Save(tf_bc_ofs);
      tf0.Save(tf_bc_ofs);
   }
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      mfem::socketstream sol_sock1(vishost, visport);
      sol_sock1.precision(8);
      sol_sock1 << "solution\n" << *mesh_fluid << v_bc << std::flush;
      mfem::socketstream sol_sock2(vishost, visport);
      sol_sock2.precision(8);
      sol_sock2 << "solution\n" << *mesh_fluid << p_bc << std::flush;
      mfem::socketstream sol_sock3(vishost, visport);
      sol_sock3.precision(8);
      sol_sock3 << "solution\n" << *mesh_fluid << tf_bc << std::flush;
   }

   // SET VALUES ON NO SLIP BOUNDARIES TO ZERO!!!
   mfem::Vector zerovector(2);
   zerovector = 0.0;
   mfem::VectorConstantCoefficient zero(zerovector);
   v_bc.ProjectBdrCoefficient(zero,vdbc_bdr_noslip);

   v0=v_bc;
   p0=p_bc;
   tf0=tf_bc;

   return true;
}


bool NurbsStokesSolver::calc_dirichletbc_fluid_cht(mfem::GridFunction &tf0,mfem::GridFunction &ts0)
{
   mfem::GridFunction tf_bc(tffes); // to calculate our gridfunction on the dirichlet boundary
   // we need grid functions to first compute the controlpoint values on the boundary, so we can project them on to our system
   // means we will build a system that needed to be solved for the desired boundary values

   // create vectors for coefficients and boundary markers

   std::vector<mfem::Array<int>> tfdbc_bdr_marker;
   std::vector<mfem::ConstantCoefficient> tfdbc_bdr_coefficient;

   for (size_t i = 0; i < user_tfdbc_bdr.Size(); i++)
   {
      tfdbc_bdr_marker.push_back(mfem::Array<int>(mesh_fluid->bdr_attributes.Max()));
      tfdbc_bdr_marker[i] = 0;
      tfdbc_bdr_marker[i][user_tfdbc_bdr[i]-1] = 1;

      mfem::ConstantCoefficient cc(user_tfdbc_bdr_values[i]);
      tfdbc_bdr_coefficient.push_back(cc);
      //std::cout << user_tdbc_bdr_values[i] << " temperature values \n";
   }

   InterfaceDirichletCoefficient ifacecoef(beta_t);
   ifacecoef.SetGridFunctionSource(ts0);
   ifacecoef.SetGridFunctionTarget(tf0);

   // TEMPERATURE
   // define rhs with the desired boundary condition values
   mfem::ConstantCoefficient One_bc(1); // coefficient for the kinematic viscosity
   mfem::LinearForm *h_bc(new mfem::LinearForm(tffes)); // define linear form for rhs
   // define bilinear form add the boundary, means the nurbs add the boundary
   mfem::BilinearForm d_bc(tffes); // define the bilinear form results in n x n matrix, we use the velocity finite element space
   h_bc->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(ifacecoef), tfiface_bdr);
   d_bc.AddBoundaryIntegrator(new mfem::MassIntegrator(One_bc),tfiface_bdr);
   for (size_t i = 0; i < tfdbc_bdr_marker.size(); i++)
   {
      h_bc->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(tfdbc_bdr_coefficient[i]),tfdbc_bdr_marker[i]); // define integrator on desired boundary
      d_bc.AddBoundaryIntegrator(new mfem::MassIntegrator(One_bc),tfdbc_bdr_marker[i]); // bilinear form (lambda*u_vector),(v_vector))
   }
   h_bc->Assemble(); // assemble the linear form (vector)   
   d_bc.Assemble(); // assemble the bilinear form (matrix)

   mfem::SparseMatrix D_BC;
   mfem::Vector T_BC, H_BC;
   
   d_bc.SetDiagonalPolicy(mfem::Matrix::DiagonalPolicy::DIAG_ZERO); // important, otherwise a different policy will be used, which results in false building of our matrix
   d_bc.FormLinearSystem(tempf_ess_tdof_list_dummy, tf_bc, *h_bc, D_BC, T_BC, H_BC); // form D_BC
   

   mfem::GMRESSolver bc_solver;
   
   bc_solver.SetAbsTol(atol);
   bc_solver.SetRelTol(rtol);
   bc_solver.SetMaxIter(maxIter);
   bc_solver.SetOperator(D_BC);
   bc_solver.SetKDim((int)maxIter/5);
   bc_solver.SetPrintLevel(3);
   // solve the system
   std::cout << "SOLVE DIRICHLET TEMPERATUREFIELD FLUID\n";
   bc_solver.Mult(H_BC, T_BC);

   // check if solver converged
   if (bc_solver.GetConverged())
   {
      std::cout << "Solver converged in " << bc_solver.GetNumIterations()
                << " iterations with a residual norm of "
                << bc_solver.GetFinalNorm() << ".\n";
   }
   else
   {
      std::cout << "Solver did not converge in " << bc_solver.GetNumIterations()
                << " iterations. Residual norm is " << bc_solver.GetFinalNorm()
                << ".\n";
   }
   //std::cout << "MINRES solver took " << chrono.RealTime() << "s.\n";

   // Save the mesh and the solution
   {
      std::ofstream mesh_ofs(jobname + "_" + std::string(meshfile_fluid));
      mesh_ofs.precision(8);
      mesh_fluid->Print(mesh_ofs);

      std::ofstream tf_bc_ofs(jobname + "_" +"sol_tf_bc.gf");
      tf_bc_ofs.precision(8);
      tf_bc.Save(tf_bc_ofs);
      tf0.Save(tf_bc_ofs);
   }
   if (false)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      mfem::socketstream sol_sock(vishost, visport);
      sol_sock.precision(8);
      sol_sock << "solution\n" << *mesh_fluid << tf_bc << std::flush;
   }

   for (size_t i = 0; i < tf0.Size(); i++)
   {
      //std::cout <<  tf0[i] << " --- " << tf_bc[i] << " \n";
      if (tf_bc[i]!=0)
      {
         tf0[i] = tf_bc[i];
         //std::cout <<  tf0[i] << " --- " << tf_bc[i] << "<--- override \n";
      }
      
   }
   //tf0=tf_bc;

   return true;
}


bool NurbsStokesSolver::calc_dirichletbc_solid(mfem::GridFunction &ts0)
{
   mfem::GridFunction ts_bc(tsfes); // to calculate our gridfunction on the dirichlet boundary

   // we need grid functions to first compute the controlpoint values on the boundary, so we can project them on to our system
   // means we will build a system that needed to be solved for the desired boundary values

   // create vectors for coefficients and boundary markers

   std::vector<mfem::Array<int>> tsdbc_bdr_marker;
   std::vector<mfem::ConstantCoefficient> tsdbc_bdr_coefficient;

   for (size_t i = 0; i < user_tsdbc_bdr.Size(); i++)
   {
      tsdbc_bdr_marker.push_back(mfem::Array<int>(mesh_solid->bdr_attributes.Max()));
      tsdbc_bdr_marker[i] = 0;
      tsdbc_bdr_marker[i][user_tsdbc_bdr[i]-1] = 1;

      mfem::ConstantCoefficient cc(user_tsdbc_bdr_values[i]);
      tsdbc_bdr_coefficient.push_back(cc);
   }

   // TEMPERATURE
   // define rhs with the desired boundary condition values
   mfem::ConstantCoefficient One_bc(1); // coefficient for the kinematic viscosity
   mfem::LinearForm *h_bc(new mfem::LinearForm(tsfes)); // define linear form for rhs
   // define bilinear form add the boundary, means the nurbs add the boundary
   mfem::BilinearForm d_bc(tsfes); // define the bilinear form results in n x n matrix, we use the velocity finite element space
   for (size_t i = 0; i < tsdbc_bdr_marker.size(); i++)
   {
      h_bc->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(tsdbc_bdr_coefficient[i]),tsdbc_bdr_marker[i]); // define integrator on desired boundary
      d_bc.AddBoundaryIntegrator(new mfem::MassIntegrator(One_bc),tsdbc_bdr_marker[i]); // bilinear form (lambda*u_vector),(v_vector))
   }
   h_bc->Assemble(); // assemble the linear form (vector)   
   d_bc.Assemble(); // assemble the bilinear form (matrix)

   mfem::SparseMatrix D_BC;
   mfem::Vector T_BC, H_BC;
   d_bc.SetDiagonalPolicy(mfem::Matrix::DiagonalPolicy::DIAG_ZERO); // important, otherwise a different policy will be used, which results in false building of our matrix
   d_bc.FormLinearSystem(temps_ess_tdof_list_dummy, ts_bc, *h_bc, D_BC, T_BC, H_BC); // form D_BC
   
   mfem::GMRESSolver bc_solver;
   
   bc_solver.SetAbsTol(atol);
   bc_solver.SetRelTol(rtol);
   bc_solver.SetMaxIter(maxIter);
   bc_solver.SetOperator(D_BC);
   bc_solver.SetKDim((int)maxIter/5);
   bc_solver.SetPrintLevel(3);

   // solve the system
   bc_solver.Mult(H_BC, T_BC);

   // check if solver converged
   if (bc_solver.GetConverged())
   {
      std::cout << "Solver converged in " << bc_solver.GetNumIterations()
                << " iterations with a residual norm of "
                << bc_solver.GetFinalNorm() << ".\n";
   }
   else
   {
      std::cout << "Solver did not converge in " << bc_solver.GetNumIterations()
                << " iterations. Residual norm is " << bc_solver.GetFinalNorm()
                << ".\n";
   }
   //std::cout << "MINRES solver took " << chrono.RealTime() << "s.\n";

   // Save the mesh and the solution
   {
      std::ofstream mesh_ofs(jobname + "_" + std::string(meshfile_solid));
      mesh_ofs.precision(8);
      mesh_solid->Print(mesh_ofs);

      std::ofstream ts_bc_ofs(jobname + "_" + "sol_ts_bc.gf");
      ts_bc_ofs.precision(8);
      ts_bc.Save(ts_bc_ofs);
      ts0.Save(ts_bc_ofs);
   }
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      mfem::socketstream sol_sock(vishost, visport);
      sol_sock.precision(8);
      sol_sock << "solution\n" << *mesh_solid << ts_bc << std::flush;
   }

   ts0=ts_bc;

   return true;
}

bool NurbsStokesSolver::calc_flowsystem_strongbc(mfem::GridFunction &v0,mfem::GridFunction &p0, mfem::GridFunction &tf0, mfem::GridFunction &v, mfem::GridFunction &p, mfem::GridFunction &tf, mfem::Coefficient &kin_vis)
{  
   //VALIDATION COUETTE FLOW
   //p0 = 100;
   //pres_ess_tdof_list = 1;

   // to set our boundary conditions we first ned to define our grid functions, so that we have something to project onto
   // we need Blockoperators to define the equation system
   // Implement Blockoperators
   
   mfem::Array<int> block_offsets(3);  // number of variables + 1
   block_offsets[0] = 0;
   block_offsets[1] = vfes->GetVSize();
   block_offsets[2] = pfes->GetVSize();
   block_offsets.PartialSum();

   /*
   std::cout << "***********************************************************\n";
   std::cout << "dim(v) = " << block_offsets[1] - block_offsets[0] << "\n";
   std::cout << "dim(p) = " << block_offsets[2] - block_offsets[1] << "\n";
   std::cout << "dim(v+p) = " << block_offsets.Last() << "\n";
   std::cout << "dim(t) = " << tfes->GetVSize() << "\n";
   std::cout << "***********************************************************\n" << std::endl;
   */

   mfem::BlockVector x_flow(block_offsets), rhs_flow(block_offsets); // blockvector for gridfunctions and our rhs
      /*
         x_flow = [ v ]      rhs_flow =   [ f ]
                  [ p ]                   [ g ]
         
         NO BLOCKVECTOR NEEDED ->   x_temperature = [ t ]  rhs_temperature = [ h ]
         
      */

   x_flow = 0.0;
   rhs_flow = 0.0;
   //x_temperature = 0.0;
   //rhs_temperature = 0.0;
   
   // make reference to block vector
   v.MakeRef(vfes, x_flow.GetBlock(0), 0);
   p.MakeRef(pfes, x_flow.GetBlock(1), 0);
   v = v0;
   p = p0;

   // Setup bilinear and linear forms

   // rhs of momentum equation
   // we don't have a source term in our equations, so we just make an zero source term
   mfem::Vector vzero(sdim);
   vzero = 0.;
   mfem::VectorConstantCoefficient vcczero(vzero); // coefficient for source term
      
   mfem::LinearForm *f(new mfem::LinearForm); // define linear form for rhs
   f->Update(vfes, rhs_flow.GetBlock(0),0);  // link to block vector and use the velocity finite element space
   f->AddDomainIntegrator(new mfem::VectorDomainLFIntegrator(vcczero));
   //f->AddBoundaryIntegrator(new mfem::VectorBoundaryLFIntegrator(vccpressure), vdbc_bdr);
   f->Assemble(); // assemble the linear form (vector)
   
   // rhs for continuity equation
   mfem::ConstantCoefficient zero(0.0); // zero source term
   mfem::LinearForm *g(new mfem::LinearForm); // define linear form for rhs
   g->Update(pfes, rhs_flow.GetBlock(1), 0); // link to block vector and use the pressure finite element space
   g->AddDomainIntegrator(new mfem::DomainLFIntegrator(zero)); // define integrator for source term -> zero in our case
   g->Assemble(); // assemble the linear form (vector)
   

   // Momentum equation
   // diffusion term
   mfem::BilinearForm a(vfes); // define the bilinear form results in n x n matrix, we use the velocity finite element space
   //mfem::ConstantCoefficient kin_vis(kin_viscosity); // coefficient for the kinematic viscosity
   
   // chose viscosity model
   CarreauModelCoefficient* temp_kin_vis_carreau;
   CarreauWLFModelCoefficient* temp_kin_vis_carreauwlf;
   PowerLawModelCoefficient* temp_kin_vis_powerlaw;
   
   if (temp_kin_vis_carreau = dynamic_cast<CarreauModelCoefficient*>(&kin_vis))
   {
      CarreauModelCoefficient kin_vis_carreau = *temp_kin_vis_carreau;
      kin_vis_carreau.SetVelocity(v0);
      a.AddDomainIntegrator(new mfem::VectorDiffusionIntegrator(kin_vis_carreau, sdim)); // bilinear form (lambda*nabla(u_vector),nabla(v_vector))
      a.Assemble();
   }else if (temp_kin_vis_carreauwlf = dynamic_cast<CarreauWLFModelCoefficient*>(&kin_vis))
   {
      CarreauWLFModelCoefficient kin_vis_carreauwlf = *temp_kin_vis_carreauwlf;
      kin_vis_carreauwlf.SetVelocity(v0);
      kin_vis_carreauwlf.SetTemperature(tf0);
      a.AddDomainIntegrator(new mfem::VectorDiffusionIntegrator(kin_vis_carreauwlf, sdim)); // bilinear form (lambda*nabla(u_vector),nabla(v_vector))
      a.Assemble();
   }
   else if (temp_kin_vis_powerlaw = dynamic_cast<PowerLawModelCoefficient*>(&kin_vis))
   {
      PowerLawModelCoefficient kin_vis_powerlaw = *temp_kin_vis_powerlaw;
      kin_vis_powerlaw.SetVelocity(v0);
      kin_vis_powerlaw.SetTemperature(tf0);
      a.AddDomainIntegrator(new mfem::VectorDiffusionIntegrator(kin_vis_powerlaw, sdim)); // bilinear form (lambda*nabla(u_vector),nabla(v_vector))
      a.Assemble();
   }
   else{
      a.AddDomainIntegrator(new mfem::VectorDiffusionIntegrator(kin_vis, sdim)); // bilinear form (lambda*nabla(u_vector),nabla(v_vector))
      a.Assemble();
   }

   //a.AddDomainIntegrator(new mfem::VectorDiffusionIntegrator(kin_vis_used, sdim)); // bilinear form (lambda*nabla(u_vector),nabla(v_vector))
   //a.Assemble(); // assemble the bilinear form (matrix)
   //a.Finalize(); not needed, will be called on form linear system

   // grad pressure term
   // define the mixed bilinear form results in n x m matrix, we use the velocity finite element space as test space and the pressure space as trial space
   mfem::MixedBilinearForm b(pfes,vfes); // (trial,test)
   // no clue right now about the right sign
   mfem::ConstantCoefficient coef_b(1.0/density); // -1 because of the sign in the equation
   b.AddDomainIntegrator(new mfem::GradientIntegrator(coef_b)); // mixed bilinear form (lambda*nabla(u),v_vector)
   //b.AddDomainIntegrator(new mfem::MixedScalarWeakGradientIntegrator(minusOne)); // mixed bilinear form (lambda*nabla(u),v_vector)
   b.Assemble(); // assemble the mixed bilinear form (matrix)
   //b.Finalize(); not needed, will be called on form linear system

   // continuity term
   // define the mixed bilinear form results in n x m matrix, we use the pressure finite element space as test space and the velocity space as trial space
   mfem::MixedBilinearForm c(vfes,pfes); // (trial,test)
   mfem::ConstantCoefficient coef_c(1.0); // +1 because of the sign in the equation
   c.AddDomainIntegrator(new mfem::VectorDivergenceIntegrator(coef_c)); // mixed bilinear form (lambda*nabla . u_vector, v)
   c.Assemble(); // assemble the mixed bilinear form (matrix)
   //c.Finalize(); not needed, will be called on form linear system


   // we need some SparseMatrix and Vector to form our linear system
   mfem::SparseMatrix A,B,C;
   mfem::Vector V, F;
   mfem::Vector P, G;
   a.SetDiagonalPolicy(mfem::Matrix::DiagonalPolicy::DIAG_ZERO); // important, otherwise a different policy will be used, which results in false building of our matrix
   a.FormLinearSystem(vel_ess_tdof_list, v0, *f, A, V, F); // form A   
   b.FormRectangularLinearSystem(pres_ess_tdof_list, vel_ess_tdof_list, p0, *f, B, P, F); // form B
   c.FormRectangularLinearSystem(vel_ess_tdof_list, pres_ess_tdof_list, v0, *g, C, V, G); // form C

   mfem::BlockOperator stokesOp(block_offsets); // Block operator to build our System for the solver

   stokesOp.SetBlock(0,0,&A);
   stokesOp.SetBlock(0,1,&B);
   stokesOp.SetBlock(1,0,&C);
   
   //mfem::TransposeOperator *Bt = NULL;
   //Bt = new mfem::TransposeOperator(&B);
   //stokesOp.SetBlock(1,0,Bt);

   mfem::StopWatch chrono; // stop watch to calc solve time
   chrono.Clear();
   chrono.Start();
   
   // SOLVER
   // setup solver

   // setup minres solver, should be enough for our linear system
   // without preconditioning
   //mfem::MINRESSolver solver; 
   mfem::GMRESSolver solver;
   //solver.iterative_mode = false;
   solver.SetAbsTol(atol);
   solver.SetRelTol(rtol);
   solver.SetMaxIter(maxIter);
   solver.SetOperator(stokesOp);
   solver.SetKDim((int)maxIter/5);
   solver.SetPrintLevel(3);
   
   //std::cout << rhs_flow.BlockSize(0)  <<"\n";
   //std::cout << rhs_flow.BlockSize(1)  <<"\n";
   //std::cout << x_flow.BlockSize(0)  <<"\n";
   //std::cout << x_flow.BlockSize(1)  <<"\n";
   
   // solve the system
   std::cout << "SOLVE FLOWFIELD \n";
   // initial guess
   V = v0;
   P = p0;
   solver.Mult(rhs_flow, x_flow);
   chrono.Stop();

   // check if solver converged
   if (solver.GetConverged())
   {
      std::cout << "GMRESSolver converged in " << solver.GetNumIterations()
                << " iterations with a residual norm of "
                << solver.GetFinalNorm() << ".\n";
   }
   else
   {
      std::cout << "GMRESSolver did not converge in " << solver.GetNumIterations()
                << " iterations. Residual norm is " << solver.GetFinalNorm()
                << ".\n";
   }
   std::cout << "GMRESSolver solver took " << chrono.RealTime() << "s.\n";

   // Save the mesh and the solution
   {
      std::ofstream mesh_ofs(jobname + "_" + std::string(meshfile_fluid));
      mesh_ofs.precision(8);
      mesh_fluid->Print(mesh_ofs);

      std::ofstream v_ofs(jobname + "_" +"sol_v.gf");
      v_ofs.precision(8);
      v.Save(v_ofs);

      std::ofstream p_ofs(jobname + "_" +"sol_p.gf");
      p_ofs.precision(8);
      p.Save(p_ofs);
   }
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      mfem::socketstream sol_sock1(vishost, visport);
      sol_sock1.precision(8);
      sol_sock1 << "solution\n" << *mesh_fluid << v << std::flush;
      mfem::socketstream sol_sock2(vishost, visport);
      sol_sock2.precision(8);
      sol_sock2 << "solution\n" << *mesh_fluid << p << std::flush;
   }

   return true;
}

bool NurbsStokesSolver::calc_temperaturesystem_strongbc_fluid(mfem::GridFunction &v0, mfem::GridFunction &tf0,mfem::GridFunction &ts0, mfem::GridFunction &v, mfem::GridFunction &tf,mfem::GridFunction &ts)
{
   // Setup bilinear and linear forms
   tf = tf0;
   // rhs for advection diffusion heat transfer
   mfem::ConstantCoefficient zero(0.0); // zero source term
   mfem::LinearForm *h(new mfem::LinearForm(tffes)); // define linear form for rhs
   h->AddDomainIntegrator(new mfem::DomainLFIntegrator(zero)); // define integrator for source term -> zero in our case
   h->Assemble(); // assemble the linear form (vector)

   // advection diffusion heat transfer
   mfem::BilinearForm d(tffes); // define the bilinear form results in n x n matrix, we use the temperature finite element space
   mfem::ConstantCoefficient temp_dcoeff(temp_diffusion_const_fluid); // coefficient for the temp_diffusion_const
   mfem::VectorGridFunctionCoefficient v_coef;
   v_coef.SetGridFunction(&v0);
   //mfem::VectorFunctionCoefficient v_coef(sdim, lambda_velocityfield);
   d.AddDomainIntegrator(new mfem::DiffusionIntegrator(temp_dcoeff)); // bilinear form (lambda*nabla(u),nabla(v))
   d.AddDomainIntegrator(new mfem::ConvectionIntegrator(v_coef,1)); // 
   //d.AddDomainIntegrator(new mfem::MixedDirectionalDerivativeIntegrator(v_coef)); // 
   d.Assemble(); // assemble the bilinear form (matrix)
   //a.Finalize(); not needed, will be called on form linear system

   // we need some SparseMatrix and Vector to form our linear system
   mfem::SparseMatrix D;
   mfem::Vector T, H;
   d.SetDiagonalPolicy(mfem::Matrix::DiagonalPolicy::DIAG_ZERO); // important, otherwise a different policy will be used, which results in false building of our matrix
   d.FormLinearSystem(tempf_ess_tdof_list, tf, *h, D, T, H); // form D
 
   //D.PrintInfo(std::cout);
   //D.PrintInfo(std::cout);
   //D.PrintMatlab(std::cout);
   
   mfem::StopWatch chrono; // stop watch to calc solve time
   chrono.Clear();
   chrono.Start();
   
   // SOLVER
   // setup solver
   //int maxIter=1000; // maximal number of iterations
   //double rtol(1.e-10); // convergence criteria
   //double atol(1.e-10); // convergence criteria

   // setup minres solver, should be enough for our linear system
   // without preconditioning
   mfem::GMRESSolver solver;
   //solver.iterative_mode = false;
   
   //mfem::MINRESSolver solver;
   solver.SetAbsTol(atol);
   solver.SetRelTol(rtol);
   solver.SetMaxIter(maxIter);
   solver.SetOperator(D);
   solver.SetKDim((int)maxIter/5);
   solver.SetPrintLevel(3);

   std::cout << "SOLVE TEMPERATUREFIELD FLUID \n";   
   // solve the system
   // initial guess
   T = tf0;
   //solver.Mult(H, tf);
   solver.Mult(H, T);
   chrono.Stop();
   //std::cout << v0;
 
   // check if solver converged
   if (solver.GetConverged())
   {
      std::cout << "GMRESSolver converged in " << solver.GetNumIterations()
                << " iterations with a residual norm of "
                << solver.GetFinalNorm() << " .\n";
   }
   else
   {
      std::cout << "GMRESSolver did not converge in " << solver.GetNumIterations()
                << " iterations. Residual norm is " << solver.GetFinalNorm()
                << ".\n";
   }
   std::cout << "GMRESSolver solver took " << chrono.RealTime() << "s.\n";

   // Save the mesh and the solution
   {
      std::ofstream mesh_ofs(jobname + "_" + std::string(meshfile_fluid));
      mesh_ofs.precision(8);
      mesh_fluid->Print(mesh_ofs);

      std::ofstream tf_ofs(jobname + "_" +"sol_tf.gf");
      tf_ofs.precision(8);
      tf.Save(tf_ofs);
   }
   
   return true;
}

bool NurbsStokesSolver::calc_temperaturesystem_strongbc_solid(mfem::GridFunction &ts0,mfem::GridFunction &tf0, mfem::GridFunction &ts,mfem::GridFunction &tf, std::vector<double> &flux)
{
   // Setup bilinear and linear forms
   ts = ts0;

   // rhs for advection diffusion heat transfer
   mfem::ConstantCoefficient zero(0.0); // zero source term
   mfem::LinearForm *h(new mfem::LinearForm(tsfes)); // define linear form for rhs
   h->AddDomainIntegrator(new mfem::DomainLFIntegrator(zero)); // define integrator for source term -> zero in our case
   
   InterfaceFluxCoefficient ifacecoef(temp_diffusion_const_fluid,temp_diffusion_const_solid,beta_q,flux);
   ifacecoef.SetGridFunctionSource(tf);
   ifacecoef.SetGridFunctionTarget(ts0);
   h->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(ifacecoef), tsiface_bdr);
   
   h->Assemble(); // assemble the linear form (vector)

   // advection diffusion heat transfer
   mfem::BilinearForm d(tsfes); // define the bilinear form results in n x n matrix, we use the temperature finite element space
   mfem::ConstantCoefficient temp_dcoeff(temp_diffusion_const_solid); // coefficient for the temp_diffusion_const
   d.AddDomainIntegrator(new mfem::DiffusionIntegrator(temp_dcoeff)); // bilinear form (lambda*nabla(u),nabla(v))
   d.Assemble(); // assemble the bilinear form (matrix)
   //a.Finalize(); not needed, will be called on form linear system

   // we need some SparseMatrix and Vector to form our linear system
   mfem::SparseMatrix D;
   mfem::Vector T, H;
   d.SetDiagonalPolicy(mfem::Matrix::DiagonalPolicy::DIAG_ZERO); // important, otherwise a different policy will be used, which results in false building of our matrix
   d.FormLinearSystem(temps_ess_tdof_list, ts, *h, D, T, H); // form D
    
   mfem::StopWatch chrono; // stop watch to calc solve time
   chrono.Clear();
   chrono.Start();
   
   // SOLVER
   mfem::GMRESSolver solver;
   //solver.iterative_mode = false;
   
   solver.SetAbsTol(atol);
   solver.SetRelTol(rtol);
   solver.SetMaxIter(maxIter);
   solver.SetOperator(D);
   solver.SetKDim((int)maxIter/5);
   solver.SetPrintLevel(3);

   std::cout << "SOLVE TEMPERATUREFIELD SOLID \n";   
   // solve the system
   // initial guess
   T = ts0;
   //solver.Mult(H, ts);
   solver.Mult(H, T);
   chrono.Stop();
   //std::cout << v0;
 
   // check if solver converged
   if (solver.GetConverged())
   {
      std::cout << "GMRESSolver converged in " << solver.GetNumIterations()
                << " iterations with a residual norm of "
                << solver.GetFinalNorm() << " .\n";
   }
   else
   {
      std::cout << "GMRESSolver did not converge in " << solver.GetNumIterations()
                << " iterations. Residual norm is " << solver.GetFinalNorm()
                << ".\n";
   }
   std::cout << "GMRESSolver solver took " << chrono.RealTime() << "s.\n";

   // Save the mesh and the solution
   {
      std::ofstream mesh_ofs(jobname + "_" +std::string(meshfile_solid));
      mesh_ofs.precision(8);
      mesh_solid->Print(mesh_ofs);

      std::ofstream ts_ofs(jobname + "_" +"sol_ts.gf");
      ts_ofs.precision(8);
      ts.Save(ts_ofs);
   }
   if (false)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      mfem::socketstream sol_sock(vishost, visport);
      sol_sock.precision(8);
      sol_sock << "solution\n" << *mesh_solid << ts << std::flush;
   }

   return true;
}

bool NurbsStokesSolver::calc_temperaturesystem_strongbc_solid_init(mfem::GridFunction &ts0,mfem::GridFunction &tf0, mfem::GridFunction &ts,mfem::GridFunction &tf)
{
   // Setup bilinear and linear forms
   ts = ts0;

   // rhs for advection diffusion heat transfer
   mfem::ConstantCoefficient zero(0.0); // zero source term
   mfem::LinearForm *h(new mfem::LinearForm(tsfes)); // define linear form for rhs
   h->AddDomainIntegrator(new mfem::DomainLFIntegrator(zero)); // define integrator for source term -> zero in our case
   /*
   InterfaceFluxCoefficient ifacecoef(temp_diffusion_const_fluid);
   ifacecoef.SetGridFunctionSource(tf);
   ifacecoef.SetGridFunctionTarget(ts0);
   h->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(ifacecoef), tsiface_bdr);
   */
   h->Assemble(); // assemble the linear form (vector)

   // advection diffusion heat transfer
   mfem::BilinearForm d(tsfes); // define the bilinear form results in n x n matrix, we use the temperature finite element space
   mfem::ConstantCoefficient temp_dcoeff(temp_diffusion_const_solid); // coefficient for the temp_diffusion_const
   d.AddDomainIntegrator(new mfem::DiffusionIntegrator(temp_dcoeff)); // bilinear form (lambda*nabla(u),nabla(v))
   d.Assemble(); // assemble the bilinear form (matrix)
   //a.Finalize(); not needed, will be called on form linear system

   // we need some SparseMatrix and Vector to form our linear system
   mfem::SparseMatrix D;
   mfem::Vector T, H;
   d.SetDiagonalPolicy(mfem::Matrix::DiagonalPolicy::DIAG_ZERO); // important, otherwise a different policy will be used, which results in false building of our matrix
   d.FormLinearSystem(temps_ess_tdof_list, ts, *h, D, T, H); // form D
    
   mfem::StopWatch chrono; // stop watch to calc solve time
   chrono.Clear();
   chrono.Start();
   
   // SOLVER
   mfem::GMRESSolver solver;
   //solver.iterative_mode = false;
   
   solver.SetAbsTol(atol);
   solver.SetRelTol(rtol);
   solver.SetMaxIter(maxIter);
   solver.SetOperator(D);
   solver.SetKDim((int)maxIter/5);
   solver.SetPrintLevel(3);

   std::cout << "SOLVE TEMPERATUREFIELD SOLID \n";   
   // solve the system
   solver.Mult(H, ts);
   chrono.Stop();
   //std::cout << v0;
 
   // check if solver converged
   if (solver.GetConverged())
   {
      std::cout << "GMRESSolver converged in " << solver.GetNumIterations()
                << " iterations with a residual norm of "
                << solver.GetFinalNorm() << " .\n";
   }
   else
   {
      std::cout << "GMRESSolver did not converge in " << solver.GetNumIterations()
                << " iterations. Residual norm is " << solver.GetFinalNorm()
                << ".\n";
   }
   std::cout << "GMRESSolver solver took " << chrono.RealTime() << "s.\n";

   // Save the mesh and the solution
   {
      std::ofstream mesh_ofs(jobname + "_" +std::string(meshfile_solid));
      mesh_ofs.precision(8);
      mesh_solid->Print(mesh_ofs);

      std::ofstream ts_ofs(jobname + "_" +"sol_ts.gf");
      ts_ofs.precision(8);
      ts.Save(ts_ofs);
   }
   if (false)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      mfem::socketstream sol_sock(vishost, visport);
      sol_sock.precision(8);
      sol_sock << "solution\n" << *mesh_solid << ts << std::flush;
   }

   return true;
}


bool NurbsStokesSolver::solve_flow(mfem::GridFunction &v0,mfem::GridFunction &p0, mfem::GridFunction &tf0, mfem::GridFunction &v, mfem::GridFunction &p, mfem::GridFunction &tf,mfem::Coefficient &kin_vis)
{  
   if (bcstrong)
   {
      calc_flowsystem_strongbc(v0, p0, tf0, v, p, tf, kin_vis);
   } else if (bcweak)
   {
      /* code */
   }
      
   return true;
}

bool NurbsStokesSolver::solve_temperature(mfem::GridFunction &v0, mfem::GridFunction &tf0, mfem::GridFunction &ts0, mfem::GridFunction &v, mfem::GridFunction &tf, mfem::GridFunction &ts, std::vector<double> &flux)
{

   if (bcstrong)
   {  
      calc_dirichletbc_fluid_cht(tf0,ts0);
      //std::cout << "\n\n ## dirichlet geht \n\n";
      calc_temperaturesystem_strongbc_fluid(v0, tf0, ts0, v, tf,ts);
      //std::cout << "\n\n ##  temp f geht \n\n";
      calc_temperaturesystem_strongbc_solid(ts0, tf0, ts, tf, flux);
      //std::cout << "\n\n ##  temp s geht \n\n";
   } else if (bcweak)
   {
      /* code */
   }
   return true;
}

double NurbsStokesSolver::error_norm_abs(mfem::GridFunction &x0, mfem::GridFunction &x1)
{
   double norm = 1.;
   mfem::GridFunction Res(x0);

   for (size_t i = 0; i < x0.Size(); i++)
   {  
      Res[i] = x1[i] - x0[i];
   }

   norm = Res.Norml2();
   //std::cout << "\n\n ##  error_norm_abs \n\n";
   return norm;
}

double NurbsStokesSolver::error_norm_rel(mfem::GridFunction &x0, mfem::GridFunction &x1)
{
   double norm = 1.;  
   mfem::GridFunction Res(x0);

   for (size_t i = 0; i < x0.Size(); i++)
   {  
      Res[i] = x1[i] - x0[i];
   }
   
   if(x1.Norml2()!=0){
      norm = Res.Norml2()/x1.Norml2();
   }
   else
   {
      norm = Res.Norml2();
   }
   
   //std::cout << "\n\n ##  error_norm_rel \n\n";
   return norm;
}

double NurbsStokesSolver::error_norm_abs_vector(std::vector<double> &x0, std::vector<double> &x1)
{
   
   double norm = 0.;
   double norm_x1 = 0.;
   std::vector<double> Res(x0.size());

   for (size_t i = 0; i < x0.size(); i++)
   {  
      Res[i] = x1[i] - x0[i];
      norm += Res[i]*Res[i];
   }
   
   norm = std::sqrt(norm);
   //std::cout << "\n\n ##  error_norm_abs_vector \n\n";
   return norm;
}

double NurbsStokesSolver::error_norm_rel_vector(std::vector<double> &x0, std::vector<double> &x1)
{
   double norm = 0.;
   double norm_x1 = 0.;
   std::vector<double> Res(x0.size());

   for (size_t i = 0; i < x0.size(); i++)
   {  
      Res[i] = x1[i] - x0[i];
      norm += Res[i]*Res[i];
      norm_x1 += x1[i] * x1[i];
   }
   if(std::sqrt(norm_x1!=0)){
      norm = std::sqrt(norm)/std::sqrt(norm_x1);
   }
   else
   {
      norm = std::sqrt(norm);
   }
   //std::cout << "\n\n ##  error_norm_rel_vector \n\n";  
   return norm;
}
;