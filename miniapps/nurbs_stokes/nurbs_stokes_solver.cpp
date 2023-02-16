#include "nurbs_stokes_solver.hpp"
#include <fstream>  // fstream for input and output in textfiles
#include <iostream> // for output to the console

#include <algorithm> // to get access to math functions
#include <cmath>     // to get access to math functions

NurbsStokesSolver::NurbsStokesSolver()
{
   v_max = 280;               // max velocity for our boundary on the inlet
   p_val = 100;           // value for pressure boundary
   kin_viscosity = 200;    // value for kinematic visosity
   temp_1 = 0;               // value for temperature
   temp_2 = 50;               // value for temperature
   temp_diffusion_const = 0.1; // value for temperature diffusion constant coefficient
   meshfile = "../../../MA/data/quad_nurbs.mesh";
   ref_levels = 0;
   order = mfem::Array<int>(4);
   order[0] = 0;                    // mesh order 
   order[1] = 0;                    // order elevation velocity
   order[2] = 0;                    // order elevation pressure
   order[3] = 0;                    // order elevation temperature
}

NurbsStokesSolver::~NurbsStokesSolver()
{}

bool NurbsStokesSolver::init()
{
  if (is_initialized)
  {
      return false; // already initialized
  }else{
      // Read the mesh from the given mesh file. 
      mesh = new mfem::Mesh(meshfile, 1, 1);
      sdim = mesh->SpaceDimension(); // get dimension in space from the mesh 1D, 2D, 3D
      order[0] = mesh->NURBSext->GetOrder();
      // mesh refinement
      if (ref_levels < 0) // if ref_level parameter is smaller then 0, set a standard refinement
      {
         ref_levels =
            (int)floor(log(500./mesh->GetNE())/log(2.)/sdim);
      }
      for (int l = 0; l < ref_levels; l++) // loop for refining of the mesh
      {
         mesh->UniformRefinement();
      }
      mesh->PrintInfo();
      std::cout << "Using Nurbs mesh." << std::endl;
      std::cout << "mesh Space Dimension " << sdim << std::endl;
      std::cout << "given Mesh Order " << order[0] << std::endl;
      std::cout << "given Order elevation velocity " << order[1] << std::endl;
      std::cout << "given Order elevation pressure" << order[2] << std::endl;
      std::cout << "given Order elevation temperature" << order[3] << std::endl;
      std::cout << " bdr attributes " << mesh->bdr_attributes.Max() << std::endl; // output the number of boundary attributes
      // fe collections
      vfec = new mfem::NURBSFECollection(order[0] + order[1]); // Pm+n
      pfec = new mfem::NURBSFECollection(order[0] + order[2]); // Pm+n
      tfec = new mfem::NURBSFECollection(order[0] + order[3]); // Pm+n
      std::cout << "velocity finite element collection Order " << vfec->GetOrder() << std::endl;
      std::cout << "pressure finite element collection Order " << pfec->GetOrder() << std::endl;
      std::cout << "temperature finite element collection Order " << tfec->GetOrder() << std::endl;
      // nurbsextensions
      vNURBSext = new mfem::NURBSExtension(mesh->NURBSext, order[0]+order[1]);
      pNURBSext = new mfem::NURBSExtension(mesh->NURBSext, order[0]+order[2]);
      tNURBSext = new mfem::NURBSExtension(mesh->NURBSext, order[0]+order[3]);

      // declaration for our finite element spaces
      vfes = new mfem::FiniteElementSpace(mesh, vNURBSext, vfec, sdim); // velocity finite element space, with dimension sdim
      pfes = new mfem::FiniteElementSpace(mesh, pNURBSext, pfec);  // pressure finite element space, with dimension 1 (scalar field)
      tfes = new mfem::FiniteElementSpace(mesh, pNURBSext, tfec);  // temperature finite element space, with dimension 1 (scalar field)
      std::cout << "vfes velocity finite element space Order " << vfes->GetMaxElementOrder() << std::endl;
      std::cout << "pfes pressure finite element space Order " << pfes->GetMaxElementOrder() << std::endl;
      std::cout << "tfes temperature finite element space Order " << tfes->GetMaxElementOrder() << std::endl;
      std::cout << "Number of finite element unknowns vfes: "
      << vfes->GetTrueVSize() << std::endl;
      std::cout << "Number of finite element unknowns pfes: "
      << pfes->GetTrueVSize() << std::endl;
      std::cout << "Number of finite element unknowns tfes: "
      << tfes->GetTrueVSize() << std::endl;

      // init the dirichlet bc
      this->init_dirichletbc();

      // to set our boundary conditions we first ned to define our grid functions, so that we have something to project onto
      // we need Blockoperators to define the equation system
      // Implement Blockoperators
      block_offsets = mfem::Array<int>(3); // number of variables + 1
      block_offsets[0] = 0;
      block_offsets[1] = vfes->GetVSize();
      block_offsets[2] = pfes->GetVSize();
      block_offsets.PartialSum();

      std::cout << "***********************************************************\n";
      std::cout << "dim(v) = " << block_offsets[1] - block_offsets[0] << "\n";
      std::cout << "dim(p) = " << block_offsets[2] - block_offsets[1] << "\n";
      std::cout << "dim(v+p) = " << block_offsets.Last() << "\n";
      std::cout << "dim(t) = " << tfes->GetVSize() << "\n";
      std::cout << "***********************************************************\n" << std::endl;

      x_flow = mfem::BlockVector(block_offsets);
      rhs_flow = mfem::BlockVector(block_offsets); // blockvector for gridfunctions and our rhs

      x_flow = 0.0;
      rhs_flow = 0.0;
      //x_temperature = 0.0;
      //rhs_temperature = 0.0;

      // make reference to block vector
      v.MakeRef(vfes, x_flow.GetBlock(0), 0);
      p.MakeRef(pfes, x_flow.GetBlock(1), 0);

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
   vdbc_bdr = mfem::Array<int>(mesh->bdr_attributes.Max()); // contains the whole boundary markers
   vdbc_bdr_noslip = mfem::Array<int>(mesh->bdr_attributes.Max()); // to select only the boundary markers for the no slip walls
   vdbc_bdr_inlet = mfem::Array<int>(mesh->bdr_attributes.Max());   // to select only the boundary markers for the inlet
   pdbc_bdr = mfem::Array<int>(mesh->bdr_attributes.Max()); // contains the whole boundary markers
   tdbc_bdr = mfem::Array<int>(mesh->bdr_attributes.Max()); // contains the whole boundary markers
   tdbc_bdr_inlet = mfem::Array<int>(mesh->bdr_attributes.Max());   // to select only the boundary markers for the inlet
   tdbc_bdr_walls= mfem::Array<int> (mesh->bdr_attributes.Max());   // to select only the boundary markers for the walls
   vdummy_bdr = mfem::Array<int>(mesh->bdr_attributes.Max()); // contains the whole boundary markers
   pdummy_bdr = mfem::Array<int>(mesh->bdr_attributes.Max()); // contains the whole boundary markers
   tdummy_bdr = mfem::Array<int>(mesh->bdr_attributes.Max()); // contains the whole boundary markers
   
   // for poiseuille flow 
   vdbc_bdr = 0;
   // the boundary for the noslip condition
   vdbc_bdr_noslip = 0;
   // the boundary at the inlet
   vdbc_bdr_inlet = 0;
   
   // pressure
   pdbc_bdr = 0;

   // temperature
   tdbc_bdr = 0;
   tdbc_bdr_inlet = 0;
   tdbc_bdr_walls = 0;

   // mark dummy
   vdummy_bdr = 0;
   pdummy_bdr = 0;
   tdummy_bdr = 0;


   // get the true dofs of the boundaries
   vfes->GetEssentialTrueDofs(vdbc_bdr, vel_ess_tdof_list);
   pfes->GetEssentialTrueDofs(pdbc_bdr, pres_ess_tdof_list);   
   tfes->GetEssentialTrueDofs(tdbc_bdr, temp_ess_tdof_list);   

   vfes->GetEssentialTrueDofs(vdummy_bdr, vel_ess_tdof_list_dummy);
   pfes->GetEssentialTrueDofs(pdummy_bdr, pres_ess_tdof_list_dummy);   
   tfes->GetEssentialTrueDofs(tdummy_bdr, temp_ess_tdof_list_dummy);

   return true;
}

bool NurbsStokesSolver::set_meshfile(std::string meshfilepath)
{
   meshfile = meshfilepath.c_str();
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

bool NurbsStokesSolver::set_order_elevation_temperature(int elevationorder)
{
   order[3] = elevationorder;
   return true;
}

bool NurbsStokesSolver::set_parameters(std::vector<double> parameters)
{
   v_max = parameters[0];               // max velocity for our boundary on the inlet
   p_val = parameters[1];           // value for pressure boundary
   kin_viscosity = parameters[2];    // value for kinematic visosity
   temp_1 = parameters[3];               // value for temperature
   temp_2 = parameters[4];               // value for temperature
   temp_diffusion_const = parameters[5]; // value for temperature diffusion constant coefficient
   return true;
}

bool NurbsStokesSolver::set_dirichletbc_velocity_noslip(std::vector<int> boundary_marker)
{
   for (size_t i = 0; i < boundary_marker.size(); i++)
   {  
      if (boundary_marker[i]==1)
      {  
         vdbc_bdr[i] = boundary_marker[i];
         vdbc_bdr_noslip[i] = boundary_marker[i];
      }
   }
   return true;
}

bool NurbsStokesSolver::set_dirichletbc_velocity_inlet(std::vector<int> boundary_marker)
{
   for (size_t i = 0; i < boundary_marker.size(); i++)
   {  
      if (boundary_marker[i]==1)
      {  
         vdbc_bdr[i] = boundary_marker[i];
         vdbc_bdr_inlet[i] = boundary_marker[i];
      }
   }
   return true;
}

bool NurbsStokesSolver::set_dirichletbc_pressure(std::vector<int> boundary_marker)
{
   for (size_t i = 0; i < boundary_marker.size(); i++)
   {  
      if (boundary_marker[i]==1)
      {  
         pdbc_bdr[i] = boundary_marker[i];
      }
   }
   return true;
}
   
bool NurbsStokesSolver::set_dirichletbc_temperature_inlet(std::vector<int> boundary_marker)
{
   for (size_t i = 0; i < boundary_marker.size(); i++)
   {  
      if (boundary_marker[i]==1)
      {  
         tdbc_bdr[i] = boundary_marker[i];
         tdbc_bdr_inlet[i] = boundary_marker[i];
      }
   }
   return true;
}

bool NurbsStokesSolver::set_dirichletbc_temperature_walls(std::vector<int> boundary_marker)
{
   for (size_t i = 0; i < boundary_marker.size(); i++)
   {  
      if (boundary_marker[i]==1)
      {  
         tdbc_bdr[i] = boundary_marker[i];
         tdbc_bdr_walls[i] = boundary_marker[i];
      }
   }
   return true;
}

bool NurbsStokesSolver::calc_dirichletbc()
{
   auto lambda_inlet = [this](const mfem::Vector &QuadraturPointPosition, mfem::Vector &VelocityValue) -> void
   {
      double h=1;
      VelocityValue[0] = 4*(h-QuadraturPointPosition(1))*QuadraturPointPosition(1)/(h*h)*v_max;
      
      //VelocityValue[0] = 0;
      //VelocityValue[1] = 4*(h-QuadraturPointPosition(1))*QuadraturPointPosition(1)/(h*h)*v_max;
      VelocityValue[1] = 0;
      //std::cout << " qp(0) = " << QuadraturPointPosition(0) << " qp(1) = " << QuadraturPointPosition(1) << std::endl;
      //std::cout << " v(0) = " << VelocityValue(0) << " v(1) = " << VelocityValue(1) << std::endl;      
      return;
   };

   // we need grid functions to first compute the controlpoint values on the boundary, so we can project them on to our system
   // means we will build a system that needed to be solved for the desired boundary values
   
   // VELOCITY
   // define rhs with the desired boundary condition values
   mfem::VectorFunctionCoefficient vfc_inlet(sdim, lambda_inlet); // function for our desired boundary condition
   mfem::LinearForm *f_bc(new mfem::LinearForm(vfes)); // define linear form for rhs
   f_bc->AddBoundaryIntegrator(new mfem::VectorBoundaryLFIntegrator(vfc_inlet),vdbc_bdr_inlet); // define integrator on desired boundary
   f_bc->Assemble(); // assemble the linear form (vector)
   // define bilinear form add the boundary, means the nurbs add the boundary
   mfem::BilinearForm a_bc(vfes); // define the bilinear form results in n x n matrix, we use the velocity finite element space
   mfem::ConstantCoefficient One_bc(1); // coefficient for the kinematic viscosity
   a_bc.AddBoundaryIntegrator(new mfem::VectorMassIntegrator(One_bc),vdbc_bdr_inlet); // bilinear form (lambda*u_vector),(v_vector))
   a_bc.Assemble(); // assemble the bilinear form (matrix)
   
   // PRESSURE
   // define rhs with the desired boundary condition values
   mfem::ConstantCoefficient pfc_outlet(p_val); // function for our desired boundary condition
   mfem::LinearForm *g_bc(new mfem::LinearForm(pfes)); // define linear form for rhs
   g_bc->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(pfc_outlet),pdbc_bdr); // define integrator on desired boundary
   g_bc->Assemble(); // assemble the linear form (vector)
   // define bilinear form add the boundary, means the nurbs add the boundary
   mfem::BilinearForm b_bc(pfes); // define the bilinear form results in n x n matrix, we use the pressure finite element space
   b_bc.AddBoundaryIntegrator(new mfem::MassIntegrator(One_bc),pdbc_bdr); // bilinear form (lambda*u_vector),(v_vector))
   b_bc.Assemble(); // assemble the bilinear form (matrix)

   // TEMPERATURE
   // define rhs with the desired boundary condition values
   mfem::ConstantCoefficient tfc_inlet(temp_1); // function for our desired boundary condition
   mfem::ConstantCoefficient tfc_walls(temp_2); // function for our desired boundary condition
   mfem::LinearForm *h_bc(new mfem::LinearForm(tfes)); // define linear form for rhs
   h_bc->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(tfc_inlet),tdbc_bdr_inlet); // define integrator on desired boundary
   h_bc->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(tfc_walls),tdbc_bdr_walls); // define integrator on desired boundary
   h_bc->Assemble(); // assemble the linear form (vector)
   // define bilinear form add the boundary, means the nurbs add the boundary
   mfem::BilinearForm d_bc(pfes); // define the bilinear form results in n x n matrix, we use the velocity finite element space
   d_bc.AddBoundaryIntegrator(new mfem::MassIntegrator(One_bc),tdbc_bdr); // bilinear form (lambda*u_vector),(v_vector))
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
   d_bc.FormLinearSystem(temp_ess_tdof_list_dummy, t_bc, *h_bc, D_BC, T_BC, H_BC); // form D_BC
   
   //std::cout << " F_BC = " << std::endl;
   //F_BC.Print(std::cout,1);

   // SOLVER
   // setup solver
   int maxIter(100); // maximal number of iterations
   double rtol(1.e-10); // convergence criteria
   double atol(1.e-10); // convergence criteria

   // setup solver
   //1
   //mfem::MINRESSolver solver;
   //2
   mfem::CGSolver bc_solver;
   mfem::Solver *bc_prec;
   bc_prec = new mfem::GSSmoother();
   bc_solver.SetPreconditioner(*bc_prec);

   bc_solver.SetAbsTol(atol);
   bc_solver.SetRelTol(rtol);
   bc_solver.SetMaxIter(maxIter);
   bc_solver.SetOperator(A_BC);
   bc_solver.SetPrintLevel(1);

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
      std::ofstream mesh_ofs("Stokes.mesh");
      mesh_ofs.precision(8);
      mesh->Print(mesh_ofs);

      std::ofstream v_bc_ofs("sol_v_bc.gf");
      v_bc_ofs.precision(8);
      v_bc.Save(v_bc_ofs);

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
      std::ofstream mesh_ofs("Stokes.mesh");
      mesh_ofs.precision(8);
      mesh->Print(mesh_ofs);

      std::ofstream p_bc_ofs("sol_p_bc.gf");
      p_bc_ofs.precision(8);
      p_bc.Save(p_bc_ofs);
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
      std::ofstream mesh_ofs("Stokes.mesh");
      mesh_ofs.precision(8);
      mesh->Print(mesh_ofs);

      std::ofstream t_bc_ofs("sol_t_bc.gf");
      t_bc_ofs.precision(8);
      t_bc.Save(t_bc_ofs);
   }
   return true;
}