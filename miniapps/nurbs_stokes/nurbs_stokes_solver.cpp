#include "nurbs_stokes_solver.hpp"
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
      std::cout << "given Order elevation pressure " << order[2] << std::endl;
      std::cout << "given Order elevation temperature " << order[3] << std::endl;
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
   MFEM_ASSERT(user_tdbc_bdr.Size() == user_tdbc_bdr_values.Size(), "number of given temperature boundaries and values does'nt match!");

   // init our boundary markers
   vdbc_bdr_all = mfem::Array<int>(mesh->bdr_attributes.Max()); // contains the whole boundary markers
   vdbc_bdr = mfem::Array<int>(mesh->bdr_attributes.Max());   // to select only the set boundary markers
   vdbc_bdr_noslip = mfem::Array<int>(mesh->bdr_attributes.Max()); // to select only the boundary markers for the no slip walls
   pdbc_bdr = mfem::Array<int>(mesh->bdr_attributes.Max()); // contains the whole boundary markers
   tdbc_bdr = mfem::Array<int>(mesh->bdr_attributes.Max()); // contains the whole boundary markers
   vdummy_bdr = mfem::Array<int>(mesh->bdr_attributes.Max()); // contains the whole boundary markers
   pdummy_bdr = mfem::Array<int>(mesh->bdr_attributes.Max()); // contains the whole boundary markers
   tdummy_bdr = mfem::Array<int>(mesh->bdr_attributes.Max()); // contains the whole boundary markers

   //set marker
   // velocity
   vdbc_bdr = 0;
   vdbc_bdr_all = 0;
   vdbc_bdr_noslip = 0;

   // pressure
   pdbc_bdr = 0;

   // temperature
   tdbc_bdr = 0;

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
   for (size_t i = 0; i < user_tdbc_bdr.Size(); i++)
   {
      tdbc_bdr[user_tdbc_bdr[i]-1] = 1;
   }
   
   // get the true dofs of the boundaries
   vfes->GetEssentialTrueDofs(vdbc_bdr_all, vel_ess_tdof_list);
   pfes->GetEssentialTrueDofs(pdbc_bdr, pres_ess_tdof_list);   
   tfes->GetEssentialTrueDofs(tdbc_bdr, temp_ess_tdof_list);

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
   
bool NurbsStokesSolver::set_dirichletbc_temperature(mfem::Array<int> boundaries,mfem::Vector boundaries_values)
{
   user_tdbc_bdr = boundaries;
   user_tdbc_bdr_values = boundaries_values;

   return true;
}

bool NurbsStokesSolver::calc_dirichletbc(mfem::GridFunction &v0, mfem::GridFunction &p0, mfem::GridFunction &t0)
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

   mfem::GridFunction v_bc(vfes),p_bc(pfes), t_bc(tfes); // to calculate our gridfunction on the dirichlet boundary

   // we need grid functions to first compute the controlpoint values on the boundary, so we can project them on to our system
   // means we will build a system that needed to be solved for the desired boundary values

   // create vectors for coefficients and boundary markers

   std::vector<mfem::Array<int>> vdbc_bdr_marker;
   std::vector<mfem::Vector> vdbc_bdr_vector;
   std::vector<mfem::VectorConstantCoefficient> vdbc_bdr_vectorcoefficient;

   for (size_t i = 0; i < user_vdbc_bdr.Size(); i++)
   {
      vdbc_bdr_marker.push_back(mfem::Array<int>(mesh->bdr_attributes.Max()));
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
      pdbc_bdr_marker.push_back(mfem::Array<int>(mesh->bdr_attributes.Max()));
      pdbc_bdr_marker[i] = 0;
      pdbc_bdr_marker[i][user_pdbc_bdr[i]-1] = 1;

      mfem::ConstantCoefficient cc(user_pdbc_bdr_values[i]);
      pdbc_bdr_coefficient.push_back(cc);
      //std::cout << user_pdbc_bdr_values[i] << " pressure values \n";
   }

   std::vector<mfem::Array<int>> tdbc_bdr_marker;
   std::vector<mfem::ConstantCoefficient> tdbc_bdr_coefficient;

   for (size_t i = 0; i < user_tdbc_bdr.Size(); i++)
   {
      tdbc_bdr_marker.push_back(mfem::Array<int>(mesh->bdr_attributes.Max()));
      tdbc_bdr_marker[i] = 0;
      tdbc_bdr_marker[i][user_tdbc_bdr[i]-1] = 1;

      mfem::ConstantCoefficient cc(user_tdbc_bdr_values[i]);
      tdbc_bdr_coefficient.push_back(cc);
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
   mfem::LinearForm *h_bc(new mfem::LinearForm(tfes)); // define linear form for rhs
   // define bilinear form add the boundary, means the nurbs add the boundary
   mfem::BilinearForm d_bc(tfes); // define the bilinear form results in n x n matrix, we use the velocity finite element space
   //h_bc->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(tfc_inlet),tdbc_bdr); // define integrator on desired boundary
   //h_bc->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(tfc_walls),tdbc_bdr); // define integrator on desired boundary
   for (size_t i = 0; i < tdbc_bdr_marker.size(); i++)
   {
      h_bc->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(tdbc_bdr_coefficient[i]),tdbc_bdr_marker[i]); // define integrator on desired boundary
      d_bc.AddBoundaryIntegrator(new mfem::MassIntegrator(One_bc),tdbc_bdr_marker[i]); // bilinear form (lambda*u_vector),(v_vector))
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
   d_bc.FormLinearSystem(temp_ess_tdof_list_dummy, t_bc, *h_bc, D_BC, T_BC, H_BC); // form D_BC
   
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
      std::ofstream mesh_ofs("Stokes.mesh");
      mesh_ofs.precision(8);
      mesh->Print(mesh_ofs);

      std::ofstream v_bc_ofs("sol_v_bc.gf");
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
      std::ofstream mesh_ofs("Stokes.mesh");
      mesh_ofs.precision(8);
      mesh->Print(mesh_ofs);

      std::ofstream p_bc_ofs("sol_p_bc.gf");
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
      std::ofstream mesh_ofs("Stokes.mesh");
      mesh_ofs.precision(8);
      mesh->Print(mesh_ofs);

      std::ofstream t_bc_ofs("sol_t_bc.gf");
      t_bc_ofs.precision(8);
      t_bc.Save(t_bc_ofs);
      t0.Save(t_bc_ofs);
   }
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      mfem::socketstream sol_sock1(vishost, visport);
      sol_sock1.precision(8);
      sol_sock1 << "solution\n" << *mesh << v_bc << std::flush;
      mfem::socketstream sol_sock2(vishost, visport);
      sol_sock2.precision(8);
      sol_sock2 << "solution\n" << *mesh << p_bc << std::flush;
      mfem::socketstream sol_sock3(vishost, visport);
      sol_sock3.precision(8);
      sol_sock3 << "solution\n" << *mesh << t_bc << std::flush;
   }

   // SET VALUES ON NO SLIP BOUNDARIES TO ZERO!!!
   mfem::Vector zerovector(2);
   zerovector = 0.0;
   mfem::VectorConstantCoefficient zero(zerovector);
   v_bc.ProjectBdrCoefficient(zero,vdbc_bdr_noslip);
   //

   v0=v_bc;
   p0=p_bc;
   t0=t_bc;

   return true;
}

bool NurbsStokesSolver::calc_flowsystem_strongbc(mfem::GridFunction &v0,mfem::GridFunction &p0, mfem::GridFunction &t0, mfem::GridFunction &v, mfem::GridFunction &p, mfem::GridFunction &t, mfem::Coefficient &kin_vis)
{  
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
   mfem::Vector vpressure(sdim);
   vpressure = p_val;
   mfem::VectorConstantCoefficient vccpressure(vpressure); // coefficient for source term
   
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
   if (temp_kin_vis_carreau = dynamic_cast<CarreauModelCoefficient*>(&kin_vis))
   {
      CarreauModelCoefficient kin_vis_carreau = *temp_kin_vis_carreau;
      kin_vis_carreau.SetVelocity(v0);
      a.AddDomainIntegrator(new mfem::VectorDiffusionIntegrator(kin_vis_carreau, sdim)); // bilinear form (lambda*nabla(u_vector),nabla(v_vector))
      a.Assemble();
   }else{
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
   b.FormRectangularLinearSystem(pres_ess_tdof_list_dummy, vel_ess_tdof_list, p0, *f, B, P, F); // form B
   c.FormRectangularLinearSystem(vel_ess_tdof_list, pres_ess_tdof_list_dummy, v0, *g, C, V, G); // form C

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
      std::ofstream mesh_ofs("Stokes.mesh");
      mesh_ofs.precision(8);
      mesh->Print(mesh_ofs);

      std::ofstream v_ofs("sol_v.gf");
      v_ofs.precision(8);
      v.Save(v_ofs);

      std::ofstream p_ofs("sol_p.gf");
      p_ofs.precision(8);
      p.Save(p_ofs);
   }
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      mfem::socketstream sol_sock1(vishost, visport);
      sol_sock1.precision(8);
      sol_sock1 << "solution\n" << *mesh << v << std::flush;
      mfem::socketstream sol_sock2(vishost, visport);
      sol_sock2.precision(8);
      sol_sock2 << "solution\n" << *mesh << p << std::flush;
   }

   return true;
}

bool NurbsStokesSolver::calc_temperaturesystem_strongbc(mfem::GridFunction &v0, mfem::GridFunction &t0, mfem::GridFunction &v, mfem::GridFunction &t)
{
   // Setup bilinear and linear forms
   t = t0;
/*
   auto lambda_velocityfield = [this](const mfem::Vector &QuadraturPointPosition, mfem::Vector &VelocityValue) -> void
   {
      double h=1;
      VelocityValue[0] = 4*(h-QuadraturPointPosition(1))*QuadraturPointPosition(1)/(h*h)*v_max;
      
      VelocityValue[0] = 0;
      //VelocityValue[1] = 4*(h-QuadraturPointPosition(1))*QuadraturPointPosition(1)/(h*h)*v_max;
      VelocityValue[1] = 0;
      std::cout << " qp(0) = " << QuadraturPointPosition(0) << " qp(1) = " << QuadraturPointPosition(1) << std::endl;
      std::cout << " v(0) = " << VelocityValue(0) << " v(1) = " << VelocityValue(1) << std::endl;      
      return;
   };
*/
   // rhs for advection diffusion heat transfer
   mfem::ConstantCoefficient zero(0.0); // zero source term
   mfem::LinearForm *h(new mfem::LinearForm(tfes)); // define linear form for rhs
   h->AddDomainIntegrator(new mfem::DomainLFIntegrator(zero)); // define integrator for source term -> zero in our case
   h->Assemble(); // assemble the linear form (vector)

   // advection diffusion heat transfer
   mfem::BilinearForm d(tfes); // define the bilinear form results in n x n matrix, we use the temperature finite element space
   mfem::ConstantCoefficient temp_dcoeff(temp_diffusion_const); // coefficient for the temp_diffusion_const
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
   d.FormLinearSystem(temp_ess_tdof_list, t, *h, D, T, H); // form D
 
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

   std::cout << "SOLVE TEMPERATUREFIELD \n";   
   // solve the system
   solver.Mult(H, t);
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
      std::ofstream mesh_ofs("Stokes.mesh");
      mesh_ofs.precision(8);
      mesh->Print(mesh_ofs);

      std::ofstream t_ofs("sol_t.gf");
      t_ofs.precision(8);
      t.Save(t_ofs);
   }
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      mfem::socketstream sol_sock(vishost, visport);
      sol_sock.precision(8);
      sol_sock << "solution\n" << *mesh << t << std::flush;
   }
   /*
   // SHEARRATE COMPUTATION FOR CHECKING
   mfem::GridFunction shearrate(tfes);
   shearrate = 0;
   //ShearRateCoefficient src;
   CarreauModelCoefficient src;
   src.SetA(200);
   src.SetB(1);
   src.SetC(1);
   src.SetVelocity(v0);

   mfem::LinearForm *rhs(new mfem::LinearForm(tfes)); // define linear form for rhs
   rhs->AddDomainIntegrator(new mfem::DomainLFIntegrator(src));
   rhs->Assemble(); // assemble the linear form (vector)
   mfem::BilinearForm a(tfes); // define the bilinear form results in n x n matrix, we use the velocity finite element space
   mfem::ConstantCoefficient One_bc(1); // coefficient for the kinematic viscosity
   a.AddDomainIntegrator(new mfem::MassIntegrator(One_bc)); // bilinear form (lambda*u_vector),(v_vector))
   a.Assemble(); // assemble the bilinear form (matrix)
   
   mfem::SparseMatrix A;
   mfem::Vector SHEARRATE, RHS;
   a.SetDiagonalPolicy(mfem::Matrix::DiagonalPolicy::DIAG_ZERO); // important, otherwise a different policy will be used, which results in false building of our matrix
   a.FormLinearSystem(temp_ess_tdof_list_dummy, shearrate, *rhs, A, SHEARRATE, RHS); // form D_BC

   solver.SetOperator(A);
   std::cout << "SOLVE SHEARRATE \n";   
   // solve the system
   solver.Mult(RHS, shearrate);
   {
      std::ofstream mesh_ofs("Stokes.mesh");
      mesh_ofs.precision(8);
      mesh->Print(mesh_ofs);

      std::ofstream shearrate_ofs("sol_shearrate.gf");
      shearrate_ofs.precision(8);
      shearrate.Save(shearrate_ofs);
   }
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      mfem::socketstream sol_sock(vishost, visport);
      sol_sock.precision(8);
      sol_sock << "solution\n" << *mesh << shearrate << std::flush;
   }
   */

   return true;
}


bool NurbsStokesSolver::solve_flow(mfem::GridFunction &v0,mfem::GridFunction &p0, mfem::GridFunction &t0, mfem::GridFunction &v, mfem::GridFunction &p, mfem::GridFunction &t,mfem::Coefficient &kin_vis)
{  
   if (bcstrong)
   {
      calc_flowsystem_strongbc(v0, p0, t0, v, p, t, kin_vis);
   } else if (bcweak)
   {
      /* code */
   }
      
   return true;
}

bool NurbsStokesSolver::solve_temperature(mfem::GridFunction &v0, mfem::GridFunction &t0, mfem::GridFunction &v, mfem::GridFunction &t)
{

   if (bcstrong)
   {
      calc_temperaturesystem_strongbc(v0, t0, v, t);
   } else if (bcweak)
   {
      /* code */
   }


   return true;
}
;