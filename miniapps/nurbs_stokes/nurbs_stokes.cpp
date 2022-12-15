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
#include <fstream>  // fstream for input and output in textfiles
#include <iostream> // for output to the console

#include <algorithm> // to get access to math functions
#include <cmath>     // to get access to math functions

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
   int ref_levels = 1;              // standard number of refinements for the mesh
   bool visualization = 1;          // bool if visualization is wanted
   mfem::Array<int> order(2);       // to store order from mesh and order elevation. the order for the finite element collections and spaces will be, velocity=order[0] + 1 + order[1], pressure=order[0] + order[1]
   order[0] = 0;                    // mesh order 
   order[1] = 0;                    // order elevation, Pn+1Pn finite element space is desired, higher order means higher order nurbs for the simulation

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

   // Read the mesh from the given mesh file. 
   mfem::Mesh *mesh = new mfem::Mesh(mesh_file, 1, 1);
   int sdim = mesh->SpaceDimension(); // get dimension in space from the mesh 1D, 2D, 3D
   order[0] = mesh->NURBSext->GetOrder();
   // Refine the mesh to increase the resolution.
   //    higher ref_levels means more refinement and more resolution
   {
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
   }
      
   std::cout << "Using Nurbs mesh." << std::endl;
   std::cout << "mesh Space Dimension " << sdim << std::endl;
   std::cout << "given Mesh Order " << order[0] << std::endl;
   std::cout << "given Order elevation " << order[1] << std::endl;
   // collection of finite elements from the same family (H1,RT,L2, nurbs,...) in multiple dimensions,
   // this is used to match the dof's of a finite element space between elements
   // https://docs.mfem.org/html/fe__coll_8hpp_source.html
   // we need an NURBSFECollection for every different component in our PDE, e.g.: velocity, pressure, ...
   mfem::FiniteElementCollection *vfec; // velocity
   mfem::FiniteElementCollection *pfec; // pressure
   vfec = new mfem::NURBSFECollection(order[0]+order[1]+1); // Pn
   pfec = new mfem::NURBSFECollection(order[0]+order[1]); // Pn-1 , pressure has one order less than velocity, we want to make an Taylor-Hood element pair
   std::cout << "velocity finite element collection Order " << vfec->GetOrder() << std::endl;
   std::cout << "pressure finite element collection Order " << pfec->GetOrder() << std::endl;

   /* left over from the nurbs miniapp, should be of no use here, delete later
   int nkv = mesh->NURBSext->GetNKV(); // get number of knot vectors
   if (order.Size() == 1)  //Size(): Return the logical size of the array. 
   {
      int tmp = order[0];
      std::cout << "nkv " << nkv << std::endl;
      order.SetSize(nkv); // Change the logical size of the array, keep existing entries. 
      order = tmp;
   }
   if (order.Size() != nkv ) { mfem::mfem_error("Wrong number of orders set."); }
   */

   // we want to use a nurbs mesh, so we need the nurbs extension from mfem
   // for every physical field in our equations we need an extension, here we also declare the order for our nurbs
   mfem::NURBSExtension *vNURBSext = NULL; // velocity
   mfem::NURBSExtension *pNURBSext = NULL; // pressure
   vNURBSext = new mfem::NURBSExtension(mesh->NURBSext, order[0]+order[1]+1);
   pNURBSext = new mfem::NURBSExtension(mesh->NURBSext, order[0]+order[1]);

   // declaration for our finite element spaces
   // here the mesh, NURBSextension, finite element collection and the space dimensions are used to define our desired fe space
   mfem::FiniteElementSpace *vfes = new mfem::FiniteElementSpace(mesh, vNURBSext, vfec, sdim); // velocity finite element space, with dimension sdim
   mfem::FiniteElementSpace *pfes = new mfem::FiniteElementSpace(mesh, pNURBSext, pfec);  // pressure finite element space, with dimension 1 (scalar field)
   std::cout << "vfes velocity finite element space Order " << vfes->GetMaxElementOrder() << std::endl;
   std::cout << "pfes pressure finite element space Order " << pfes->GetMaxElementOrder() << std::endl;
   std::cout << "Number of finite element unknowns vfes: "
     << vfes->GetTrueVSize() << std::endl;
   std::cout << "Number of finite element unknowns pfes: "
     << pfes->GetTrueVSize() << std::endl;
   
   //    Determine the list of true (i.e. conforming) essential boundary dofs.
   //    The boundary conditions are defined by marking all
   //    the boundary attributes from the mesh as essential (Dirichlet) and
   //    converting them to a list of true dofs.
   
   // Dirichlet
   mfem::Array<int> vdbc_bdr(mesh->bdr_attributes.Max()); // contains the whole boundary markers
   mfem::Array<int> vdbc_bdr_noslip(mesh->bdr_attributes.Max()); // to select only the boundary markers for the no slip walls
   mfem::Array<int> vdbc_bdr_inlet(mesh->bdr_attributes.Max());   // to select only the boundary markers for the inlet
   mfem::Array<int> pdbc_bdr(mesh->bdr_attributes.Max()); // contains the whole boundary markers
   mfem::Array<int> vdummy_bdr(mesh->bdr_attributes.Max()); // contains the whole boundary markers
   mfem::Array<int> pdummy_bdr(mesh->bdr_attributes.Max()); // contains the whole boundary markers
   
   std::cout << " bdr attributes " << mesh->bdr_attributes.Max() << std::endl; // output the number of boundary attributes
   // Neumann
   // not used here
   
   // for poiseuille flow 
   // mark all used boundary attributes for velocity
   vdbc_bdr = 0; vdbc_bdr[2] = 1; vdbc_bdr[0] = 1; vdbc_bdr[3] = 1;

   // the boundary attribute 0 and 2 is used for the noslip condition
   //vdbc_bdr_noslip = 0; //vdbc_bdr_noslip[2] = 1; vdbc_bdr_noslip[0] = 1; not needed if already above marked 
   // the boundary attribute 3 is used for a constant velocity at the inlet
   vdbc_bdr_inlet = 0; vdbc_bdr_inlet[3] = 1; 
   
   // mark all used boundary attributes for pressure
   pdbc_bdr = 0; pdbc_bdr[1] = 1;

   // mark dummy
   vdummy_bdr = 0;
   pdummy_bdr = 0;

   // get the true dofs of the boundaries
   mfem::Array<int> vel_ess_tdof_list;
   mfem::Array<int> pres_ess_tdof_list;
   mfem::Array<int> vel_ess_tdof_list_dummy;
   mfem::Array<int> pres_ess_tdof_list_dummy;
   
   vfes->GetEssentialTrueDofs(vdbc_bdr, vel_ess_tdof_list);
   pfes->GetEssentialTrueDofs(pdbc_bdr, pres_ess_tdof_list);   

   vfes->GetEssentialTrueDofs(vdummy_bdr, vel_ess_tdof_list_dummy);
   pfes->GetEssentialTrueDofs(pdummy_bdr, pres_ess_tdof_list_dummy);   

   // to set our boundary conditions we first ned to define our grid functions, so that we have something to project onto
   // we need Blockoperators to define the equation system
   // Implement Blockoperators
   mfem::Array<int> block_offsets(3); // number of variables + 1
   block_offsets[0] = 0;
   block_offsets[1] = vfes->GetVSize();
   block_offsets[2] = pfes->GetVSize();
   block_offsets.PartialSum();

   std::cout << "***********************************************************\n";
   std::cout << "dim(v) = " << block_offsets[1] - block_offsets[0] << "\n";
   std::cout << "dim(p) = " << block_offsets[2] - block_offsets[1] << "\n";
   std::cout << "dim(v+p) = " << block_offsets.Last() << "\n";
   std::cout << "***********************************************************\n" << std::endl;

   mfem::BlockVector x(block_offsets), rhs(block_offsets); // blockvector for gridfunctions and our rhs
   /*
      x = [ v ]      rhs = [ f ]
          [ p ]            [ g ]
      
   */

   x = 0.0;
   rhs = 0.0;

   // initial guess set to be exact solution
   mfem::GridFunction v, p;
   v.MakeRef(vfes, x.GetBlock(0), 0);
   p.MakeRef(pfes, x.GetBlock(1), 0);
   
   // available lambda functions for our wanted boundary conditions
   /*
   auto lambda_noslip = [](const mfem::Vector &QuadraturPointPosition, mfem::Vector &VelocityValue) -> void
   {
      VelocityValue[0] = 0;
      VelocityValue[1] = 0;
      //std::cout << " v(0) = " << QuadraturPointPosition(0) << " v(1) = " << QuadraturPointPosition(1) << std::endl;
      //std::cout << " x(0) = " << VelocityValue(0) << " x(1) = " << VelocityValue(1) << std::endl;      
      return;
   };
   */

   auto lambda_inlet = [&v_max](const mfem::Vector &QuadraturPointPosition, mfem::Vector &VelocityValue) -> void
   {
      double h=1;
      VelocityValue[0] = 4*(h-QuadraturPointPosition(1))*QuadraturPointPosition(1)/(h*h)*v_max;
      
      //VelocityValue[0] = 0;
      //VelocityValue[1] = v_max;
      VelocityValue[1] = 0;
      std::cout << " qp(0) = " << QuadraturPointPosition(0) << " qp(1) = " << QuadraturPointPosition(1) << std::endl;
      std::cout << " v(0) = " << VelocityValue(0) << " v(1) = " << VelocityValue(1) << std::endl;      
      return;
   };

   
   // we need grid functions to first compute the controlpoint values on the boundary, so we can project them on to our system
   // means we will build a system that needed to be solved for the desired boundary values
   
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
   
  // define rhs with the desired boundary condition values
   mfem::ConstantCoefficient pfc_outlet(p_val); // function for our desired boundary condition
   mfem::LinearForm *g_bc(new mfem::LinearForm(pfes)); // define linear form for rhs
   g_bc->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(pfc_outlet),pdbc_bdr); // define integrator on desired boundary
   g_bc->Assemble(); // assemble the linear form (vector)
   // define bilinear form add the boundary, means the nurbs add the boundary
   mfem::BilinearForm b_bc(pfes); // define the bilinear form results in n x n matrix, we use the velocity finite element space
   //mfem::ConstantCoefficient One_bc(1); // coefficient for the kinematic viscosity
   b_bc.AddBoundaryIntegrator(new mfem::MassIntegrator(One_bc),pdbc_bdr); // bilinear form (lambda*u_vector),(v_vector))
   b_bc.Assemble(); // assemble the bilinear form (matrix)

   mfem::SparseMatrix A_BC, B_BC;
   mfem::Vector V_BC, F_BC;
   mfem::Vector P_BC, G_BC;
   a_bc.SetDiagonalPolicy(mfem::Matrix::DiagonalPolicy::DIAG_ZERO); // important, otherwise a different policy will be used, which results in false building of our matrix
   a_bc.FormLinearSystem(vel_ess_tdof_list_dummy, v, *f_bc, A_BC, V_BC, F_BC); // form A_BC
   b_bc.SetDiagonalPolicy(mfem::Matrix::DiagonalPolicy::DIAG_ZERO); // important, otherwise a different policy will be used, which results in false building of our matrix
   b_bc.FormLinearSystem(pres_ess_tdof_list_dummy, p, *g_bc, B_BC, P_BC, G_BC); // form B_BC
   
   //std::cout << " F_BC = " << std::endl;
   //F_BC.Print(std::cout,1);

   // SOLVER
   // setup solver
   int maxIter(100); // maximal number of iterations
   double rtol(1.e-10); // convergence criteria
   double atol(1.e-10); // convergence criteria

   // setup minres solver, should be enough for our linear system
   // without preconditioning
   mfem::MINRESSolver solver; 
   solver.SetAbsTol(atol);
   solver.SetRelTol(rtol);
   solver.SetMaxIter(maxIter);
   solver.SetOperator(A_BC);
   solver.SetPrintLevel(1);

   // solve the system
   solver.Mult(F_BC, V_BC);

   // check if solver converged
   if (solver.GetConverged())
   {
      std::cout << "MINRES converged in " << solver.GetNumIterations()
                << " iterations with a residual norm of "
                << solver.GetFinalNorm() << ".\n";
   }
   else
   {
      std::cout << "MINRES did not converge in " << solver.GetNumIterations()
                << " iterations. Residual norm is " << solver.GetFinalNorm()
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
      v.Save(v_bc_ofs);

      //std::ofstream p_ofs("sol_p.gf");
      //p_ofs.precision(8);
      //p.Save(p_ofs);
   }

   //std::cout << "v_bc " << v_bc << "\n";
   
   solver.SetOperator(B_BC);

   // solve the system
   solver.Mult(G_BC, P_BC);

   // check if solver converged
   if (solver.GetConverged())
   {
      std::cout << "MINRES converged in " << solver.GetNumIterations()
                << " iterations with a residual norm of "
                << solver.GetFinalNorm() << ".\n";
   }
   else
   {
      std::cout << "MINRES did not converge in " << solver.GetNumIterations()
                << " iterations. Residual norm is " << solver.GetFinalNorm()
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
      p.Save(p_bc_ofs);
   }



   //std::cout << "v " << v << "\n";
   
   //std::cout << "p " << p << "\n";
   

   // projecting our computed controll points onto our grid functions for our pde system

   //v.ProjectBdrCoefficient(vfc_noslip,vdbc_bdr_invert); // we project the coefficient to the desired boundary
   //std::cout << " v = " << v << std::endl;

   // !!! i need to overwork this to only set a pressure at an vertex
   //p.ProjectBdrCoefficient(pdbcCoef, pdbc_bdr_invert); // we project the coefficient to the desired boundary
   //std::cout << " p = " << p << std::endl;

   // Setup bilinear and linear forms

   // rhs of momentum equation
   // we don't have a source term in our equations, so we just make an zero source term
   mfem::Vector vzero(sdim);
   vzero = 0.;
   mfem::VectorConstantCoefficient vcczero(vzero); // coefficient for source term

   mfem::LinearForm *f(new mfem::LinearForm); // define linear form for rhs
   f->Update(vfes, rhs.GetBlock(0),0);  // link to block vector and use the velocity finite element space
   f->AddDomainIntegrator(new mfem::VectorDomainLFIntegrator(vcczero)); // define integrator for source term -> zero in our case
   f->Assemble(); // assemble the linear form (vector)
   
   // rhs for continuity equation
   mfem::ConstantCoefficient zero(0.0); // zero source term
   mfem::LinearForm *g(new mfem::LinearForm); // define linear form for rhs
   g->Update(pfes, rhs.GetBlock(1), 0); // link to block vector and use the pressure finite element space
   g->AddDomainIntegrator(new mfem::DomainLFIntegrator(zero)); // define integrator for source term -> zero in our case
   g->Assemble(); // assemble the linear form (vector)
   
   // Momentum equation
   // diffusion term
   mfem::BilinearForm a(vfes); // define the bilinear form results in n x n matrix, we use the velocity finite element space
   mfem::ConstantCoefficient kin_vis(kin_viscosity); // coefficient for the kinematic viscosity
   a.AddDomainIntegrator(new mfem::VectorDiffusionIntegrator(kin_vis,sdim)); // bilinear form (lambda*nabla(u_vector),nabla(v_vector))
   a.Assemble(); // assemble the bilinear form (matrix)
   //a.Finalize(); not needed, will be called on form linear system

   // grad pressure term
   // define the mixed bilinear form results in n x m matrix, we use the velocity finite element space as test space and the pressure space as trial space
   mfem::MixedBilinearForm b(pfes,vfes); // (trial,test)
   mfem::ConstantCoefficient minusOne(-1.0); // -1 because of the sign in the equation
   b.AddDomainIntegrator(new mfem::GradientIntegrator(minusOne)); // mixed bilinear form (lambda*nabla(u),v_vector)
   b.Assemble(); // assemble the mixed bilinear form (matrix)
   //b.Finalize(); not needed, will be called on form linear system

   // continuity term
   // define the mixed bilinear form results in n x m matrix, we use the pressure finite element space as test space and the velocity space as trial space
   mfem::MixedBilinearForm c(vfes,pfes); // (trial,test)
   mfem::ConstantCoefficient One(1.0); // +1 because of the sign in the equation
   c.AddDomainIntegrator(new mfem::VectorDivergenceIntegrator(One)); // mixed bilinear form (lambda*nabla . u_vector, v)
   c.Assemble(); // assemble the mixed bilinear form (matrix)
   //c.Finalize(); not needed, will be called on form linear system

   // we need some SparseMatrix and Vector to form our linear system
   mfem::SparseMatrix A,B,C;
   mfem::Vector V, F;
   mfem::Vector P, G;
   a.SetDiagonalPolicy(mfem::Matrix::DiagonalPolicy::DIAG_ZERO); // important, otherwise a different policy will be used, which results in false building of our matrix
   a.FormLinearSystem(vel_ess_tdof_list, v, *f, A, V, F); // form A
   b.FormRectangularLinearSystem(pres_ess_tdof_list, vel_ess_tdof_list, p, *f, B, P, G); // form B
   c.FormRectangularLinearSystem(vel_ess_tdof_list, pres_ess_tdof_list, v, *g, C, V, F); // form C

   // Setup stokes operator
   /*
      S = [ A    B ] [ u ] = [ f ]
          [ C    0 ] [ p ] = [ g ]
   */

   // for this code f,g = 0

   mfem::BlockOperator stokesOp(block_offsets); // Block operator to build our System for the solver

   stokesOp.SetBlock(0,0,&A);
   stokesOp.SetBlock(0,1,&B);
   stokesOp.SetBlock(1,0,&C);

   
   mfem::StopWatch chrono; // stop watch to calc solve time
   chrono.Clear();
   chrono.Start();
   
   
   // SOLVER
   // setup solver
   maxIter=1000; // maximal number of iterations
   /*double rtol(1.e-10); // convergence criteria
   double atol(1.e-10); // convergence criteria

   // setup minres solver, should be enough for our linear system
   // without preconditioning
   mfem::MINRESSolver solver; 
   solver.SetAbsTol(atol);
   solver.SetRelTol(rtol);
   solver.SetMaxIter(maxIter);
   solver.SetOperator(stokesOp);
   solver.SetPrintLevel(1);
   */
   solver.SetMaxIter(maxIter);
   solver.SetOperator(stokesOp);

   // solve the system
   solver.Mult(rhs, x);
   chrono.Stop();

   // check if solver converged
   if (solver.GetConverged())
   {
      std::cout << "MINRES converged in " << solver.GetNumIterations()
                << " iterations with a residual norm of "
                << solver.GetFinalNorm() << ".\n";
   }
   else
   {
      std::cout << "MINRES did not converge in " << solver.GetNumIterations()
                << " iterations. Residual norm is " << solver.GetFinalNorm()
                << ".\n";
   }
   std::cout << "MINRES solver took " << chrono.RealTime() << "s.\n";

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
   
   // stream solution to glvis
   // doesn't work properly right now
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      mfem::socketstream v_sock(vishost, visport);
      v_sock.precision(8);
      v_sock << "solution\n" << mesh << v << "window_title 'Velocity'" << std::endl;
      mfem::socketstream p_sock(vishost, visport);
      p_sock.precision(8);
      p_sock << "solution\n" << mesh << p << "window_title 'Pressure'" << std::endl;
   }

   // print solution vectors
   //std::cout << "v " << v << "\n";
   //std::cout << "p " << p << "\n";
   
   return 0;
}