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
   double v_max = 280;               // max velocity for our boundary on the inlet
   double p_val = 100;           // value for pressure boundary
   double kin_viscosity = 200;    // value for kinematic visosity
   double temp_1 = 0;               // value for temperature
   double temp_2 = 50;               // value for temperature
   double temp_diffusion_const = 0.1; // value for temperature diffusion constant coefficient
   const char *mesh_file = "../../../MA/data/quad_nurbs.mesh";  //our standard test mesh
   int ref_levels = 1;              // standard number of refinements for the mesh
   bool visualization = 1;          // bool if visualization is wanted
   mfem::Array<int> order(2);       // to store order from mesh and order elevation. the order for the finite element collections and spaces will be, velocity=order[0] + 1 + order[1], pressure=order[0] + order[1]
   order[0] = 0;                    // mesh order will be set from given mesh 
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
   mfem::FiniteElementCollection *tfec; // temperature
   vfec = new mfem::NURBSFECollection(order[0]+order[1]+1); // Pn
   pfec = new mfem::NURBSFECollection(order[0]+order[1]); // Pn-1 , pressure has one order less than velocity, we want to make an Taylor-Hood element pair
   tfec = new mfem::NURBSFECollection(order[0]+order[1]+1); // Pn
   std::cout << "velocity finite element collection Order " << vfec->GetOrder() << std::endl;
   std::cout << "pressure finite element collection Order " << pfec->GetOrder() << std::endl;
   std::cout << "temperature finite element collection Order " << tfec->GetOrder() << std::endl;

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
   mfem::NURBSExtension *tNURBSext = NULL; // temperature
   vNURBSext = new mfem::NURBSExtension(mesh->NURBSext, order[0]+order[1]+1);
   pNURBSext = new mfem::NURBSExtension(mesh->NURBSext, order[0]+order[1]);
   tNURBSext = new mfem::NURBSExtension(mesh->NURBSext, order[0]+order[1]+1);
   
   // declaration for our finite element spaces
   // here the mesh, NURBSextension, finite element collection and the space dimensions are used to define our desired fe space
   mfem::FiniteElementSpace *vfes = new mfem::FiniteElementSpace(mesh, vNURBSext, vfec, sdim); // velocity finite element space, with dimension sdim
   mfem::FiniteElementSpace *pfes = new mfem::FiniteElementSpace(mesh, pNURBSext, pfec);  // pressure finite element space, with dimension 1 (scalar field)
   mfem::FiniteElementSpace *tfes = new mfem::FiniteElementSpace(mesh, pNURBSext, tfec);  // temperature finite element space, with dimension 1 (scalar field)
   std::cout << "vfes velocity finite element space Order " << vfes->GetMaxElementOrder() << std::endl;
   std::cout << "pfes pressure finite element space Order " << pfes->GetMaxElementOrder() << std::endl;
   std::cout << "tfes temperature finite element space Order " << tfes->GetMaxElementOrder() << std::endl;
   std::cout << "Number of finite element unknowns vfes: "
     << vfes->GetTrueVSize() << std::endl;
   std::cout << "Number of finite element unknowns pfes: "
     << pfes->GetTrueVSize() << std::endl;
   std::cout << "Number of finite element unknowns tfes: "
     << tfes->GetTrueVSize() << std::endl;

   //    Determine the list of true (i.e. conforming) essential boundary dofs.
   //    The boundary conditions are defined by marking all
   //    the boundary attributes from the mesh as essential (Dirichlet) and
   //    converting them to a list of true dofs.
   
   // Dirichlet
   mfem::Array<int> vdbc_bdr(mesh->bdr_attributes.Max()); // contains the whole boundary markers
   mfem::Array<int> vdbc_bdr_noslip(mesh->bdr_attributes.Max()); // to select only the boundary markers for the no slip walls
   mfem::Array<int> vdbc_bdr_inlet(mesh->bdr_attributes.Max());   // to select only the boundary markers for the inlet
   mfem::Array<int> pdbc_bdr(mesh->bdr_attributes.Max()); // contains the whole boundary markers
   mfem::Array<int> tdbc_bdr(mesh->bdr_attributes.Max()); // contains the whole boundary markers
   mfem::Array<int> tdbc_bdr_inlet(mesh->bdr_attributes.Max());   // to select only the boundary markers for the inlet
   mfem::Array<int> tdbc_bdr_walls(mesh->bdr_attributes.Max());   // to select only the boundary markers for the walls
   mfem::Array<int> vdummy_bdr(mesh->bdr_attributes.Max()); // contains the whole boundary markers
   mfem::Array<int> pdummy_bdr(mesh->bdr_attributes.Max()); // contains the whole boundary markers
   mfem::Array<int> tdummy_bdr(mesh->bdr_attributes.Max()); // contains the whole boundary markers
   
   std::cout << " bdr attributes " << mesh->bdr_attributes.Max() << std::endl; // output the number of boundary attributes
   // Neumann
   // not used here
   
   // for poiseuille flow 
   // mark all used boundary attributes for velocity
   vdbc_bdr = 0; vdbc_bdr[2] = 1; vdbc_bdr[0] = 1; vdbc_bdr[3] = 1;

   // the boundary attribute 0 and 2 is used for the noslip condition
   // vdbc_bdr_noslip = 0; vdbc_bdr_noslip[2] = 1; vdbc_bdr_noslip[0] = 1; //not needed if already above marked 
   // the boundary attribute 3 is used for a constant velocity at the inlet
   vdbc_bdr_inlet = 0; vdbc_bdr_inlet[3] = 1; 
   
   // mark all used boundary attributes for pressure
   pdbc_bdr = 0; pdbc_bdr[1] = 0;

   // temperature
   // mark all used boundary attributes for temperature
   tdbc_bdr = 0; tdbc_bdr[2] = 1; tdbc_bdr[3] = 1; tdbc_bdr[0] = 1;
   tdbc_bdr_inlet = 0; tdbc_bdr_inlet[3] = 1;
   tdbc_bdr_walls = 0; tdbc_bdr_walls[2] = 1; tdbc_bdr_walls[0] = 1;

   // mark dummy
   vdummy_bdr = 0;
   pdummy_bdr = 0;
   tdummy_bdr = 0;

   // get the true dofs of the boundaries
   mfem::Array<int> vel_ess_tdof_list;
   mfem::Array<int> pres_ess_tdof_list;
   mfem::Array<int> temp_ess_tdof_list;
   mfem::Array<int> vel_ess_tdof_list_dummy;
   mfem::Array<int> pres_ess_tdof_list_dummy;
   mfem::Array<int> temp_ess_tdof_list_dummy;
   
   vfes->GetEssentialTrueDofs(vdbc_bdr, vel_ess_tdof_list);
   pfes->GetEssentialTrueDofs(pdbc_bdr, pres_ess_tdof_list);   
   tfes->GetEssentialTrueDofs(tdbc_bdr, temp_ess_tdof_list);   

   vfes->GetEssentialTrueDofs(vdummy_bdr, vel_ess_tdof_list_dummy);
   pfes->GetEssentialTrueDofs(pdummy_bdr, pres_ess_tdof_list_dummy);   
   tfes->GetEssentialTrueDofs(tdummy_bdr, temp_ess_tdof_list_dummy);

   // to set our boundary conditions we first ned to define our grid functions, so that we have something to project onto
   // we need Blockoperators to define the equation system
   // Implement Blockoperators
   mfem::Array<int> block_offsets(4); // number of variables + 1
   block_offsets[0] = 0;
   block_offsets[1] = vfes->GetVSize();
   block_offsets[2] = pfes->GetVSize();
   block_offsets[3] = tfes->GetVSize();
   block_offsets.PartialSum();

   std::cout << "***********************************************************\n";
   std::cout << "dim(v) = " << block_offsets[1] - block_offsets[0] << "\n";
   std::cout << "dim(p) = " << block_offsets[2] - block_offsets[1] << "\n";
   std::cout << "dim(t) = " << block_offsets[3] - block_offsets[2] << "\n";
   std::cout << "dim(v+p) = " << block_offsets.Last() << "\n";
   std::cout << "***********************************************************\n" << std::endl;

   mfem::BlockVector x(block_offsets), rhs(block_offsets); // blockvector for gridfunctions and our rhs
   /*
      x = [ v ]      rhs = [ f ]
          [ p ]            [ g ]
          [ t ]            [ h ]
      
   */

   x = 0.0;
   rhs = 0.0;

   // initial guess set to be exact solution
   mfem::GridFunction v, p, t;
   v.MakeRef(vfes, x.GetBlock(0), 0);
   p.MakeRef(pfes, x.GetBlock(1), 0);
   t.MakeRef(tfes, x.GetBlock(2), 0);
   
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
   a_bc.FormLinearSystem(vel_ess_tdof_list_dummy, v, *f_bc, A_BC, V_BC, F_BC); // form A_BC
   b_bc.SetDiagonalPolicy(mfem::Matrix::DiagonalPolicy::DIAG_ZERO); // important, otherwise a different policy will be used, which results in false building of our matrix
   b_bc.FormLinearSystem(pres_ess_tdof_list_dummy, p, *g_bc, B_BC, P_BC, G_BC); // form B_BC
   d_bc.SetDiagonalPolicy(mfem::Matrix::DiagonalPolicy::DIAG_ZERO); // important, otherwise a different policy will be used, which results in false building of our matrix
   d_bc.FormLinearSystem(temp_ess_tdof_list_dummy, t, *h_bc, D_BC, T_BC, H_BC); // form D_BC
   
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
      v.Save(v_bc_ofs);

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
      p.Save(p_bc_ofs);
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
      t.Save(t_bc_ofs);
   }

   //std::cout << "v " << v << "\n";
   
   //std::cout << "p " << p << "\n";
   
   //std::cout << "t " << t << "\n";

   // set all values which are too small zero
   /*
   for (size_t i = 0; i < v.Size(); i++)
   {
      if (std::abs(v(i))<1.e-6)
      {
         v(i)=0;
      }
   }
   */

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
   
   // rhs for advection diffusion heat transfer
   //mfem::ConstantCoefficient zero(0.0); // zero source term
   mfem::LinearForm *h(new mfem::LinearForm); // define linear form for rhs
   h->Update(tfes, rhs.GetBlock(2), 0); // link to block vector and use the temperature finite element space
   h->AddDomainIntegrator(new mfem::DomainLFIntegrator(zero)); // define integrator for source term -> zero in our case
   h->Assemble(); // assemble the linear form (vector)

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

   // advection diffusion heat transfer
   mfem::BilinearForm d(tfes); // define the bilinear form results in n x n matrix, we use the temperature finite element space
   mfem::ConstantCoefficient temp_dcoeff(temp_diffusion_const); // coefficient for the temp_diffusion_const
   mfem::VectorGridFunctionCoefficient v_coef;
   v_coef.SetGridFunction(&v);
   //d.AddDomainIntegrator(new mfem::DiffusionIntegrator(temp_dcoeff)); // bilinear form (lambda*nabla(u),nabla(v))
   d.AddDomainIntegrator(new mfem::ConvectionIntegrator(v_coef,1)); // 
   //d.AddDomainIntegrator(new mfem::MixedDirectionalDerivativeIntegrator(v_coef)); // 
   d.Assemble(); // assemble the bilinear form (matrix)
   //a.Finalize(); not needed, will be called on form linear system

   // we need some SparseMatrix and Vector to form our linear system
   mfem::SparseMatrix A,B,C,D;
   mfem::Vector V, F;
   mfem::Vector P, G;
   mfem::Vector T, H;
   a.SetDiagonalPolicy(mfem::Matrix::DiagonalPolicy::DIAG_ZERO); // important, otherwise a different policy will be used, which results in false building of our matrix
   a.FormLinearSystem(vel_ess_tdof_list, v, *f, A, V, F); // form A   
   b.FormRectangularLinearSystem(pres_ess_tdof_list, vel_ess_tdof_list, p, *f, B, P, F); // form B
   c.FormRectangularLinearSystem(vel_ess_tdof_list, pres_ess_tdof_list, v, *g, C, V, G); // form C
   d.SetDiagonalPolicy(mfem::Matrix::DiagonalPolicy::DIAG_ZERO); // important, otherwise a different policy will be used, which results in false building of our matrix
   d.FormLinearSystem(temp_ess_tdof_list, t, *h, D, T, H); // form D
   
   /*
   a.FormSystemMatrix(vel_ess_tdof_list, A);
   b.FormRectangularSystemMatrix(pres_ess_tdof_list, vel_ess_tdof_list, B); // form B
   c.FormRectangularSystemMatrix(vel_ess_tdof_list, pres_ess_tdof_list, C); // form C
   d.FormSystemMatrix(temp_ess_tdof_list, D); // form D
   */

   // Setup stokes operator
   /*
      S = [ A    B    0] [ u ] = [ f ]
          [ C    0    0] [ p ] = [ g ]
          [ 0    0    D] [ t ] = [ h ]
   */
      
   // for this code f,g = 0

   mfem::BlockOperator stokesOp(block_offsets); // Block operator to build our System for the solver

   stokesOp.SetBlock(0,0,&A);
   stokesOp.SetBlock(0,1,&B);
   stokesOp.SetBlock(1,0,&C);
   stokesOp.SetBlock(2,2,&D);
   
   mfem::StopWatch chrono; // stop watch to calc solve time
   chrono.Clear();
   chrono.Start();
   
   // 9. Setup a very simple block diagonal preconditioner
   /*
      P = [ diag(A)             0        0]
          [ 0                   I        0]
          [ 0                   0        I]
   */


   mfem::BilinearForm a_prec(vfes); // define the bilinear form results in n x n matrix, we use the velocity finite element space
   a_prec.AddDomainIntegrator(new mfem::VectorDiffusionIntegrator(kin_vis,sdim)); // bilinear form (lambda*nabla(u_vector),nabla(v_vector))
   a_prec.Assemble(); // assemble the bilinear form (matrix)
   mfem::SparseMatrix A_DSmoother;
   a_prec.SetDiagonalPolicy(mfem::Matrix::DiagonalPolicy::DIAG_KEEP);
   a_prec.FormSystemMatrix(vel_ess_tdof_list_dummy, A_DSmoother);
   
   /*
   mfem::BilinearForm d_prec(tfes); // define the bilinear form results in n x n matrix, we use the temperature finite element space
   d_prec.AddDomainIntegrator(new mfem::DiffusionIntegrator(temp_dcoeff)); // bilinear form (lambda*nabla(u),nabla(v))
   d_prec.AddDomainIntegrator(new mfem::ConvectionIntegrator(v_coef,1)); // 
   d_prec.Assemble(); // assemble the bilinear form (matrix)
   mfem::SparseMatrix D_DSmoother;
   d_prec.SetDiagonalPolicy(mfem::Matrix::DiagonalPolicy::DIAG_KEEP);
   d_prec.FormSystemMatrix(temp_ess_tdof_list, D_DSmoother); // form D
   */

   mfem::BlockDiagonalPreconditioner stokesPrec(block_offsets);
   mfem::Solver *invM1,*invM2;

   invM1 = new mfem::DSmoother(A_DSmoother);
   invM1->iterative_mode = false;

   //invM2 = new mfem::DSmoother(D_DSmoother);
   //invM2->iterative_mode = false;

   mfem::Vector id_vec(b.Width());
   id_vec = 1.0;
   mfem::SparseMatrix id_mat(id_vec);

   mfem::Vector id_vec2(d.Width());
   id_vec2 = 1.0;
   mfem::SparseMatrix id_mat2(id_vec2);

   stokesPrec.SetDiagonalBlock(0, invM1);
   stokesPrec.SetDiagonalBlock(1, &id_mat);
   stokesPrec.SetDiagonalBlock(2, &id_mat2);

   //Wrapper to ensure mean of pressure is zero after application of preconditioner
   //mfem::BlockOrthoSolver stokesPrec_wrap(block_offsets);
   //stokesPrec_wrap.SetSolver(stokesPrec);

   // SOLVER
   // setup solver
   double lmaxIter=1000; // maximal number of iterations
   double lrtol(1.e-10); // convergence criteria
   double latol(1.e-10); // convergence criteria

   // setup solver
   // 1
   mfem::MINRESSolver lsolver; 
   
   lsolver.SetAbsTol(latol);
   lsolver.SetRelTol(lrtol);
   lsolver.SetMaxIter(lmaxIter);
   //lsolver.SetPreconditioner(stokesPrec_wrap);
   // operator
   lsolver.SetOperator(stokesOp);
   lsolver.SetPrintLevel(1);
   
   // solve the system
   // Wrapper to ensure mean of pressure is zero after each iteration of MINRES
   //mfem::BlockOrthoSolver solver_wrap(block_offsets);
   //solver_wrap.SetSolver(lsolver);
   //solver_wrap.Mult(rhs, x);
   
   //lsolver.Mult(rhs, x);
      
   for (size_t i = 0; i < 5; i++)
   {
      //H.Print(std::cout,1);
      lsolver.SetPrintLevel(0);
      lsolver.Mult(rhs, x);
      //update D
      mfem::GridFunction t_copy(t);
      //std::cout << "t " << t << "\n";
      //H.Print(std::cout,1);     
      std::cout << "Update D for the " << i << "th time \n";
      H = 0.0;
      D = 0.0;
      h->Update();
      h->Assemble(); // assemble the linear form (vector)
      d.Update();
      d.Assemble();
      d.FormLinearSystem(temp_ess_tdof_list, t, *h, D, T, H); // form D
      //H.Print(std::cout,1);
      //std::cout << "t " << t << "\n";
      t=t_copy;
      std::cout << "******** t " << t << "\n";
   }
   chrono.Stop();

   // check if solver converged
   if (lsolver.GetConverged())
   {
      std::cout << "Solver converged in " << lsolver.GetNumIterations()
                << " iterations with a residual norm of "
                << lsolver.GetFinalNorm() << ".\n";
   }
   else
   {
      std::cout << "Solver did not converge in " << lsolver.GetNumIterations()
                << " iterations. Residual norm is " << lsolver.GetFinalNorm()
                << ".\n";
   }
   std::cout << "Solver took " << chrono.RealTime() << "s.\n";

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

      std::ofstream t_ofs("sol_t.gf");
      t_ofs.precision(8);
      t.Save(t_ofs);
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
      mfem::socketstream t_sock(vishost, visport);
      t_sock.precision(8);
      t_sock << "solution\n" << mesh << t << "window_title 'Temperature'" << std::endl;
   }

   // print solution vectors
   //std::cout << "v " << v << "\n";
   //std::cout << "p " << p << "\n";
   
   return 0;
}