
#include "nurbs_stokes_weak_bc.hpp"
#include "mfem.hpp"
#include <fstream>
#include <iostream>

#include <algorithm>
#include <cmath>

//using namespace std; // bad idea!!!!
//using namespace mfem;

int main(int argc, char *argv[])
{
   // 0. Setup
   double v_max = 28;
   double pdbc_val = 100;
   double kin_viscosity = 1;
   double sigma = 1.0;
   double kappa = 10;
   // 1. Parse command-line options.
   //const char *mesh_file = "../../../MA/data/pipe-nurbs-boundary-test_2.mesh";
   //const char *mesh_file = "../../../MA/data/quad_nurbs_2.mesh";
   const char *mesh_file = "../../../MA/data/quad_nurbs.mesh";

   int ref_levels = 0;
   bool visualization = 1;
   mfem::Array<int> order(1);
   order[0] = 3;

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

   // 2. Read the mesh from the given mesh file. We can handle nurbs meshes with
   //    the code.
   mfem::Mesh *mesh = new mfem::Mesh(mesh_file, 1, 1);
   int sdim = mesh->Dimension();

   // 3. Refine the mesh to increase the resolution. In this example we do
   //    'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
   //    largest number that gives a final mesh with no more than 50,000
   //    elements.
   //    higher ref_levels means more refinement and more resolution
   {
      if (ref_levels < 0)
      {
         ref_levels =
            (int)floor(log(500./mesh->GetNE())/log(2.)/sdim);
      }

      for (int l = 0; l < ref_levels; l++)
      {
         mesh->UniformRefinement();
      }
      mesh->PrintInfo();
   }

   // 4. Define a finite element space on the mesh. Here we use use an isogeometric space.

   mfem::NURBSExtension *vNURBSext = NULL;
   mfem::NURBSExtension *pNURBSext = NULL;
   
   std::cout << "Using Nurbs mesh." << std::endl;
   std::cout << "mesh Dimension " << sdim << std::endl;
   std::cout << "given Order " << order[0] << std::endl;
   // we need an NURBSFECollection for every different component in our PDE, e.g.: velocity, pressure, ...
   mfem::FiniteElementCollection *vfec; // velocity
   mfem::FiniteElementCollection *pfec; // pressure
   vfec = new mfem::NURBSFECollection(order[0]);
   pfec = new mfem::NURBSFECollection(order[0]-1);
   std::cout << "velocity fec Order " << vfec->GetOrder() << std::endl;
   std::cout << "pressure fec Order " << pfec->GetOrder() << std::endl;

   int nkv = mesh->NURBSext->GetNKV();
   if (order.Size() == 1)
   {
      int tmp = order[0];
      order.SetSize(nkv);
      order = tmp;
   }

   if (order.Size() != nkv ) { mfem::mfem_error("Wrong number of orders set."); }
   vNURBSext = new mfem::NURBSExtension(mesh->NURBSext, order[0]);
   pNURBSext = new mfem::NURBSExtension(mesh->NURBSext, order[0]-1);

   mfem::FiniteElementSpace *vfes = new mfem::FiniteElementSpace(mesh, vNURBSext, vfec, sdim);
   mfem::FiniteElementSpace *pfes = new mfem::FiniteElementSpace(mesh, pNURBSext, pfec);
   std::cout << "vfes Order " << vfes->GetMaxElementOrder() << std::endl;
   std::cout << "pfes Order " << pfes->GetMaxElementOrder() << std::endl;
   std::cout << "Number of finite element unknowns vfes: "
     << vfes->GetTrueVSize() << std::endl;
   std::cout << "Number of finite element unknowns pfes: "
     << pfes->GetTrueVSize() << std::endl;
   
   
   // 5. Determine the list of true (i.e. conforming) essential boundary dofs.
   //    The boundary conditions are defined by marking all
   //    the boundary attributes from the mesh as essential (Dirichlet) and
   //    converting them to a list of true dofs.
   
   // Dirichlet
   mfem::Array<int> vdbc_bdr(mesh->bdr_attributes.Max());
   mfem::Array<int> vdbc_bdr_noslip(mesh->bdr_attributes.Max());
   mfem::Array<int> vdbc_bdr_inlet(mesh->bdr_attributes.Max());
   mfem::Array<int> pdbc_bdr(mesh->bdr_attributes.Max());
   std::cout << " bdr attributes " << mesh->bdr_attributes.Max() << std::endl;
   // Neumann
   //Array<int> vnbc_bdr(mesh->bdr_attributes.Max());
   //Array<int> pnbc_bdr(mesh->bdr_attributes.Max());
   
   // we assume that the boundary attribute 1,2 are dirchlet, we want to use either pressure or velocity or a mix
   //vdbc_bdr = 0; vdbc_bdr[0] = 1; vdbc_bdr[2] = 1;
   vdbc_bdr = 0; vdbc_bdr[2] = 1; vdbc_bdr[0] = 1; vdbc_bdr[3] = 1;  
   vdbc_bdr_noslip = 0; vdbc_bdr_noslip[2] = 1; vdbc_bdr_noslip[0] = 1; 
   vdbc_bdr_inlet = 0; vdbc_bdr_inlet[3] = 1;  // poiseuille flow
   
   pdbc_bdr = 0; pdbc_bdr[1] = 1;

   mfem::Array<int> vel_ess_tdof_list;
   mfem::Array<int> pres_ess_tdof_list;
   //vfes->GetEssentialTrueDofs(vdbc_bdr, vel_ess_tdof_list);
   vfes->GetEssentialTrueDofs(vdbc_bdr, vel_ess_tdof_list);
   pfes->GetEssentialTrueDofs(pdbc_bdr, pres_ess_tdof_list);   

   // 7. Implement Blockoperators
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

   mfem::BlockVector x(block_offsets), rhs(block_offsets);
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


   auto lambda_noslip = [&sdim](const mfem::Vector &QuadraturPointPosition, mfem::Vector &VelocityValue) -> void
   {
      VelocityValue[0] = 0;
      VelocityValue[1] = 0;
      
      std::cout << " v(0) = " << QuadraturPointPosition(0) << " v(1) = " << QuadraturPointPosition(1) << std::endl;
      std::cout << " x(0) = " << VelocityValue(0) << " x(1) = " << VelocityValue(1) << std::endl;      
      return;
   };

   auto lambda_inlet = [&v_max, &sdim](const mfem::Vector &QuadraturPointPosition, mfem::Vector &VelocityValue) -> void
   {
      /*
      v=4(h-y)y/h²*v_max
      */
      double h=1;
      //VelocityValue[0] = 4*(h-QuadraturPointPosition(1))*QuadraturPointPosition(1)/(h*h)*v_max;
      VelocityValue[0] = v_max;
      VelocityValue[1] = 0;
      
      std::cout << " qp(0) = " << QuadraturPointPosition(0) << " qp(1) = " << QuadraturPointPosition(1) << std::endl;
      std::cout << " v(0) = " << VelocityValue(0) << " v(1) = " << VelocityValue(1) << std::endl;      
      return;
   };

   mfem::VectorFunctionCoefficient vfc_noslip(sdim, lambda_noslip);
   mfem::VectorFunctionCoefficient vfc_inlet(sdim, lambda_inlet);
   v.ProjectBdrCoefficient(vfc_noslip,vdbc_bdr_noslip);
   v.ProjectBdrCoefficient(vfc_inlet,vdbc_bdr_inlet);
   std::cout << " v = " << v << std::endl;

   //ConstantCoefficient vdbcCoef(vdbc_val_2_1);
   //v.ProjectCoefficient(vdbcCoef, vdbc_bdr);
   //std::cout << " p = " << p << std::endl;
   mfem::ConstantCoefficient pdbcCoef(pdbc_val);
   p.ProjectBdrCoefficient(pdbcCoef, pdbc_bdr);
   std::cout << " p = " << p << std::endl;

   // Setup bilinear and linear forms

   // rhs of momentum equation
   // LinearForm f(vfes);  
   mfem::Vector vzero(sdim);
   vzero = 0.;
   mfem::VectorConstantCoefficient vcczero(vzero);

   mfem::LinearForm *f(new mfem::LinearForm);
   f->Update(vfes, rhs.GetBlock(0),0);   
   f->AddDomainIntegrator(new mfem::VectorDomainLFIntegrator(vcczero));
   f->AddBdrFaceIntegrator(new VectorDGDirichletLFIntegrator(vfc_noslip,sigma,kappa,sdim));
   f->AddBdrFaceIntegrator(new VectorDGDirichletLFIntegrator(vfc_inlet,sigma,kappa,sdim));
   f->Assemble();
   
   // rhs for continuity equation
   mfem::ConstantCoefficient zero(0.0);
   mfem::LinearForm *g(new mfem::LinearForm);
   g->Update(pfes, rhs.GetBlock(1), 0);
   g->AddDomainIntegrator(new mfem::DomainLFIntegrator(zero));
   g->Assemble();
   
   // Momentum equation
   // diffusion term
   mfem::BilinearForm a(vfes);
   mfem::ConstantCoefficient kin_vis(kin_viscosity);
   a.AddDomainIntegrator(new mfem::VectorDiffusionIntegrator(kin_vis,sdim));
   a.AddBdrFaceIntegrator(new VectorDGDiffusionIntegrator(sigma,kappa,sdim));
   a.Assemble();
   a.Finalize();

   // grad pressure term
   mfem::MixedBilinearForm b(pfes,vfes); // (trial,test)
   mfem::ConstantCoefficient minusOne(-1.0);
   b.AddDomainIntegrator(new mfem::GradientIntegrator(minusOne));
   b.AddBdrFaceIntegrator(new DGAvgNormalJumpIntegrator(sdim));
   b.Assemble();
   b.Finalize();

   // continuity term
   mfem::MixedBilinearForm c(vfes,pfes); // (trial,test)
   mfem::ConstantCoefficient One(1.0);
   c.AddDomainIntegrator(new mfem::VectorDivergenceIntegrator(One));
   c.Assemble();
   c.Finalize();

   //std::cout << " v = " << v << std::endl;

   /*
   mfem::SparseMatrix A,B,C;
   mfem::Vector V, F;
   mfem::Vector P, G;
   a.SetDiagonalPolicy(mfem::Matrix::DiagonalPolicy::DIAG_ZERO);
   //a.EliminateVDofs(vel_ess_tdof_list, diag);
   //a.Finalize();
   //a.FormSystemMatrix(vel_ess_tdof_list, A);
   //b.FormRectangularSystemMatrix(pres_ess_tdof_list, vel_ess_tdof_list, B);
   a.FormLinearSystem(vel_ess_tdof_list, v, *f, A, V, F);
   b.FormRectangularLinearSystem(pres_ess_tdof_list, vel_ess_tdof_list, p, *f, B, P, G);
   c.FormRectangularLinearSystem(vel_ess_tdof_list, pres_ess_tdof_list, v, *g, C, V, F);
   */

   // Setup stokes operator
   /*
      S = [ A    B ] [ u ] = [ f ]
          [ C    0 ] [ p ] = [ g ]
   */

   // for this code f,g = 0

   mfem::BlockOperator stokesOp(block_offsets);
   //mfem::SparseMatrix A,B,C;
   mfem::SparseMatrix &A(a.SpMat());
   mfem::SparseMatrix &B(b.SpMat());
   mfem::SparseMatrix &C(c.SpMat());
 

   //A.PrintInfo(std::cout);
   //A.PrintMatlab(std::cout);
   //B.PrintInfo(std::cout);
   //B.PrintMatlab(std::cout);

   stokesOp.SetBlock(0,0,&A);
   stokesOp.SetBlock(0,1,&B);
   stokesOp.SetBlock(1,0,&C);
   
   //stokesOp.SetBlock(1,0,&C);

   
   // 9. SOLVER
   // setup solver
   int maxIter(1000);
   double rtol(1.e-12);
   double atol(1.e-12);

   mfem::StopWatch chrono;
   chrono.Clear();
   chrono.Start();
   
   /*
   mfem::Solver *J_solver;
   mfem::Solver *J_prec;
   J_prec = new mfem::DSmoother(1);
   mfem::MINRESSolver *J_minres = new mfem::MINRESSolver;
   J_minres->SetRelTol(rtol);
   J_minres->SetAbsTol(0.0);
   J_minres->SetMaxIter(10);
   J_minres->SetPrintLevel(-1);
   J_minres->SetPreconditioner(*J_prec);
   J_solver = J_minres;
   
   mfem::NewtonSolver solver;
   solver.SetSolver(*J_solver);
   solver.SetAbsTol(atol);
   solver.SetRelTol(rtol);
   solver.SetMaxIter(maxIter);
   solver.SetOperator(stokesOp);
   solver.SetPrintLevel(1);
   */
   
   mfem::MINRESSolver solver;
   solver.SetAbsTol(atol);
   solver.SetRelTol(rtol);
   solver.SetMaxIter(maxIter);
   solver.SetOperator(stokesOp);
   solver.SetPrintLevel(1);
   

   // solve the system
   solver.Mult(rhs, x);
   chrono.Stop();

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
   
   std::cout << " v.sol = " << v << std::endl;
   std::cout << " p.sol = " << p << std::endl;
   std::cout << " vfes = " << vfes->GetTrueVSize() << std::endl;
   std::cout << " pfes = " << pfes->GetTrueVSize() << std::endl;

   //F.Print(std::cout);
   //std::cout << " RHS F= " << std::endl;
   //rhs.GetBlock(0).Print(std::cout);
   //std::cout << " RHS G= " << std::endl;
   //rhs.GetBlock(1).Print(std::cout);

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


   // 15. Free the used memory.
   /*delete a;
   delete A;
   delete b;
   delete B;
   delete vfes;
   delete pfes;
   delete mesh;
   delete pres_ess_tdof_list;
   delete vel_ess_tdof_list;*/
   vel_ess_tdof_list = 0;
   pres_ess_tdof_list = 0;


   /*   
   A.PrintInfo(std::cout);
   A.PrintMatlab(std::cout);
   B.PrintInfo(std::cout);
   B.PrintMatlab(std::cout);
   */
   return 0;
}