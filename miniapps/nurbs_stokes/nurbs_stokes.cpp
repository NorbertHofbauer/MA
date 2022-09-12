//                          MFEM - NURBS Version
//
// Compile with: make nurbs_stokes
//
// Sample runs:  nurbs_stokes -m ../../mesh/pipe-nurbs-boundary-test_2.mesh
//
// Description:  nurbs stokes solver
//
// TODO: - for building an MPI application we need to integrate the "Par" everywhere!
//       - degree elevation?



#include "mfem.hpp"
#include <fstream>
#include <iostream>

#include <algorithm>
#include <cmath>

//using namespace std; // bad idea!!!!
using namespace mfem;

int main(int argc, char *argv[])
{
   // 0. Setup
   double vdbc_val1 = 5e3;
   //double vdbc_val2 = 22;
   double pdbc_val1 = 1e5;
   //double pdbc_val2 = 22;
   double vnbc_val = 0;
   double pnbc_val = 0;

   // 1. Parse command-line options.
   const char *mesh_file = "../../../MA/mesh/pipe-nurbs-boundary-test_2.mesh";
   
   int ref_levels = 1;
   bool visualization = 1;
   Array<int> order(1);
   order[0] = 2;

   OptionsParser args(argc, argv);
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
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   int dim = mesh->Dimension();

   // 3. Refine the mesh to increase the resolution. In this example we do
   //    'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
   //    largest number that gives a final mesh with no more than 50,000
   //    elements.
   //    higher ref_levels means more refinement and more resolution
   {
      if (ref_levels < 0)
      {
         ref_levels =
            (int)floor(log(500./mesh->GetNE())/log(2.)/dim);
      }

      for (int l = 0; l < ref_levels; l++)
      {
         mesh->UniformRefinement();
      }
      mesh->PrintInfo();
   }

   // 4. Define a finite element space on the mesh. Here we use use an isogeometric space.


   NURBSExtension *vNURBSext = NULL;
   NURBSExtension *pNURBSext = NULL;
   
   std::cout << "Using Nurbs mesh." << std::endl;
   std::cout << "mesh Dimension " << dim << std::endl;
   std::cout << "given Order " << order[0] << std::endl;
   // we need an NURBSFECollection for every different component in our PDE, e.g.: velocity, pressure, ...
   FiniteElementCollection *vfec; // velocity
   FiniteElementCollection *pfec; // pressure
   vfec = new NURBSFECollection(order[0]+1);
   pfec = new NURBSFECollection(order[0]);
   std::cout << "velocity fec Order " << vfec->GetOrder() << std::endl;
   std::cout << "pressure fec Order " << pfec->GetOrder() << std::endl;

   int nkv = mesh->NURBSext->GetNKV();
   if (order.Size() == 1)
   {
      int tmp = order[0];
      order.SetSize(nkv);
      order = tmp;
   }

   if (order.Size() != nkv ) { mfem_error("Wrong number of orders set."); }
   vNURBSext = new NURBSExtension(mesh->NURBSext, order[0]+1);
   pNURBSext = new NURBSExtension(mesh->NURBSext, order[0]);

   FiniteElementSpace *vfes = new FiniteElementSpace(mesh, vNURBSext, vfec, dim);
   FiniteElementSpace *pfes = new FiniteElementSpace(mesh, pNURBSext, pfec);
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
   Array<int> vdbc_bdr(mesh->bdr_attributes.Max());
   Array<int> pdbc_bdr(mesh->bdr_attributes.Max());
   // Neumann
   Array<int> vnbc_bdr(mesh->bdr_attributes.Max());
   Array<int> pnbc_bdr(mesh->bdr_attributes.Max());
   
   // we assume that the boundary attribute 1,2 are dirchlet, we want to use either pressure or velocity or a mix
   vdbc_bdr = 0; vdbc_bdr[0] = 1; 
   //vdbc_bdr[2] = 1; vdbc_bdr[3] = 1;  //velocitiy boundary attribute 1,2,3
   pdbc_bdr = 0; pdbc_bdr[1] = 1;  //pressure boundary attribute 4
   
   // we assume that the boundary attribute 3,4 are neumann, we want to use either pressure or velocity or a mix
   //vnbc_bdr = 0; vnbc_bdr[2] = 0;  //velocitiy boundary attribute 3
   //pnbc_bdr = 0; pnbc_bdr[3] = 0;  //pressure boundary attribute 4

   // For a continuous basis the linear system must be modified to enforce an
   // essential (Dirichlet) boundary condition.
   Array<int> vel_ess_tdof_list;
   Array<int> pres_ess_tdof_list;
   vfes->GetEssentialTrueDofs(vdbc_bdr, vel_ess_tdof_list);
   pfes->GetEssentialTrueDofs(pdbc_bdr, pres_ess_tdof_list);   

   // 6. Setup the various coefficients needed for the PDE and the
   //    various boundary conditions. In general these coefficients could be
   //    functions of position or constants. For functions use lambda functions!!!

   // first try with constant coefficients

   // construct array for Coefficients to combine the boundary values
   // DIRICHLET
   // Velocity
   /*VectorArrayCoefficient vac_vdbc(2);
   Vector vec_vdbc(mesh->bdr_attributes.Max());
   //vec_vdbc = 0.0;
   vec_vdbc(0) = vdbc_val1;
   vec_vdbc(1) = 0.0;
   vec_vdbc(2) = 0.0;
   vac_vdbc.Set(0,new PWConstCoefficient(vec_vdbc));
   */


   // pressure
   
   //VectorArrayCoefficient vac_pdbc(0);
   //Vector vec_pdbc(mesh->bdr_attributes.Max());
   //vec_pdbc = 0.0;
   //vec_pdbc(0) = pdbc_val1;
   //vec_pdbc(1) = pdbc_val2;
   //vac_pdbc.Set(0,new PWConstCoefficient(vec_pdbc));
   


   // NEUMANN
   // Velocity
   //ConstantCoefficient vnbcCoef(vnbc_val);
   // pressure
   //ConstantCoefficient pnbcCoef(pnbc_val);



   // 7. Implement Blockoperators
   Array<int> block_offsets(3); // number of variables + 1
   block_offsets[0] = 0;
   block_offsets[1] = vfes->GetVSize();
   block_offsets[2] = pfes->GetVSize();
   block_offsets.PartialSum();

   std::cout << "***********************************************************\n";
   std::cout << "dim(v) = " << block_offsets[1] - block_offsets[0] << "\n";
   std::cout << "dim(p) = " << block_offsets[2] - block_offsets[1] << "\n";
   std::cout << "dim(v+p) = " << block_offsets.Last() << "\n";
   std::cout << "***********************************************************\n" << std::endl;

   BlockVector x(block_offsets), rhs(block_offsets);
   /*
      x = [ v ]      rhs = [ f ]
          [ p ]            [ g ]
      
   */


   // initial guess set to be exact solution
   GridFunction v, p;
   v.MakeRef(vfes, x.GetBlock(0), 0);
   p.MakeRef(pfes, x.GetBlock(1), 0);
   

   // Setup bilinear and linear forms

   // rhs of momentum equation
   // LinearForm f(vfes);  
   Vector vzero(dim);
   vzero = 0.;
   VectorConstantCoefficient vcczero(vzero);

   LinearForm *f(new LinearForm);
   f->Update(vfes, rhs.GetBlock(0),0);   
   f->AddDomainIntegrator(new VectorDomainLFIntegrator(vcczero));
   f->Assemble();
   
   // rhs for continuity equation
   ConstantCoefficient zero(0.0);
   LinearForm *g(new LinearForm);
   g->Update(pfes, rhs.GetBlock(1), 0);
   g->AddDomainIntegrator(new DomainLFIntegrator(zero));
   g->Assemble();
   
   // Momentum equation
   // diffusion term
   BilinearForm a(vfes);
   a.AddDomainIntegrator(new VectorDiffusionIntegrator(dim));
   a.Assemble();
   a.Finalize();

   // grad pressure term
   MixedBilinearForm b(pfes,vfes); // (trial,test)
   ConstantCoefficient minusOne(-1.0);
   b.AddDomainIntegrator(new TransposeIntegrator(new VectorDivergenceIntegrator(minusOne)));
   b.Assemble();
   b.Finalize();

   SparseMatrix A,B;
   a.FormSystemMatrix(vel_ess_tdof_list, A);
   b.FormRectangularSystemMatrix(vel_ess_tdof_list, pres_ess_tdof_list, B);


   // Setup stokes operator
   /*
      S = [ A    B ] [ u ] = [ f ]
          [ B^T  0 ] [ p ] = [ g ]
   */

   // for this code f,g = 0

   BlockOperator stokesOp(block_offsets);

   //SparseMatrix &A(a.SpMat());
   //SparseMatrix &B(b.SpMat());
   B.EnsureMultTranspose();
   TransposeOperator Bt(&B);


   stokesOp.SetBlock(0,0,&A);
   stokesOp.SetBlock(0,1,&B);
   stokesOp.SetBlock(1,0,&Bt);
   
   // 8. Define the solution vectors v,p as a finite element grid function
   //    corresponding to the fespaces. Initialize v,p with initial guess of zero,
   //    which satisfies the boundary conditions.
   //GridFunction v(vfes);
   //v = 0.0;
   ///GridFunction p(pfes);
   //p = 0.0;

   // boundary conditions have to be set!
   //v.ProjectCoefficient(velocity);
   // Set the Dirichlet values in the solution vector
   //v.ProjectBdrCoefficient(vac_vdbc, vdbc_bdr);
   ConstantCoefficient vdbcCoef(vdbc_val1);
   v.ProjectCoefficient(vdbcCoef, vdbc_bdr);
   ConstantCoefficient pdbcCoef(pdbc_val1);
   p.ProjectCoefficient(pdbcCoef, pdbc_bdr);

   // 9. SOLVER
   // setup MINRES solver
   int maxIter(30000);
   double rtol(1.e-10);
   double atol(1.e-10);

   StopWatch chrono;
   chrono.Clear();
   chrono.Start();
   MINRESSolver solver;
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

   
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream v_sock(vishost, visport);
      v_sock.precision(8);
      v_sock << "solution\n" << mesh << v << "window_title 'Velocity'" << std::endl;
      socketstream p_sock(vishost, visport);
      p_sock.precision(8);
      p_sock << "solution\n" << mesh << p << "window_title 'Pressure'" << std::endl;
   }



   // 8. Set up the linear form b(.) which corresponds to the right-hand side of
   //    the FEM linear system.
   
   //b->Assemble();


   // 9. Set up the bilinear form a(.,.) on the finite element space
   //    corresponding to the left hand side.
   //BilinearForm *a = new BilinearForm(fespace);

   //std::cout << "using DiffusionIntegrator."<< std::endl;
   
   //a->Assemble();

   //SparseMatrix A;
   //Vector B, X;
   //a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);
   //std::cout << "Size of linear system: " << A.Height() << std::endl;

//#ifndef MFEM_USE_SUITESPARSE
   // 10. Define a simple Jacobi preconditioner and use it to
   //     solve the system A X = B with PCG.
//   GSSmoother M(A);
//   PCG(A, M, B, X, 1, 200, 1e-12, 0.0);
//#else
   // 10. If MFEM was compiled with SuiteSparse, use UMFPACK to solve the system.
//   UMFPackSolver umf_solver;
//   umf_solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
//   umf_solver.SetOperator(A);
//   umf_solver.Mult(B, X);
//#endif

   // 11. Recover the solution as a finite element grid function.
//   a->RecoverFEMSolution(X, *b, x);

   // 12. Save the refined mesh and the solution. This output can be viewed later
   //     using GLVis: "glvis -m refined.mesh -g sol.gf".
//   ofstream mesh_ofs("refined.mesh");
//   mesh_ofs.precision(8);
//   mesh->Print(mesh_ofs);
//   ofstream sol_ofs("sol.gf");
//   sol_ofs.precision(8);
//   x.Save(sol_ofs);

   // 13. Send the solution by socket to a GLVis server.
//   if (visualization)
//   {
//      char vishost[] = "localhost";
//      int  visport   = 19916;
//      socketstream sol_sock(vishost, visport);
//      sol_sock.precision(8);
//      sol_sock << "solution\n" << *mesh << x << flush;
//   }

   // 14. Save data in the VisIt format
//   VisItDataCollection visit_dc("Example1", mesh);
//   visit_dc.RegisterField("solution", &x);
//   visit_dc.Save();

   // 15. Free the used memory.
//   delete a;
//  delete b;
// delete fespace;
//   if (own_fec) { delete fec; }
//   delete mesh;

   return 0;
}