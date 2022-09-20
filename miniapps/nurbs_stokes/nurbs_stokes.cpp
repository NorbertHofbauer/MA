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

void vec_up(const Vector &, Vector &);

int main(int argc, char *argv[])
{
   // 0. Setup
   double vdbc_val_0_0 = 0;
   double vdbc_val_0_1 = 0;
   double vdbc_val_2_0 = 0;
   double vdbc_val_2_1 = 5;
   double pdbc_val = 55;
   //double pdbc_val2 = 22;
   //double vnbc_val = 0;
   //double pnbc_val = 0;

   /*
   auto lambda_up = [vdbc_val_2_0,vdbc_val_2_1] (const Vector &bdr_up)->std::vector<double>
   {
      std::vector<double> vec_up(3);
      vec_up(0)=vdbc_val_2_0;
      vec_up(1)=vdbc_val_2_1;
      vec_up(2)=0.0;
      std::cout << " x(0) = " << bdr_up(0) << " x(1) = " << bdr_up(1) << " x(2) = " << bdr_up(2) << std::endl;
      return vec_up;
   };
   */

   // 1. Parse command-line options.
   //const char *mesh_file = "../../../MA/mesh/pipe-nurbs-boundary-test_2.mesh";
   const char *mesh_file = "../../../MA/mesh/quad_nurbs.mesh";

   int ref_levels = 2;
   bool visualization = 1;
   Array<int> order(1);
   order[0] = 0;

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

   NURBSExtension *vNURBSext = NULL;
   NURBSExtension *pNURBSext = NULL;
   
   std::cout << "Using Nurbs mesh." << std::endl;
   std::cout << "mesh Dimension " << sdim << std::endl;
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

   FiniteElementSpace *vfes = new FiniteElementSpace(mesh, vNURBSext, vfec, sdim);
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
   std::cout << " bdr attributes " << mesh->bdr_attributes.Max() << std::endl;
   // Neumann
   //Array<int> vnbc_bdr(mesh->bdr_attributes.Max());
   //Array<int> pnbc_bdr(mesh->bdr_attributes.Max());
   
   // we assume that the boundary attribute 1,2 are dirchlet, we want to use either pressure or velocity or a mix
   //vdbc_bdr = 0; vdbc_bdr[0] = 1; vdbc_bdr[2] = 1;
   vdbc_bdr = 0; vdbc_bdr[2] = 1; 
   pdbc_bdr = 0; pdbc_bdr[1] = 1;

   Array<int> vel_ess_tdof_list;
   Array<int> pres_ess_tdof_list;
   //vfes->GetEssentialTrueDofs(vdbc_bdr, vel_ess_tdof_list);
   vfes->GetEssentialTrueDofs(vdbc_bdr, vel_ess_tdof_list);
   pfes->GetEssentialTrueDofs(pdbc_bdr, pres_ess_tdof_list);   

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
   Vector vzero(sdim);
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
   a.AddDomainIntegrator(new VectorDiffusionIntegrator(sdim));
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
   b.FormRectangularSystemMatrix(pres_ess_tdof_list, vel_ess_tdof_list, B);
   
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
   
   std::cout << " v = " << v << std::endl;

   auto lambda_up = [&sdim](const Vector &ControlPointIn, Vector &ControlPointOut) -> void
   {
      ControlPointOut[0] = 777777;
      ControlPointOut[1] = 999999;
      
      std::cout << " v(0) = " << ControlPointIn(0) << " v(1) = " << ControlPointIn(1) << std::endl;
      std::cout << " x(0) = " << ControlPointOut(0) << " x(1) = " << ControlPointOut(1) << std::endl;      
      return;
   };

   //VectorFunctionCoefficient vfc_up(sdim, lambda_up);

   //Array<int> pwattr(1);
   //pwattr[1] = 1;
   //Array<VectorCoefficient> pwvcar(1);
   //pwvcar[0] = vfc_up;

   //PWVectorCoefficient pwvc(sdim,pwattr,pwvcar);

   //std::cout << " jkl = " << jkl << std::endl;

   VectorFunctionCoefficient vfc_up(sdim, lambda_up);
   v.ProjectBdrCoefficient(vfc_up,vdbc_bdr);
   //v.ProjectCoefficient(pwvc);
   std::cout << " v = " << v << std::endl;


   //ConstantCoefficient vdbcCoef(vdbc_val_2_1);
   //v.ProjectCoefficient(vdbcCoef, vdbc_bdr);
   std::cout << " p = " << p << std::endl;
   ConstantCoefficient pdbcCoef(pdbc_val);
   p.ProjectBdrCoefficient(pdbcCoef, pdbc_bdr);
   std::cout << " p = " << p << std::endl;


   std::cout << " x.size = " << " 0 to " << x.BlockSize(0)-1 << std::endl; 
   for (int i = 0; i < x.BlockSize(0); i++) {
      std::cout << x(i) << std::endl;
   }

   std::cout << " x.size = " << x.BlockSize(0) << " to " << x.BlockSize(0)+x.BlockSize(1)-1 << std::endl; 
   for (int i = x.BlockSize(0); i < x.BlockSize(0)+x.BlockSize(1); i++) {
      std::cout << x(i) << std::endl;
   }

   int getArrayLength = sizeof(vel_ess_tdof_list) / sizeof(int);
   std::cout << "sizeof(vel_ess_tdof_list) " << getArrayLength << std::endl;
   for (int i = 0; i < getArrayLength; i++) {
      std::cout << vel_ess_tdof_list[i] << std::endl;
   }

   getArrayLength = sizeof(pres_ess_tdof_list) / sizeof(int);
   std::cout << "sizeof(pres_ess_tdof_list) " << getArrayLength << std::endl;
   for (int i = 0; i < getArrayLength; i++) {
      std::cout << pres_ess_tdof_list[i] << std::endl;
   }

   // 9. SOLVER
   // setup MINRES solver
   int maxIter(300);
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
   
   std::cout << " v.sol = " << v << std::endl;
   std::cout << " p.sol = " << p << std::endl;
   std::cout << " vfes = " << vfes->GetTrueVSize() << std::endl;
   std::cout << " pfes = " << pfes->GetTrueVSize() << std::endl;


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

   return 0;
}