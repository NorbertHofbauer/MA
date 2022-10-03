//                          MFEM Example 1 - NURBS Version
//
// Compile with: make nurbs_ex1
//
// Sample runs:  nurbs_ex1 -m ../../data/square-nurbs.mesh -o 2 -no-ibp
//               nurbs_ex1 -m ../../data/square-nurbs.mesh -o 2 --weak-bc
//               nurbs_ex1 -m ../../data/cube-nurbs.mesh -o 2 -no-ibp
//               nurbs_ex1 -m ../../data/pipe-nurbs-2d.mesh -o 2 -no-ibp
//               nurbs_ex1 -m ../../data/square-disc-nurbs.mesh -o -1
//               nurbs_ex1 -m ../../data/disc-nurbs.mesh -o -1
//               nurbs_ex1 -m ../../data/pipe-nurbs.mesh -o -1
//               nurbs_ex1 -m ../../data/beam-hex-nurbs.mesh -pm 1 -ps 2
//
// Description:  This example code demonstrates the use of MFEM to define a
//               simple finite element discretization of the Laplace problem
//               -Delta u = 1 with homogeneous Dirichlet boundary conditions.
//               The boundary conditions can be enforced either strongly or weakly.
//               Specifically, we discretize using a FE space of the specified
//               order, or if order < 1 using an isoparametric/isogeometric
//               space (i.e. quadratic for quadratic curvilinear mesh, NURBS for
//               NURBS mesh, etc.)
//
//               The example highlights the use of mesh refinement, finite
//               element grid functions, as well as linear and bilinear forms
//               corresponding to the left-hand side and right-hand side of the
//               discrete linear system. We also cover the explicit elimination
//               of essential boundary conditions, static condensation, and the
//               optional connection to the GLVis tool for visualization.

#include "mfem.hpp"
#include <fstream>
#include <iostream>

#include <algorithm>
#include <cmath>

using namespace std; // bad idea!!!!
using namespace mfem;

double lambda(const Vector &x);

/** Class for integrating the bilinear form a(u,v) := (Q Laplace u, v) where Q
    can be a scalar coefficient. */
class Diffusion2Integrator: public BilinearFormIntegrator
{
private:
#ifndef MFEM_THREAD_SAFE
   Vector shape,laplace;
#endif
   Coefficient *Q;

public:
   /// Construct a diffusion integrator with coefficient Q = 1
   Diffusion2Integrator() { Q = NULL; }

   /// Construct a diffusion integrator with a scalar coefficient q
   Diffusion2Integrator (Coefficient &q) : Q(&q) { }

   /** Given a particular Finite Element
       computes the element stiffness matrix elmat. */
   virtual void AssembleElementMatrix(const FiniteElement &el,
                                      ElementTransformation &Trans,
                                      DenseMatrix &elmat)
   {
      int nd = el.GetDof();
      int dim = el.GetDim();
      double w;

#ifdef MFEM_THREAD_SAFE
      Vector shape(nd);
      Vector laplace(nd);
#else
      shape.SetSize(nd);
      laplace.SetSize(nd);
#endif
      elmat.SetSize(nd);

      const IntegrationRule *ir = IntRule;
      if (ir == NULL)
      {
         int order;
         if (el.Space() == FunctionSpace::Pk)
         {
            order = 2*el.GetOrder() - 2;
         }
         else
         {
            order = 2*el.GetOrder() + dim - 1;
         }

         if (el.Space() == FunctionSpace::rQk)
         {
            ir = &RefinedIntRules.Get(el.GetGeomType(),order);
         }
         else
         {
            ir = &IntRules.Get(el.GetGeomType(),order);
         }
      }

      elmat = 0.0;
      for (int i = 0; i < ir->GetNPoints(); i++)
      {
         const IntegrationPoint &ip = ir->IntPoint(i);
         Trans.SetIntPoint(&ip);
         w = -ip.weight * Trans.Weight();

         el.CalcShape(ip, shape);
         el.CalcPhysLaplacian(Trans, laplace);

         if (Q)
         {
            w *= Q->Eval(Trans, ip);
         }

         for (int jj = 0; jj < nd; jj++)
         {
            for (int ii = 0; ii < nd; ii++)
            {
               elmat(ii, jj) += w*shape(ii)*laplace(jj);
            }
         }
      }
   }

};

int main(int argc, char *argv[])
{
   // 1. Parse command-line options.
   //const char *mesh_file = "../../data/star.mesh";
   const char *mesh_file = "../../../MA/mesh/pipe-nurbs-boundary-test_2.mesh";
   const char *per_file  = "none";
   int ref_levels = 0;
   Array<int> master(0);
   Array<int> slave(0);
   bool static_cond = false;
   bool visualization = 1;
   bool ibp = 1;
   bool strongBC = 1;
   double kappa = -1;
   Array<int> order(1);
   order[0] = 1;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&ref_levels, "-r", "--refine",
                  "Number of times to refine the mesh uniformly, -1 for auto.");
   args.AddOption(&per_file, "-p", "--per",
                  "Periodic BCS file.");
   args.AddOption(&master, "-pm", "--master",
                  "Master boundaries for periodic BCs");
   args.AddOption(&slave, "-ps", "--slave",
                  "Slave boundaries for periodic BCs");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree) or -1 for"
                  " isoparametric space.");
   args.AddOption(&ibp, "-ibp", "--ibp", "-no-ibp",
                  "--no-ibp",
                  "Selects the standard weak form (IBP) or the nonstandard (NO-IBP).");
   args.AddOption(&strongBC, "-sbc", "--strong-bc", "-wbc",
                  "--weak-bc",
                  "Selects strong or weak enforcement of Dirichlet BCs.");
   args.AddOption(&kappa, "-k", "--kappa",
                  "Sets the SIPG penalty parameters, should be positive."
                  " Negative values are replaced with (order+1)^2.");
   args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc",
                  "--no-static-condensation", "Enable static condensation.");
   args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                  "--no-visualization",
                  "Enable or disable GLVis visualization.");
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   if (!strongBC & (kappa < 0))
   {
      kappa = 4*(order.Max()+1)*(order.Max()+1);
   }
   args.PrintOptions(cout);

   // 2. Read the mesh from the given mesh file. We can handle triangular,
   //    quadrilateral, tetrahedral, hexahedral, surface and volume meshes with
   //    the same code.
   Mesh *mesh = new Mesh(mesh_file, 1, 1);
   int dim = mesh->Dimension();

   // 3. Refine the mesh to increase the resolution. In this example we do
   //    'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
   //    largest number that gives a final mesh with no more than 50,000
   //    elements.
   //# higher ref_levels means more refinement and more resolution
   {
      if (ref_levels < 0)
      {
         ref_levels =
            (int)floor(log(5000./mesh->GetNE())/log(2.)/dim);
      }

      for (int l = 0; l < ref_levels; l++)
      {
         mesh->UniformRefinement();
      }
      mesh->PrintInfo();
   }

   // 4. Define a finite element space on the mesh. Here we use continuous
   //    Lagrange finite elements of the specified order. If order < 1, we
   //    instead use an isoparametric/isogeometric space.


   FiniteElementCollection *fec;
   NURBSExtension *NURBSext = NULL;
   int own_fec = 0;

   if (mesh->NURBSext)
   {
      cout << "Using Nurbs mesh." << endl;
      cout << "given Order" << order[0] << endl;
      fec = new NURBSFECollection(order[0]);
      own_fec = 1;

      int nkv = mesh->NURBSext->GetNKV();
      if (order.Size() == 1)
      {
         int tmp = order[0];
         order.SetSize(nkv);
         order = tmp;
      }

      if (order.Size() != nkv ) { mfem_error("Wrong number of orders set."); }
      NURBSext = new NURBSExtension(mesh->NURBSext, order);

   }
   else if (order[0] == -1) // Isoparametric
   {
    //# ignore for nurbs mesh
     cout << "isoparametric?" << endl;
      if (mesh->GetNodes())
      {
         fec = mesh->GetNodes()->OwnFEC();
         own_fec = 0;
         cout << "Using isoparametric FEs: " << fec->Name() << endl;
      }
      else
      {
         cout <<"Mesh does not have FEs --> Assume order 1.\n";
         fec = new H1_FECollection(1, dim);
         own_fec = 1;
      }
   }
   else
   {
      //# ignore for nurbs mesh
      if (order.Size() > 1) { cout <<"Wrong number of orders set, needs one.\n"; }
      fec = new H1_FECollection(abs(order[0]), dim);
      own_fec = 1;
   }

   FiniteElementSpace *fespace = new FiniteElementSpace(mesh, NURBSext, fec);
   cout << "Number of finite element unknowns: "
        << fespace->GetTrueVSize() << endl;

  //# check if no use of integration by parts is possible for the current code
   if (!ibp)
   {
      if (!mesh->NURBSext)
      {
         cout << "No integration by parts requires a NURBS mesh."<< endl;
         return 2;
      }
      if (mesh->NURBSext->GetNP()>1)
      {
         cout << "No integration by parts requires a NURBS mesh, with only 1 patch."<<
              endl;
         cout << "A C_1 discretisation is required."<< endl;
         cout << "Currently only C_0 multipatch coupling implemented."<< endl;
         return 3;
      }
      if (order[0]<2)
      {
         cout << "No integration by parts requires at least quadratic NURBS."<< endl;
         cout << "A C_1 discretisation is required."<< endl;
         return 4;
      }
   }

   // 5. Determine the list of true (i.e. conforming) essential boundary dofs.
   //    In this example, the boundary conditions are defined by marking all
   //    the boundary attributes from the mesh as essential (Dirichlet) and
   //    converting them to a list of true dofs.

   // first set markers for the boundary attributes

   Array<int> dbc_bdr(mesh->bdr_attributes.Max());
   Array<int> nbc1_bdr(mesh->bdr_attributes.Max());
   Array<int> nbc2_bdr(mesh->bdr_attributes.Max());
   Array<int> rbc_bdr(mesh->bdr_attributes.Max());

   // we assume that the boundary attribute 1,2 are dirchlet and 3,4 are neumann
   dbc_bdr = 0; dbc_bdr[0] = 1; dbc_bdr[1] = 1;  //dirichlet
   //dbc_bdr = 0; dbc_bdr[1] = 1;  //dirichlet
   nbc1_bdr = 0; nbc1_bdr[2] = 1;  //neumann1
   nbc2_bdr = 0; nbc2_bdr[3] = 1;  //neumann2
   // we will lay a robin bc on top of the dirchlet...totally nonsense ^^
   rbc_bdr = 0; rbc_bdr[0] = 1; //robin

   // For a continuous basis the linear system must be modified to enforce an
   // essential (Dirichlet) boundary condition.
   Array<int> ess_tdof_list(0);
   fespace->GetEssentialTrueDofs(dbc_bdr, ess_tdof_list);

   // 5. Setup the various coefficients needed for the Laplace operator and the
   //    various boundary conditions. In general these coefficients could be
   //    functions of position but here we use only constants.

   // values should be defined at top of class main, but for now it doesn't really matter
   double mat_val = 1.0;
   double dbc1_val = -999.0;
   double dbc2_val = 111.0;
   double nbc1_val = 0;
   double nbc2_val = 0;
   double rbc_a_val = 0;
   double rbc_b_val = 0;

   ConstantCoefficient matCoef(mat_val);
   ConstantCoefficient nbc1Coef(nbc1_val);
   ConstantCoefficient nbc2Coef(nbc2_val);
   ConstantCoefficient rbcACoef(rbc_a_val);
   ConstantCoefficient rbcBCoef(rbc_b_val);


   // okay lets try lambda functions

   auto lambda = [] (const Vector &first)->double
   {
      cout << " x(0) = " << first(0) << " x(1) = " << first(1) << " x(2) = " << first(2) << endl;
      return first(0);
   };
   
   FunctionCoefficient lambdaCoef(lambda);

   // Since the n.Grad(u) terms arise by integrating -Div(m Grad(u)) by parts we
   // must introduce the coefficient 'm' into the boundary conditions.
   // Therefore, in the case of the Neumann BC, we actually enforce m n.Grad(u)
   // = m g rather than simply n.Grad(u) = g.
   ProductCoefficient m_nbc1Coef(matCoef, nbc1Coef);
   ProductCoefficient m_nbc2Coef(matCoef, nbc2Coef);
   ProductCoefficient m_rbcACoef(matCoef, rbcACoef);
   ProductCoefficient m_rbcBCoef(matCoef, rbcBCoef);

   // construct array for Coefficients to combine the boundary values
   // dirichlet
   VectorArrayCoefficient vac_dbc(1);
   Vector vec_dbc(mesh->bdr_attributes.Max());
   vec_dbc = 0.0;
   vec_dbc(0) = dbc1_val;
   vec_dbc(1) = dbc2_val;
   vac_dbc.Set(0,new PWConstCoefficient(vec_dbc));
   
   //cout << " vac_dbc = " << vac_dbc.GetCoeff(0) << endl;

   // 7. Define the solution vector x as a finite element grid function
   //    corresponding to fespace. Initialize x with initial guess of zero,
   //    which satisfies the boundary conditions.
   GridFunction x(fespace);
   x = 0.0;

   // 6. Set up the linear form b(.) which corresponds to the right-hand side of
   //    the FEM linear system, which in this case is (1,phi_i) where phi_i are
   //    the basis functions in the finite element fespace.
   //ConstantCoefficient one(1.0);
   //ConstantCoefficient zero(0.0);

   LinearForm *b = new LinearForm(fespace);

   // not needed because the system is -laplace(u)=0
   //b->AddDomainIntegrator(new DomainLFIntegrator(one));

   // Set the Dirichlet values in the solution vector
   x.ProjectBdrCoefficient(vac_dbc, dbc_bdr);
   x.ProjectBdrCoefficient(lambdaCoef, rbc_bdr);
   //cout << " x = " << x << endl;

   // Add the desired value for n.Grad(u) on the Neumann boundary
   b->AddBoundaryIntegrator(new BoundaryLFIntegrator(m_nbc1Coef), nbc1_bdr);
   //b->AddBoundaryIntegrator(new BoundaryLFIntegrator(m_nbc2Coef), nbc2_bdr);
   b->AddBoundaryIntegrator(new BoundaryLFIntegrator(lambdaCoef), nbc2_bdr);

   // robin
   b->AddBoundaryIntegrator(new BoundaryLFIntegrator(m_rbcBCoef), rbc_bdr);


   b->Assemble();


   // 8. Set up the bilinear form a(.,.) on the finite element space
   //    corresponding to the Laplacian operator -Delta, by adding the Diffusion
   //    domain integrator.
   BilinearForm *a = new BilinearForm(fespace);

   cout << "using DiffusionIntegrator."<< endl;
   //a->AddDomainIntegrator(new DiffusionIntegrator(one));
   a->AddDomainIntegrator(new DiffusionIntegrator(matCoef));

   // Add a Mass Integrator on the Robin Boundary
   a->AddDomainIntegrator(new MassIntegrator(m_rbcACoef),rbc_bdr);

   // 9. Assemble the bilinear form and the corresponding linear system,
   //    applying any necessary transformations such as: eliminating boundary
   //    conditions, applying conforming constraints for non-conforming AMR,
   //    static condensation, etc.
   if (static_cond) { a->EnableStaticCondensation(); }
   a->Assemble();

   SparseMatrix A;
   Vector B, X;
   a->FormLinearSystem(ess_tdof_list, x, *b, A, X, B);

   cout << "Size of linear system: " << A.Height() << endl;

#ifndef MFEM_USE_SUITESPARSE
   // 10. Define a simple Jacobi preconditioner and use it to
   //     solve the system A X = B with PCG.
   GSSmoother M(A);
   PCG(A, M, B, X, 1, 200, 1e-12, 0.0);
#else
   // 10. If MFEM was compiled with SuiteSparse, use UMFPACK to solve the system.
   UMFPackSolver umf_solver;
   umf_solver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_METIS;
   umf_solver.SetOperator(A);
   umf_solver.Mult(B, X);
#endif

   // 11. Recover the solution as a finite element grid function.
   a->RecoverFEMSolution(X, *b, x);

   // 12. Save the refined mesh and the solution. This output can be viewed later
   //     using GLVis: "glvis -m refined.mesh -g sol.gf".
   ofstream mesh_ofs("refined.mesh");
   mesh_ofs.precision(8);
   mesh->Print(mesh_ofs);
   ofstream sol_ofs("sol.gf");
   sol_ofs.precision(8);
   x.Save(sol_ofs);

   // 13. Send the solution by socket to a GLVis server.
   if (visualization)
   {
      char vishost[] = "localhost";
      int  visport   = 19916;
      socketstream sol_sock(vishost, visport);
      sol_sock.precision(8);
      sol_sock << "solution\n" << *mesh << x << flush;
   }

   // 14. Save data in the VisIt format
   VisItDataCollection visit_dc("Example1", mesh);
   visit_dc.RegisterField("solution", &x);
   visit_dc.Save();

   // 15. Free the used memory.
   delete a;
   delete b;
   delete fespace;
   if (own_fec) { delete fec; }
   delete mesh;

   return 0;
}


//double lambda(const Vector &x)
//{  
   //double n = (x(0)*x(0)+x(1)*x(1)+x(2)*x(2))/100;
//   double n = x(0);
//   cout << " x(0) = " << x(0) << " x(1) = " << x(1) << " x(2) = " << x(2) << endl;
//   return n;
//}