#pragma once

#include "mfem.hpp"

class MixedMassIntegrator : public mfem::BilinearFormIntegrator
{

private:
   int vdim;
   mfem::Vector shape, te_shape, vec;
   mfem::DenseMatrix partelmat;
   int Q_order;

protected:
   mfem::Coefficient *Q;
   // PA extension
   mfem::Vector pa_data;
   const mfem::DofToQuad *trial_maps, *test_maps; ///< Not owned
   const mfem::GeometricFactors *geom;  ///< Not owned
   int dim, ne, nq;
   int trial_dofs1D, test_dofs1D, quad1D;
public:
   /// Construct an integrator with coefficient 1.0
   MixedMassIntegrator()
      : vdim(-1), Q_order(0), Q(NULL), trial_maps{NULL}, test_maps{NULL}, geom{NULL} { }
   /** Construct an integrator with scalar coefficient q.  If possible, save
diagonal block    memory by using a scalar integrator since the resulting matrix is block
       diagonal with the same  repeated. */
   MixedMassIntegrator(mfem::Coefficient &q, int qo = 0)
      : vdim(-1), Q_order(qo), Q(&q), trial_maps{NULL}, test_maps{NULL}, geom{NULL} { }
   MixedMassIntegrator(mfem::Coefficient &q, const mfem::IntegrationRule *ir)
      : mfem::BilinearFormIntegrator(ir), vdim(-1), Q_order(0), Q(&q), trial_maps{NULL}, test_maps{NULL}, geom{NULL} { }
   
   virtual void AssembleElementMatrix2(const mfem::FiniteElement &trial_fe,
                                       const mfem::FiniteElement &test_fe,
                                       mfem::ElementTransformation &Trans,
                                       mfem::DenseMatrix &elmat);

   /*
   virtual void AssemblePA(const mfem::FiniteElementSpace &trial_fes,
                           const mfem::FiniteElementSpace &test_fes);

   virtual void AddMultPA(const mfem::Vector &x, mfem::Vector &y) const;
   virtual void AddMultTransposePA(const mfem::Vector &x, mfem::Vector &y) const;
   */
   using BilinearFormIntegrator::AssemblePA;
   
   static const mfem::IntegrationRule &GetRule(const mfem::FiniteElement &trial_fe,
                                         const mfem::FiniteElement &test_fe,
                                         mfem::ElementTransformation &Trans);
};

void MixedMassIntegrator::AssembleElementMatrix2(
   const mfem::FiniteElement &trial_fe, const mfem::FiniteElement &test_fe,
   mfem::ElementTransformation &Trans, mfem::DenseMatrix &elmat)
{
   int tr_nd = trial_fe.GetDof();
   int te_nd = test_fe.GetDof();
   int dim = test_fe.GetDim();
   double norm;

   // If vdim is not set, set it to the space dimension
   vdim = (vdim == -1) ? Trans.GetSpaceDim() : vdim;
   
   // If vdim is not set, set it to the space dimension
   //vdim = (vdim == -1) ? test_fe.GetDim() : vdim;
   
   //elmat.SetSize(te_nd*vdim, tr_nd*vdim);
   elmat.SetSize(te_nd*vdim, tr_nd);
   shape.SetSize(tr_nd);
   te_shape.SetSize(te_nd);
   partelmat.SetSize(te_nd, tr_nd);

   const mfem::IntegrationRule *ir = IntRule;
   if (ir == NULL)
   {
      int order = (trial_fe.GetOrder() + test_fe.GetOrder() +
                   Trans.OrderW() + Q_order);

      if (trial_fe.Space() == mfem::FunctionSpace::rQk)
      {
         ir = &mfem::RefinedIntRules.Get(trial_fe.GetGeomType(), order);
      }
      else
      {
         ir = &mfem::IntRules.Get(trial_fe.GetGeomType(), order);
      }
   }

   elmat = 0.0;
   for (int s = 0; s < ir->GetNPoints(); s++)
   {
      const mfem::IntegrationPoint &ip = ir->IntPoint(s);
      trial_fe.CalcShape(ip, shape);
      test_fe.CalcShape(ip, te_shape);

      Trans.SetIntPoint(&ip);
      norm = ip.weight * Trans.Weight();

      MultVWt(te_shape, shape, partelmat);

      if (Q)
      {
         norm *= Q->Eval(Trans, ip);
      }
      partelmat *= norm;
      for (int k = 0; k < vdim; k++)
      {
         //elmat.AddMatrix(partelmat, te_nd*k, tr_nd*k);
         elmat.AddMatrix(partelmat, te_nd*k, 0);
      }
   }
}

const mfem::IntegrationRule &MixedMassIntegrator::GetRule(const mfem::FiniteElement
                                                   &trial_fe,
                                                   const mfem::FiniteElement &test_fe,
                                                   mfem::ElementTransformation &Trans)
{
   int order = Trans.OrderGrad(&trial_fe) + test_fe.GetOrder() + Trans.OrderJ();
   return mfem::IntRules.Get(trial_fe.GetGeomType(), order);
}
