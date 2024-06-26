
#pragma once

#include "mfem.hpp"

// Vector DG diffusion integrator
class VectorDGDiffusionIntegrator : public mfem::BilinearFormIntegrator
{
protected:
   mfem::Coefficient *Q;
   mfem::MatrixCoefficient *MQ;
   double sigma, kappa;
   int vdim;

   // these are not thread-safe!
   mfem::Vector shape1, shape2, dshape1dn, dshape2dn, nor, nh, ni;
   mfem::DenseMatrix jmat, dshape1, dshape2, mq, adjJ;

public:
   VectorDGDiffusionIntegrator(const double s, const double k,
                               const int vdim_) : Q(NULL), MQ(NULL), sigma(s), kappa(k), vdim(vdim_) { }
   VectorDGDiffusionIntegrator(mfem::Coefficient &q, const double s, const double k,
                               const int vdim_) : Q(&q), MQ(NULL), sigma(s), kappa(k), vdim(vdim_) { }
   VectorDGDiffusionIntegrator(mfem::MatrixCoefficient &q, const double s,
                               const double k, const int vdim_) : Q(NULL), MQ(&q), sigma(s), kappa(k),
      vdim(vdim_) { }

   using BilinearFormIntegrator::AssembleFaceMatrix;
   virtual void AssembleFaceMatrix( const mfem::FiniteElement &el1,
                                    const mfem::FiniteElement &el2,
                                    mfem::FaceElementTransformations &Trans,
                                    mfem::DenseMatrix &full_elmat);
};

void VectorDGDiffusionIntegrator::AssembleFaceMatrix(
   const mfem::FiniteElement &el1, const mfem::FiniteElement &el2,
   mfem::FaceElementTransformations &Trans, mfem::DenseMatrix &full_elmat)
{
   int dim, ndof1, ndof2, ndofs;
   bool kappa_is_nonzero = (kappa != 0.);
   double w, wq = 0.0;

   dim = el1.GetDim();
   ndof1 = el1.GetDof();

   nor.SetSize(dim);
   nh.SetSize(dim);
   ni.SetSize(dim);
   adjJ.SetSize(dim);
   if (MQ)
   {
      mq.SetSize(dim);
   }

   shape1.SetSize(ndof1);
   dshape1.SetSize(ndof1, dim);
   dshape1dn.SetSize(ndof1);
   if (Trans.Elem2No >= 0)
   {
      ndof2 = el2.GetDof();
      shape2.SetSize(ndof2);
      dshape2.SetSize(ndof2, dim);
      dshape2dn.SetSize(ndof2);
   }
   else
   {
      ndof2 = 0;
   }

   ndofs = ndof1 + ndof2;
   mfem::DenseMatrix elmat;
   elmat.SetSize(ndofs);
   elmat = 0.0;
   if (kappa_is_nonzero)
   {
      jmat.SetSize(ndofs);
      jmat = 0.;
   }

   const mfem::IntegrationRule *ir = IntRule;
   if (ir == NULL)
   {
      // a simple choice for the integration order; is this OK?
      int order;
      if (ndof2)
      {
         order = 2 * std::max(el1.GetOrder(), el2.GetOrder());
      }
      else
      {
         order = 2 * el1.GetOrder();
      }
      ir = &mfem::IntRules.Get(Trans.GetGeometryType(), order);
   }

   // assemble: < {(Q \nabla u).n},[v] >      --> elmat
   //           kappa < {h^{-1} Q} [u],[v] >  --> jmat
   for (int p = 0; p < ir->GetNPoints(); p++)
   {
      const mfem::IntegrationPoint &ip = ir->IntPoint(p);

      // Set the integration point in the face and the neighboring elements
      Trans.SetAllIntPoints(&ip);

      // Access the neighboring elements' integration points
      // Note: eip2 will only contain valid data if Elem2 exists
      const mfem::IntegrationPoint &eip1 = Trans.GetElement1IntPoint();
      const mfem::IntegrationPoint &eip2 = Trans.GetElement2IntPoint();

      if (dim == 1)
      {
         nor(0) = 2 * eip1.x - 1.0;
      }
      else
      {
         CalcOrtho(Trans.Jacobian(), nor);
      }

      el1.CalcShape(eip1, shape1);
      el1.CalcDShape(eip1, dshape1);
      w = ip.weight / Trans.Elem1->Weight();
      if (ndof2)
      {
         w /= 2;
      }
      if (!MQ)
      {
         if (Q)
         {
            w *= Q->Eval(*Trans.Elem1, eip1);
         }
         ni.Set(w, nor);
      }
      else
      {
         nh.Set(w, nor);
         MQ->Eval(mq, *Trans.Elem1, eip1);
         mq.MultTranspose(nh, ni);
      }
      CalcAdjugate(Trans.Elem1->Jacobian(), adjJ);
      adjJ.Mult(ni, nh);
      if (kappa_is_nonzero)
      {
         wq = ni * nor;
      }
      // Note: in the jump term, we use 1/h1 = |nor|/det(J1) which is
      // independent of Loc1 and always gives the size of element 1 in
      // direction perpendicular to the face. Indeed, for linear transformation
      //     |nor|=measure(face)/measure(ref. face),
      //   det(J1)=measure(element)/measure(ref. element),
      // and the ratios measure(ref. element)/measure(ref. face) are
      // compatible for all element/face pairs.
      // For example: meas(ref. tetrahedron)/meas(ref. triangle) = 1/3, and
      // for any tetrahedron vol(tet)=(1/3)*height*area(base).
      // For interior faces: q_e/h_e=(q1/h1+q2/h2)/2.

      dshape1.Mult(nh, dshape1dn);
      for (int i = 0; i < ndof1; i++)
         for (int j = 0; j < ndof1; j++)
         {
            elmat(i, j) += shape1(i) * dshape1dn(j);
         }

      if (ndof2)
      {
         el2.CalcShape(eip2, shape2);
         el2.CalcDShape(eip2, dshape2);
         w = ip.weight / 2 / Trans.Elem2->Weight();
         if (!MQ)
         {
            if (Q)
            {
               w *= Q->Eval(*Trans.Elem2, eip2);
            }
            ni.Set(w, nor);
         }
         else
         {
            nh.Set(w, nor);
            MQ->Eval(mq, *Trans.Elem2, eip2);
            mq.MultTranspose(nh, ni);
         }
         CalcAdjugate(Trans.Elem2->Jacobian(), adjJ);
         adjJ.Mult(ni, nh);
         if (kappa_is_nonzero)
         {
            wq += ni * nor;
         }

         dshape2.Mult(nh, dshape2dn);

         for (int i = 0; i < ndof1; i++)
            for (int j = 0; j < ndof2; j++)
            {
               elmat(i, ndof1 + j) += shape1(i) * dshape2dn(j);
            }

         for (int i = 0; i < ndof2; i++)
            for (int j = 0; j < ndof1; j++)
            {
               elmat(ndof1 + i, j) -= shape2(i) * dshape1dn(j);
            }

         for (int i = 0; i < ndof2; i++)
            for (int j = 0; j < ndof2; j++)
            {
               elmat(ndof1 + i, ndof1 + j) -= shape2(i) * dshape2dn(j);
            }
      }

      if (kappa_is_nonzero)
      {
         // only assemble the lower triangular part of jmat
         wq *= kappa;
         for (int i = 0; i < ndof1; i++)
         {
            const double wsi = wq * shape1(i);
            for (int j = 0; j <= i; j++)
            {
               jmat(i, j) += wsi * shape1(j);
            }
         }
         if (ndof2)
         {
            for (int i = 0; i < ndof2; i++)
            {
               const int i2 = ndof1 + i;
               const double wsi = wq * shape2(i);
               for (int j = 0; j < ndof1; j++)
               {
                  jmat(i2, j) -= wsi * shape1(j);
               }
               for (int j = 0; j <= i; j++)
               {
                  jmat(i2, ndof1 + j) += wsi * shape2(j);
               }
            }
         }
      }
   }

   // elmat := -elmat + sigma*elmat^t + jmat
   if (kappa_is_nonzero)
   {
      for (int i = 0; i < ndofs; i++)
      {
         for (int j = 0; j < i; j++)
         {
            double aij = elmat(i, j), aji = elmat(j, i), mij = jmat(i, j);
            elmat(i, j) = sigma * aji - aij + mij;
            elmat(j, i) = sigma * aij - aji + mij;
         }
         elmat(i, i) = (sigma - 1.) * elmat(i, i) + jmat(i, i);
      }
   }
   else
   {
      for (int i = 0; i < ndofs; i++)
      {
         for (int j = 0; j < i; j++)
         {
            double aij = elmat(i, j), aji = elmat(j, i);
            elmat(i, j) = sigma * aji - aij;
            elmat(j, i) = sigma * aij - aji;
         }
         elmat(i, i) *= (sigma - 1.);
      }
   }

   // populate full matrix following github issue #2909
   full_elmat.SetSize(vdim*(ndof1 + ndof2));
   full_elmat = 0.0;
   for (int d=0; d<vdim; ++d)
   {
      for (int j=0; j<ndofs; ++j)
      {
         int jj = (j < ndof1) ? j + d*ndof1 : j - ndof1 + d*ndof2 + vdim*ndof1;
         for (int i=0; i<ndofs; ++i)
         {
            int ii = (i < ndof1) ? i + d*ndof1 : i - ndof1 + d*ndof2 + vdim*ndof1;
            full_elmat(ii, jj) += elmat(i, j);
         }
      }
   }
};

// Vector version of DGDirichletLFIntegrator
class VectorDGDirichletLFIntegrator : public mfem::LinearFormIntegrator
{
protected:
   mfem::VectorCoefficient *uD;
   mfem::Coefficient *Q;
   mfem::MatrixCoefficient *MQ;
   double sigma, kappa;
   int vdim;

   // these are not thread-safe!
   mfem::Vector shape, dshape_dn, nor, nh, ni;
   mfem::DenseMatrix dshape, mq, adjJ;

public:
   VectorDGDirichletLFIntegrator(mfem::VectorCoefficient &u, const double s,
                                 const double k, const int vdim_)
      : uD(&u), Q(NULL), MQ(NULL), sigma(s), kappa(k), vdim(vdim_) { }
   VectorDGDirichletLFIntegrator(mfem::VectorCoefficient &u, mfem::Coefficient &q,
                                 const double s, const double k, const int vdim_)
      : uD(&u), Q(&q), MQ(NULL), sigma(s), kappa(k), vdim(vdim_) { }
   VectorDGDirichletLFIntegrator(mfem::VectorCoefficient &u, mfem::MatrixCoefficient &q,
                                 const double s, const double k, const int vdim_)
      : uD(&u), Q(NULL), MQ(&q), sigma(s), kappa(k), vdim(vdim_) { }

   virtual void AssembleRHSElementVect(const mfem::FiniteElement &el,
                                       mfem::ElementTransformation &Tr,
                                       mfem::Vector &elvect);
   virtual void AssembleRHSElementVect(const mfem::FiniteElement &el,
                                       mfem::FaceElementTransformations &Tr,
                                       mfem::Vector &full_elvect);

   using LinearFormIntegrator::AssembleRHSElementVect;
};

void VectorDGDirichletLFIntegrator::AssembleRHSElementVect(
   const mfem::FiniteElement &el, mfem::ElementTransformation &Tr, mfem::Vector &elvect)
{
   mfem::mfem_error("DGDirichletLFIntegrator::AssembleRHSElementVect");
}

void VectorDGDirichletLFIntegrator::AssembleRHSElementVect(
   const mfem::FiniteElement &el, mfem::FaceElementTransformations &Tr,mfem:: Vector &full_elvect)
{
   int dim, ndof;
   bool kappa_is_nonzero = (kappa != 0.);
   double w;

   dim = el.GetDim();
   ndof = el.GetDof();

   nor.SetSize(dim);
   nh.SetSize(dim);
   ni.SetSize(dim);
   adjJ.SetSize(dim);
   if (MQ)
   {
      mq.SetSize(dim);
   }

   shape.SetSize(ndof);
   dshape.SetSize(ndof, dim);
   dshape_dn.SetSize(ndof);

   mfem::Vector elvect;
   elvect.SetSize(ndof);
   elvect = 0.0;

   mfem::Vector u_DIR;
   u_DIR.SetSize(vdim);
   u_DIR = 0.0;

   full_elvect.SetSize(vdim*ndof);
   full_elvect = 0.0;

   const mfem::IntegrationRule *ir = IntRule;
   if (ir == NULL)
   {
      // a simple choice for the integration order; is this OK?
      int order = 2*el.GetOrder() + 2;
      ir = &mfem::IntRules.Get(Tr.GetGeometryType(), order);
   }

   for (int i = 0; i < vdim; i++)
   {
      // u_DIR = 0.0;
      for (int p = 0; p < ir->GetNPoints(); p++)
      {
         const mfem::IntegrationPoint &ip = ir->IntPoint(p);

         // Set the integration point in the face and the neighboring element
         Tr.SetAllIntPoints(&ip);

         // Access the neighboring element's integration point
         const mfem::IntegrationPoint &eip = Tr.GetElement1IntPoint();

         if (dim == 1)
         {
            nor(0) = 2*eip.x - 1.0;
         }
         else
         {
            CalcOrtho(Tr.Jacobian(), nor);
         }

         el.CalcShape(eip, shape);
         el.CalcDShape(eip, dshape);

         // compute uD through the face transformation
         uD->Eval(u_DIR,Tr,ip);
         w = ip.weight * u_DIR(i) / Tr.Elem1->Weight();
         if (!MQ)
         {
            if (Q)
            {
               w *= Q->Eval(*Tr.Elem1, eip);
            }
            ni.Set(w, nor);
         }
         else
         {
            nh.Set(w, nor);
            MQ->Eval(mq, *Tr.Elem1, eip);
            mq.MultTranspose(nh, ni);
         }
         CalcAdjugate(Tr.Elem1->Jacobian(), adjJ);
         adjJ.Mult(ni, nh);

         dshape.Mult(nh, dshape_dn);
         elvect.Add(sigma, dshape_dn);

         if (kappa_is_nonzero)
         {
            elvect.Add(kappa*(ni*nor), shape); // this is +=, so zero out after each loop
         }
      }
      // copy into full vector
      mfem::Array<int> curr_dofs(ndof);
      curr_dofs = 0;
      for (int j = 0; j < ndof; j++)
      {
         int index = j +
                     i*ndof; // i is loop over vdim the unknown solution vector dimension
         curr_dofs[j] = index;
      }
      full_elvect.SetSubVector(curr_dofs,elvect);

      elvect = 0.0;
   }
}

// boundary terms from grad pressure term
/*
   < {p},[v]*nor >
*/
class DGAvgNormalJumpIntegrator : public mfem::BilinearFormIntegrator
{
private:
   const int vdim;

public:
   DGAvgNormalJumpIntegrator(const int& vdim_) : vdim(vdim_) {};

   void AssembleFaceMatrix(const mfem::FiniteElement &tr_fe1,
                           const mfem::FiniteElement &tr_fe2,
                           const mfem::FiniteElement &te_fe1,
                           const mfem::FiniteElement &te_fe2,
                           mfem::FaceElementTransformations &T,
                           mfem::DenseMatrix &elmat);
};

void DGAvgNormalJumpIntegrator::AssembleFaceMatrix(const mfem::FiniteElement &tr_fe1,
                                                   const mfem::FiniteElement &tr_fe2,
                                                   const mfem::FiniteElement &te_fe1,
                                                   const mfem::FiniteElement &te_fe2,
                                                   mfem::FaceElementTransformations &T,
                                                   mfem::DenseMatrix &elmat)
{
   // test space here is the velocity (vector space), trial space is pressure (scalar space)
   int dim = tr_fe1.GetDim();
   int tr_ndof1, te_ndof1, tr_ndof2, te_ndof2, tr_ndofs, te_ndofs;
   tr_ndof1 = tr_fe1.GetDof();
   te_ndof1 = te_fe1.GetDof();

   if (T.Elem2No >= 0)
   {
      tr_ndof2 = tr_fe2.GetDof();
      te_ndof2 = te_fe2.GetDof();
   }
   else
   {
      tr_ndof2 = 0;
      te_ndof2 = 0;
   }

   tr_ndofs = tr_ndof1 + tr_ndof2;
   te_ndofs = te_ndof1 + te_ndof2;
   elmat.SetSize(te_ndofs*vdim, tr_ndofs);
   elmat = 0.0;

   mfem::Vector ortho(dim), nor(dim);
   mfem::Vector tr_s1(tr_ndof1);
   mfem::Vector tr_s2(tr_ndof2);
   mfem::Vector te_s1(te_ndof1);
   mfem::Vector te_s2(te_ndof2);

   const mfem::IntegrationRule *ir = IntRule;
   if (ir == NULL)
   {
      int order;
      if (tr_ndof2)
      {
         order = 2*(std::max(tr_fe1.GetOrder(), tr_fe2.GetOrder()) + std::max(te_fe1.GetOrder(),
                                                                    te_fe2.GetOrder())) + 2;
      }
      else
      {
         order = 2*(tr_fe1.GetOrder() + te_fe1.GetOrder()) + 2;
      }
      ir = &mfem::IntRules.Get(T.GetGeometryType(), order);
   }

   // elmat = [ A11   A12 ]
   //         [ A21   A22 ]
   // where the blocks corresponds to the terms in the face integral < {p},[v]*nor > from
   // the different elements and trial/test space, i.e.
   // A11 : terms from element 1 test and element 1 trial space
   // A21 : terms from element 2 test and element 1 trial space
   mfem::DenseMatrix A11(te_ndof1*vdim, tr_ndof1);
   mfem::DenseMatrix A12(te_ndof1*vdim, tr_ndof2);
   mfem::DenseMatrix A21(te_ndof2*vdim, tr_ndof1);
   mfem::DenseMatrix A22(te_ndof2*vdim, tr_ndof2);
   double w, detJ;
   for (int n=0; n<ir->GetNPoints(); n++)
   {
      const mfem::IntegrationPoint &ip = ir->IntPoint(n);
      T.SetAllIntPoints(&ip);
      const mfem::IntegrationPoint &eip1 = T.GetElement1IntPoint();
      const mfem::IntegrationPoint &eip2 = T.GetElement2IntPoint();

      // normal
      if (dim == 1)
      {
         nor(0) = 2*eip1.x - 1.0;
      }
      else
      {
         CalcOrtho(T.Jacobian(), nor);
      }

      // normalize nor, see example 18 in general it is not a unit normal
      double normag = 0;
      for (int ii = 0; ii < nor.Size(); ii++)
      {
         normag += nor(ii)*nor(ii);
      }
      normag = sqrt(normag);
      nor *= 1/normag;

      // below if statement needed so on the boundary {p} = p (definition of {} operator)
      if (T.Elem2No >= 0)
      {
         w = ip.weight;
      }
      else
      {
         w = ip.weight*2;
      }
      detJ = T.Weight();

      // calc shape functions in element 1 at current integration point
      tr_fe1.CalcShape(eip1, tr_s1);
      te_fe1.CalcShape(eip1, te_s1);

      // form A11
      for (int d = 0; d<vdim; d++)
      {
         for (int i = 0; i < te_ndof1; i++)
         {
            for (int j = 0; j < tr_ndof1; j++)
            {
               A11(i + te_ndof1*d,j) += 0.5*tr_s1(j)*(te_s1(i)*nor(d))*w*detJ;
            }
         }
      }

      // if element 2 exists form the rest of the blocks
      if (tr_ndof2)
      {
         // calc shape functions in element 2 at current integration point
         tr_fe2.CalcShape(eip2, tr_s2);
         te_fe2.CalcShape(eip2, te_s2);

         // form A12
         for (int d = 0; d<vdim; d++)
         {
            for (int i = 0; i < te_ndof1; i++)
            {
               for (int j = 0; j < tr_ndof2; j++)
               {
                  A12(i + te_ndof1*d,j) += 0.5*tr_s2(j)*(te_s1(i)*nor(d))*w*detJ;
               }
            }
         }

         // form A21
         for (int d = 0; d<vdim; d++)
         {
            for (int i = 0; i < te_ndof2; i++)
            {
               for (int j = 0; j < tr_ndof1; j++)
               {
                  A21(i + te_ndof2*d,j) += -1.0*0.5*tr_s1(j)*(te_s2(i)*nor(d))*w*detJ;
               }
            }
         }

         // form A22
         for (int d = 0; d<vdim; d++)
         {
            for (int i = 0; i < te_ndof2; i++)
            {
               for (int j = 0; j < tr_ndof2; j++)
               {
                  A22(i + te_ndof2*d,j) += -1.0*0.5*tr_s2(j)*(te_s2(i)*nor(d))*w*detJ;
               }
            }
         }
      }
   }

   // populate elmat with the blocks computed above
   elmat.AddMatrix(A11,0,0);
   if (tr_ndof2)
   {
      elmat.AddMatrix(A12, 0, tr_ndof1);
      elmat.AddMatrix(A21, vdim*te_ndof1, 0);
      elmat.AddMatrix(A22, vdim*te_ndof1, tr_ndof1);
   }
}


// boundary terms from grad pressure term
/*
   < {p},[v]*nor > + kappa <{h^{-1} Q} [p],[v]> + sigma <[p],{(Q grad(v)).n}>
*/
class BIntegrator : public mfem::BilinearFormIntegrator
{
private:
   const int vdim;
   mfem::Coefficient *Q;
   mfem::MatrixCoefficient *MQ;
   double sigma, kappa;
   mfem::DenseMatrix jmat,smat,adjJ, te_ds1;
   mfem::Vector ni, nh;

public:
   BIntegrator(const double s, const double k, const int vdim_) : Q(NULL), MQ(NULL), sigma(s), kappa(k), vdim(vdim_) { }
   //BIntegrator(mfem::Coefficient &q, const double s, const double k, const int vdim_) : Q(&q), MQ(NULL), sigma(s), kappa(k), vdim(vdim_) { }
   //BIntegrator(mfem::MatrixCoefficient &q, const double s, const double k, const int vdim_) : Q(NULL), MQ(&q), sigma(s), kappa(k), vdim(vdim_) { }

   void AssembleFaceMatrix(const mfem::FiniteElement &tr_fe1,
                           const mfem::FiniteElement &tr_fe2,
                           const mfem::FiniteElement &te_fe1,
                           const mfem::FiniteElement &te_fe2,
                           mfem::FaceElementTransformations &T,
                           mfem::DenseMatrix &elmat);
};
void BIntegrator::AssembleFaceMatrix(const mfem::FiniteElement &tr_fe1,
                                     const mfem::FiniteElement &tr_fe2,
                                     const mfem::FiniteElement &te_fe1,
                                     const mfem::FiniteElement &te_fe2,
                                     mfem::FaceElementTransformations &T,
                                     mfem::DenseMatrix &elmat)
{

bool kappa_is_nonzero = (kappa != 0.);
double w, wq = 0.0;


// test space here is the velocity (vector space), trial space is pressure (scalar space)
   int dim = tr_fe1.GetDim();
   int tr_ndof1, te_ndof1, tr_ndof2, te_ndof2, tr_ndofs, te_ndofs;
   tr_ndof1 = tr_fe1.GetDof();
   te_ndof1 = te_fe1.GetDof();

   nh.SetSize(dim);
   ni.SetSize(dim);
   adjJ.SetSize(dim);

   if (T.Elem2No >= 0)
   {
      tr_ndof2 = tr_fe2.GetDof();
      te_ndof2 = te_fe2.GetDof();
   }
   else
   {
      tr_ndof2 = 0;
      te_ndof2 = 0;
   }

   tr_ndofs = tr_ndof1 + tr_ndof2;
   te_ndofs = te_ndof1 + te_ndof2;
   elmat.SetSize(te_ndofs*vdim, tr_ndofs);
   elmat = 0.0;
   if (kappa_is_nonzero)
   {
      jmat.SetSize(te_ndofs*vdim, tr_ndofs);
      jmat = 0.;
   }
   smat.SetSize(te_ndofs*vdim, tr_ndofs); // for sigma term
   smat = 0.;

   mfem::Vector ortho(dim), nor(dim);
   mfem::Vector tr_s1(tr_ndof1);
   mfem::Vector tr_s2(tr_ndof2);
   mfem::Vector te_s1(te_ndof1);
   mfem::Vector te_s2(te_ndof2);   
   te_ds1.SetSize(te_ndofs,vdim);
   mfem::Vector te_ds1n(te_ndof1);

   const mfem::IntegrationRule *ir = IntRule;
   if (ir == NULL)
   {
      int order;
      if (tr_ndof2)
      {
         order = 2*(std::max(tr_fe1.GetOrder(), tr_fe2.GetOrder()) + std::max(te_fe1.GetOrder(),
                                                                    te_fe2.GetOrder())) + 2;
      }
      else
      {
         order = 2*(tr_fe1.GetOrder() + te_fe1.GetOrder()) + 2;
      }
      ir = &mfem::IntRules.Get(T.GetGeometryType(), order);
   }

   // elmat = [ A11   A12 ]
   //         [ A21   A22 ]
   // where the blocks corresponds to the terms in the face integral < {p},[v]*nor > from
   // the different elements and trial/test space, i.e.
   // A11 : terms from element 1 test and element 1 trial space
   // A21 : terms from element 2 test and element 1 trial space
   mfem::DenseMatrix A11(te_ndof1*vdim, tr_ndof1);
   mfem::DenseMatrix A12(te_ndof1*vdim, tr_ndof2);
   mfem::DenseMatrix A21(te_ndof2*vdim, tr_ndof1);
   mfem::DenseMatrix A22(te_ndof2*vdim, tr_ndof2);

   for (int n=0; n<ir->GetNPoints(); n++)
   {
      const mfem::IntegrationPoint &ip = ir->IntPoint(n);
      T.SetAllIntPoints(&ip);
      const mfem::IntegrationPoint &eip1 = T.GetElement1IntPoint();
      const mfem::IntegrationPoint &eip2 = T.GetElement2IntPoint();

      // normal
      if (dim == 1)
      {
         nor(0) = 2*eip1.x - 1.0;
      }
      else
      {
         CalcOrtho(T.Jacobian(), nor);
      }

      // normalize nor, see example 18 in general it is not a unit normal
      double normag = 0;
      for (int ii = 0; ii < nor.Size(); ii++)
      {
         normag += nor(ii)*nor(ii);
      }
      normag = sqrt(normag);
      nor *= 1/normag;

      // below if statement needed so on the boundary {p} = p (definition of {} operator)
      if (T.Elem2No >= 0)
      {  
         std::cout <<"nope thats not good, should not exist, you have to take a look into it";
         w = ip.weight/2/T.Weight();
      }
      else
      {
         w = ip.weight/T.Weight();
      }

      // calc ni and wq
      ni.Set(w, nor);

      CalcAdjugate(T.Elem1->Jacobian(), adjJ);
      adjJ.Mult(ni, nh);
      
      // calc shape functions in element 1 at current integration point
      tr_fe1.CalcShape(eip1, tr_s1);
      te_fe1.CalcShape(eip1, te_s1);
      // calc derivative of shape functions in element 1 at current integration point in test space
      te_fe1.CalcDShape(eip1, te_ds1);
      
      te_ds1.Mult(nh,te_ds1n);
      
      if (kappa_is_nonzero)
      {
         wq = ni * nor;
      }

      // form A11
      for (int d = 0; d<vdim; d++)
      {
         for (int i = 0; i < te_ndof1; i++)
         {
            for (int j = 0; j < tr_ndof1; j++)
            {
               A11(i + te_ndof1*d,j) += tr_s1(j)*(te_s1(i)*nor(d))*w;
            }
         }
      }

      // if element 2 exists form the rest of the blocks, should never exist for the pure boundary integral
      if (tr_ndof2)
      {
         std::cout <<"nope thats not good, should not exist, you have to take a look into it";
      }

      if (kappa_is_nonzero)
      {
         wq *= kappa;
         for (int d = 0; d<vdim; d++)
         {
            for (int i = 0; i < te_ndof1; i++)
            {
               for (int j = 0; j < tr_ndof1; j++)
               {
                  const double wsi = wq*tr_s1(j);
                  jmat(i + te_ndof1*d,j) += wsi*te_s1(i);
               }
            }
         }

         if (tr_ndof2) // if element 2 exists form the rest of the blocks, should never exist for the pure boundary integral
         {
            std::cout <<"nope thats not good, should not exist, you have to take a look into it";
         }
      }

      // calc sigma term
      
      for (int d = 0; d<vdim; d++)
      {
         for (int i = 0; i < te_ndof1; i++)
         {
            for (int j = 0; j < tr_ndof1; j++)
            {
               smat(i + te_ndof1*d,j) += tr_s1(j)*te_ds1n(i);
            }
         }
      }

      if (tr_ndof2) // if element 2 exists form the rest of the blocks, should never exist for the pure boundary integral
      {
         std::cout <<"nope thats not good, should not exist, you have to take a look into it";
      }
   }

   // A11 := +A11 + smat + jmat <- if kappa_is_nonzero
   // A11 := +A11 + smat
   if (kappa_is_nonzero)
   {
      A11.AddMatrix(jmat,0,0);
   }
   A11.AddMatrix(smat,0,0);

   // populate elmat with the blocks computed above
   elmat.AddMatrix(A11,0,0);
   if (tr_ndof2) // should not exist for pure boundary integrator
   {
      std::cout << "nooooooooooo";
   }
}