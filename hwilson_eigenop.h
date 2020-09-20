// -*- C++ -*-
// $Id: hwilson_eigenop.h, v0.9 2020-07-22$
/*! \file
 * \brief The H-Wilson eigensystem.
 *
 * Calculates the eigensystem 
 */

#ifndef __Hwilson_eigenOP_h__
#define __Hwilson_eigenOP_h__

#include "inline_eigen_maker.h"
#include "actions/ferm/linop/dslash_w.h"

#ifdef BUILD_QUDA
#include "quda_utils.h"
#endif

namespace Chroma 
{ 

  class HwilsonEigenOperator: public EigenOperator<LatticeFermion>
  {
  public:
     
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;
    
    HwilsonEigenOperator(GroupXML_t &fermact,multi1d<LatticeColorMatrix> &u,int neig=0):EigenOperator<LatticeFermion>(fermact,u,neig)
    {
         //read param;
         std::istringstream  is(fermact.xml);
         XMLReader  paramtop(is);
         Real kappa;
         read(paramtop,"Kappa",kappa);
         rho = 4-0.5/kappa;
         
#ifdef BUILD_QUDA
         read(paramtop,"use_gpu",use_gpu);
         if(use_gpu>0)
               set_BasicQudaParam(u,quda_inv_param,DEFAULT,false,false,RECONS_12);
         quda_inv_param.kappa=toDouble(kappa);
#endif         

         AnisoParam_t anisoParam;
         D.create(fs,anisoParam);
    }
     
    void operator()(LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const
    {
#ifdef BUILD_QUDA
         if(use_gpu==true)
         {
#ifndef BUILD_QUDA_DEVIFACE_SPINOR           
           void* spinorIn =(void *)&(psi.elem(0).elem(0).elem(0).real());
           void* spinorOut =(void *)&(chi.elem(0).elem(0).elem(0).real());
#else
           void* spinorIn = GetMemoryPtr( psi.getId() );
           void* spinorOut = GetMemoryPtr( chi.getId() );
#endif        
//             QDPIO::cout << "(GPU)" << "\n";
             ApplyHWilsonQuda(spinorOut,spinorIn,&quda_inv_param);
         }
         else
#endif
         {
             LatticeFermion tmp,tmp2;
             D(tmp,psi,isign);
    	     chi=Gamma(15)*(psi-0.5/(4-rho)*tmp);
         }
    }

    const Real& Rho() const {return rho;}
    
    const FermState<T,P,Q>& getFermState() const {return *fs;}

    int create_eigen(InlineEigenMakerParams &Params)
    {
#ifdef BUILD_QUDA
	    int hw_size=this->size();
            quda_inv_param.dslash_type=QUDA_OVERLAP_WILSON_DSLASH;
            quda_inv_param.kappa=toDouble(1.0/(8-2*rho));
            quda_inv_param.eigen_size=hw_size;
            quda_inv_param.eigen_cut=0.3;
            quda_inv_param.krylov_space=hw_size+50;

	    std::vector<void *> evec(hw_size);
	    std::vector<double> eval(hw_size);
	    EigenOperator<LatticeFermion> &es=*this;
            for(int i=0;i<hw_size;i++)
            {
#ifndef BUILD_QUDA_DEVIFACE_SPINOR           
                     void* spinorIn =(void *)&(es[i].vec.elem(0).elem(0).elem(0).real());
#else
                     void* spinorIn = GetMemoryPtr( es[i].vec.getId() );
#endif
                     evec[i]=spinorIn;
            }
            ov_d=newOverlapQuda(&quda_inv_param,evec.data(),eval.data());
	    for(int i=0;i<hw_size;i++)
		     es[i].val=eval[i];
#endif          
    }     

  protected:    
#ifdef BUILD_QUDA
    bool use_gpu;
    QudaInvertParam quda_inv_param;
    void* ov_d;
#endif      
    Real rho;
    WilsonDslash D;
  };

}
#endif
