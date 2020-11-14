// -*- C++ -*-
// $Id: hwilson_eigenop.h, v0.9 2020-07-22$
/*! \file
 * \brief The H-Wilson eigensystem.
 *
 * Calculates the eigensystem 
 */

#ifndef __Overlap_eigenOP_h__
#define __Overlap_eigenOP_h__

#include "hwilson_eigenop.h"
#include "chebyshev_coeff.h"

//#define _time_debug_

#define _chey_size_ 30

namespace Chroma 
{ 

class Overlap_param
{
	public:
	std::vector<void *> hw_evec;
	std::vector<double> hw_eval;
	std::vector<std::vector<double> > coef;
	std::vector<int> hw_size;
};


class OverlapEigenOperator: public EigenOperator<LatticeFermion>
{
public:
	typedef LatticeFermion               T;
	typedef multi1d<LatticeColorMatrix>  P;
	typedef multi1d<LatticeColorMatrix>  Q;
	double accuracy;
	std::string topocharge_para;
	Real rho;
	Real kappa;

	Real get_rho(void) {return rho;}

	OverlapEigenOperator(GroupXML_t &fermact,multi1d<LatticeColorMatrix> &u,int neig=0):EigenOperator<LatticeFermion>(fermact,u,neig)
	{   
		std::istringstream  is(fermact.xml); 
		XMLReader  paramtop(is);
		read(paramtop,"Kappa",kappa);	
		rho = 4-0.5/kappa;
		read(paramtop,"HW_eigen",HW_eigen);
		read(paramtop,"Noeigen_hw",Noeigen_hw);
		if(Noeigen_hw>0){
			if(TheNamedObjMap::Instance().check(HW_eigen)){
				Hw=(HwilsonEigenOperator*)TheNamedObjMap::Instance().getData<EigenOperator<LatticeFermion>* >(HW_eigen);
				if(Hw->size()<Noeigen_hw){
					QDPIO::cerr << "Size of " << HW_eigen << " is not large enough (" <<
					Hw->size() <<" v.s. " << Noeigen_hw <<")." << std::endl;
					QDP_abort(1);
				}
				read(paramtop,"accuracy",accuracy);
				read(paramtop,"topocharge_para",topocharge_para);
				QDPIO::cout << "read param finished" << std::endl; fflush(stdout);
			}
			else{
				QDPIO::cerr << HW_eigen << " for overlap is not found." << std::endl;
				QDP_abort(1);
			}
		}
			
		coef.resize(_chey_size_);
		EigenOperator<LatticeFermion> &es=*Hw;
		if(Noeigen_hw!=0){
			cut_sq=pow(toDouble(real(es[Noeigen_hw-1].val/(1+8*kappa))),2);
			LatticeFermion tmp1=es[0].vec,tmp2;
			StopWatch snoop;
			snoop.reset();
			snoop.start();
			for(int i=0;i<Noeigen_hw;i++){
				Complex inner = innerProduct(es[i].vec, tmp1);
				tmp2 -= inner*es[i].vec;
				tmp2 += inner*es[i].vec;
			}
			snoop.stop();
			double time_l=snoop.getTimeInSeconds()/Noeigen_hw;
			snoop.reset();
			snoop.start();
			for(int i=0;i<10;i++){
				HwSq_scaled(tmp1,tmp2,cut_sq);
				tmp1=tmp2;
			}
			snoop.stop();
			double time_h=snoop.getTimeInSeconds()/10;
			print0("high op=%13.3f, low op=%13.4f\n",time_h,time_l);
               
			for(int i=0;i<_chey_size_;i++){
				double prec=0.2*exp(-1.0*i);
				coef[i].set_HwSize(Noeigen_hw,cut_sq,prec,time_h,time_l);
				double cut_sq_tmp=pow(toDouble(real(es[coef[i].get_HwSize()-1].val/(1+8*kappa))),2);
				coef[i].run(cut_sq_tmp,prec);
				print0("prec=%13.3e, cut=%13.3e, size=%6d, order=%6d\n",prec,cut_sq_tmp,coef[i].get_HwSize(),coef[i].size());
			}
		}
		else{
			cut_sq=8.0e-4;
			QDPIO::cout << "WARNING: the sign function will not be well approximated between 0 and 8.0e-4." << std::endl;
			for(int i=0;i<_chey_size_;i++){
				double prec=0.2*exp(-1.0*i);
				coef[i].set_HwSize(Noeigen_hw,cut_sq,prec,0.0,0.0);
				coef[i].run(cut_sq,prec);
			}
		}
		read(paramtop,"use_gpu",use_gpu);
#ifdef BUILD_QUDA
		if(use_gpu>0)
			set_BasicQudaParam(u,quda_inv_param,DEFAULT,false,false,RECONS_12);
		quda_inv_param.kappa=toDouble(kappa);
		quda_inv_param.dslash_type=QUDA_OVERLAP_WILSON_DSLASH;
		quda_inv_param.eigen_size=Noeigen_hw;
		quda_inv_param.eigen_cut=0.3;
		quda_inv_param.krylov_space=Noeigen_hw+50;
		ov_d=newOverlapQuda(&quda_inv_param);
#endif         	 
	}

	void PackforQUDA(Overlap_param &ov_param){
		ov_param.hw_evec.resize(Noeigen_hw);
		ov_param.hw_eval.resize(Noeigen_hw);
		ov_param.coef.resize(_chey_size_);
		ov_param.hw_size.resize(_chey_size_);
		EigenOperator<LatticeFermion> &es=*Hw;
		for(int i=0;i<Noeigen_hw;i++){
#ifndef BUILD_QUDA_DEVIFACE_SPINOR           
			void* spinorIn =(void *)&(es[i].vec.elem(0).elem(0).elem(0).real());
#else
			void* spinorIn = GetMemoryPtr( es[i].vec.getId() );
#endif        
			ov_param.hw_evec[i]=spinorIn;
			ov_param.hw_eval[i]=toDouble(real(es[i].val));
		}
		for(int i=0;i<_chey_size_;i++){
			ov_param.coef[i]=coef[i];
			ov_param.hw_size[i]=coef[i].get_HwSize();
		}
	}

	void HwSq_scaled(LatticeFermion& chi, const LatticeFermion& psi, double cut) const
	{
		T tmp1,tmp2;
		EigenOperator<LatticeFermion> &es=*Hw;
		es(tmp1,psi,PLUS);
		es(tmp2,tmp1,PLUS);
		double sc1=2/((1+8*toDouble(kappa))*(1+8*toDouble(kappa))*(1-cut));
		double sc2=(1+cut)/(1-cut);
		chi = sc1*tmp2-sc2*psi;
	}
    
	void general_dov(LatticeFermion& chi, const LatticeFermion& psi, double k0, double k1,double k2,double prec) const
	{
#ifdef BUILD_QUDA
		if(use_gpu==true){
#ifndef BUILD_QUDA_DEVIFACE_SPINOR           
			void* spinorIn =(void *)&(psi.elem(0).elem(0).elem(0).real());
			void* spinorOut =(void *)&(chi.elem(0).elem(0).elem(0).real());
#else
			void* spinorIn = GetMemoryPtr( psi.getId() );
			void* spinorOut = GetMemoryPtr( chi.getId() );
#endif        
			Overlap_param ov_param;
			PackforQUDA(ov_param);
			ApplyOverlapQuda(spinorOut,spinorIn,&quda_inv_param,k0,k1,k2,prec,ov_d);
		}
		else     
#endif    	    	
		{
			T low=zero;
			T high=zero;
			T tmp=psi;
			EigenOperator<T> &es=*Hw;
			int is=-(int)log(prec/0.2);
			if(is<0)is=0;
			if(is>coef.size()-1)is=coef.size()-1;

			StopWatch snoop;
#ifdef _time_debug_
			snoop.reset();
			snoop.start();
#endif		
			Complex inner;
			for(int i=0;i<coef[is].get_HwSize();i++){
				inner = innerProduct(es[i].vec, tmp);
				Real n0=norm2(tmp),n1=norm2(es[i].vec);
				tmp -= inner*es[i].vec;
				int sign=(toDouble(real(es[i].val))>0)?1:-1;
				low += sign*inner*es[i].vec;
			}	
#ifdef _time_debug_    	    	
			snoop.stop();
			double time_l=snoop.getTimeInSeconds();
			snoop.reset();
			snoop.start();
#endif     	    	
			double cut_sq_tmp=pow(toDouble(real(es[coef[is].get_HwSize()-1].val/(1+8*kappa))),2);
			T bn2=zero,bn1=zero,bn0=zero;
			T *ptmp,*pbn2=&bn2,*pbn1=&bn1,*pbn0=&bn0;	
			for(int i=coef[is].size()-1;i>=1;i--){
				if(i<coef[is].size()-1)HwSq_scaled(high,*pbn1,cut_sq_tmp);
				*pbn2 = 2.0*high - *pbn0 + coef[is][i]*tmp;
				ptmp=pbn0;
				pbn0=pbn1;
				pbn1=pbn2;
				pbn2=ptmp;
			}
			HwSq_scaled(high,*pbn1,cut_sq_tmp);
			*pbn1 = high - *pbn0 + coef[is][0]*tmp;
			es(high,*pbn1,PLUS);
			snoop.stop();	
#ifdef _time_debug_       	    	
			double time_h=snoop.getTimeInSeconds();
			print0("Overlap timer: low %13.4f, high %13.4f\n",time_l, time_h);
#endif     	    	
			tmp = high/(1+8*toDouble(kappa)) + low;
			chi = k0*psi+k1*(Gamma(15)*tmp)+k2*tmp;
		}
	}
    
	void eps(LatticeFermion& chi, const LatticeFermion& psi,double prec) const
	{
		general_dov(chi,psi,0.0,0.0,1.0,prec);
	}
    	
	void eps5(LatticeFermion& chi, const LatticeFermion& psi,double prec) const
	{
		general_dov(chi,psi,0.0,1.0,0.0,prec);
	}
    	 
	void operator()(LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const 
	{
		general_dov(chi,psi,toDouble(rho),toDouble(rho),0.0,accuracy);
	}
    	
	int create_eigen(InlineEigenMakerParams &Params){
		Complex ctemp0;
		T vectemp;
		EigenOperator<T> &es=*this;
		int conv =  arnoldi_eigensystem(Params,kappa,-1);
		reconstruct();
		for(int i=0;i<Noeigen;i++){
			es(vectemp,es[i].vec,PLUS);
			ctemp0 = innerProduct(es[i].vec,vectemp);
			es[i].val.elem().elem().elem().real() = toDouble(real(ctemp0));
			es[i].val.elem().elem().elem().imag() = toDouble(imag(ctemp0));
		}
		return Noeigen;
	}

	~OverlapEigenOperator()
	{
#ifdef BUILD_QUDA
		if(use_gpu) destroyOverlapQuda(ov_d);
#endif
	}

protected:
	bool use_gpu;     
#ifdef BUILD_QUDA
	QudaInvertParam quda_inv_param;
	void* ov_d;
#endif          	        
	std::string HW_eigen;
	int Noeigen_hw;
	HwilsonEigenOperator* Hw;
	multi1d<chebyshev_coef> coef;
	double cut_sq;
};

}//end namespace

#endif
