// -*- C++ -*-
// $Id: inline_eigen_maker.h, v0.9 2020-07-22$
/*! \file
 * \brief Inline construction of the H-Wilson and Overlap eigensystem, or just load the existed one.
 *
 * Calculates the eigensystem 
 */

#ifndef __inline_EigenMaker_h__
#define __inline_EigenMaker_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "linearop.h"
#include "meas/inline/io/named_objmap.h"
#include "io/xml_group_reader.h"
#include "actions/ferm/fermstates/fermstates.h"
#include "actions/ferm/fermbcs/fermbc_factory_w.h"
#include <stdio.h>
#include "util/ferm/diractodr.h"
#include "util/ferm/transf.h"
#include "Eigen/Dense"



namespace Chroma 
{ 

	/*! \ingroup inlinehadron */
	namespace InlineEigenMakerEnv 
	{
		extern const std::string name;
		bool registerAll();
	}

	//! Parameter structure
	/*! \ingroup inlinehadron */
	struct InlineEigenMakerParams 
	{
		// Default constructor
		InlineEigenMakerParams(void);
		// Construct from XML
		InlineEigenMakerParams(XMLReader& xml_in, const std::string& path);
		// Write out the configuration of the parameters
		void write(XMLWriter& xml_out, const std::string& path);

		// How often should I be called during a gauge evolution?
		unsigned long frequency;

		// Holds the non-lattice parameters
		struct Param_t
		{
			GroupXML_t      fermact; 
			int  Noeigen;   /*<! the amount of the eigenvectors*/
			int gpu_level; /* 0 for pure CPU, 1 for CPU(vector)+GPU(calc), 2 for pure GPU*/
			bool check_residual; /* flag to check the residual of the data*/
			std::string  filename;          /*!< bb output file name */

			int  extra_space;   /*<! the extra space needed by this maker*/
			double  chebyshev_cut;   /*<! the upper band of the eigenvalues*/
			int  chebyshev_order;   /*<! the maxinum length of t*/
			bool use_ckpoint; /* flag to enable the checkpointing*/

			bool file_exist;
		} param;


		// Holds the names of the lattice objects
		struct NamedObject_t
		{
			std::string eigen_id;  /*!< output eigen id */
			std::string gauge_id;     /*!< input gauge id */
		} named_obj;

		std::string xml_file;  // Alternate XML file pattern
	}; // end of struct InlineEigenMakerParams


	//! Inline measurement of hadron spectrum
	/*! \ingroup inlinehadron */
	class InlineEigenMaker : public AbsInlineMeasurement 
	{
		public:
			// Default destructor
			~InlineEigenMaker() {}
			// Constructor from param struct: copies param struct
			InlineEigenMaker(const InlineEigenMakerParams& p) : params(p) {}
			// Copy constructor
			InlineEigenMaker(const InlineEigenMaker& p) : params(p.params) {}

			// Getter for measurement frequency
			unsigned long getFrequency() const {return params.frequency;}

			//! Sets up the XML and invokes func, which does the acual work
			void operator()( const unsigned long update_no, XMLWriter& xml_out); 

		protected:
			//! Does the actual work
			void func(const unsigned long update_no, XMLWriter& xml_out); 

		private:
			//! The parameter structure; holds names of props, gauge field, XML, etc.
			InlineEigenMakerParams params;
	}; // end of class InlineEigenMaker

	//-----------------------------------//
	//-----------------------------------//

	template <typename T>
		struct EigenPair{
			T vec;
			Complex val;
			Real residual;
		};

	template <typename T>
		class EigenOperator: public UnprecLinearOperator<T, 
		multi1d<LatticeColorMatrix>,
		multi1d<LatticeColorMatrix> >, public multi1d<EigenPair<T> >
	{
		public:
			typedef multi1d<LatticeColorMatrix>  P;
			typedef multi1d<LatticeColorMatrix>  Q;    

			EigenOperator(GroupXML_t &fermact,multi1d<LatticeColorMatrix> &u,int neig):multi1d<EigenPair<T> >(neig)
			{
				std::istringstream  is(fermact.xml);
				XMLReader  paramtop(is);
				read(paramtop, "FermAct", fermact_id);
				bc=new SimpleFermBC<T,P,Q>(SimpleFermBCParams(paramtop, "FermionBC"));
				fs=new SimpleFermState<T,P,Q>(bc,u);
				KYToDRMat();
			}

			//eigen maker start			
			void remap(LatticeFermion& chi, const LatticeFermion& psi){
				EigenOperator<T> &es=*this;
				es(chi,psi,PLUS);
				if(chirality == 0) es(chi,chi,PLUS);
				else chi = (chi+chirality*(Gamma(15)*chi))/2;
				double factor = 1/(max_eigen-cut_eigen);
				chi = 2*factor*chi - factor*(max_eigen+cut_eigen)*psi;	
			}

			void chebyshev_op(LatticeFermion& chi, const LatticeFermion& psi){
				int degree = order;
				LatticeFermion rt0;
				LatticeFermion rt1;
				LatticeFermion rt2;
				rt0 = psi;
				LatticeFermion *t0 = &rt0;
				LatticeFermion *t1 = &rt1;
				LatticeFermion *t2 = &rt2;
				remap(*t1,*t0);
				while(degree > 1){
					remap(*t2,*t1);
					*t2 = 2 * (*t2) - *t0;
    				LatticeFermion *tmp = t1;
    				t1 = t2;
    				t2 = t0;
    				t0 = tmp;
    				--degree;
				}
				chi = *t1;
			}

			int arnoldi_eigensystem(InlineEigenMakerParams &Params, Real kappa, int chiral=0)
    		{	
				StopWatch shijian;
    			shijian.reset();
   	 			shijian.start();

				chirality = chiral;
				if(chirality == 0) error = 1e-14;
				else error = 1e-11;
        		dim = Params.param.extra_space +  Params.param.Noeigen;
        		Noeigen = Params.param.Noeigen;
				if(chirality == 0){
					cut_eigen = pow(Params.param.chebyshev_cut,2);
					max_eigen = pow(1+8*toDouble(kappa),2);
				}
				else{
					cut_eigen = Params.param.chebyshev_cut;
					max_eigen = 8 - 1/toDouble(kappa);
				}
				order = Params.param.chebyshev_order;

        		TOLERANCE = 1.0e-15;
				LOCK_THRESHOLD = 1.0e-20;
        		Maxiteration = 10000;
        		MINFILTER = 5;
        		symflag = true;

        		H.resize(dim,dim);
        		for(int i=0;i<dim;i++){
            		for(int j=0;j<dim;j++){
                		H(i,j)=0;
            		}
        		}
        		Qv.resize(dim);
        		valtemp.resize(dim);
        		Index.resize(dim);
        		for(int i=0;i<dim;i++){
            		Index[i]=i;
        		}
        		random(b);
				if(chirality != 0) b=(b+chirality*(Gamma(15)*b))/2;

				EigenOperator<T> &es=*this;
				int converge,iteration,iNlocked,iNtolock,iNbuffer;
				iNlocked = 0;
				iNbuffer = 0;
				double normb,proof;
				ctemp1 = 1;
				ctemp2 = 1;
				bool yn = false;
				int nvec = 0;
				converge = arnoldi_factorization(nvec,dim);
				shijian.stop();
				QDPIO::cout << "initial arnoldi factorization" << ": total time = "
                			<< shijian.getTimeInSeconds()
                			<< " secs" << std::endl;
				if(converge == dim){
					iteration = dim - nvec;
					for(;;){
						//main loop

						//step1
						shijian.reset();
						shijian.start();
						take_matrix(dim-iNlocked,dim-iNlocked,iNlocked,iNlocked);
						eigensystem();
						for(int i=0;i<dim-iNlocked;i++){
							valtemp[i] = Htemp(i,i);
							Index[i] = i;
						}
						eigsort(0,dim-1-iNlocked);
						for(int i=0; i<dim-iNlocked; i++){ 
      					//if(Layout::primaryNode() && chirality != 0) printf("e[%4d] = %18.10le + i %18.10le \t %d\n", i, valtemp[Index[i]].real(), valtemp[Index[i]].imag(), Index[i]);
						}
						shijian.stop();
						QDPIO::cout << "main loop step1" << ": total time = "
                					<< shijian.getTimeInSeconds()
                					<< " secs" << std::endl;

						//step2
						shijian.reset();
						shijian.start();
						normb = toDouble(RealD(sqrt(norm2(b))));
						converge = 0;
						iNtolock = 0;
						for(int i=0;i<Noeigen-iNlocked;i++){
							proof = normb * abs(Temp(dim-1-iNlocked,Index[i]));
							if(proof < error*abs(valtemp[Index[i]])) converge++;
							if(proof < LOCK_THRESHOLD) iNtolock++;
							//if(Layout::primaryNode() && chirality != 0) printf("res[%4d] = %18.10le crit = %18.10le normb = %18.10le error = %18.10le\n", i, proof, error * abs(valtemp[Index[i]]), normb, error);
						}
						if(Layout::primaryNode() && chirality != 0) printf("Converged=%10d, Iteration=%10d\n",converge+iNlocked,iteration);
						if(converge == Noeigen-iNlocked){
							converge += iNlocked;
							break;
						}
						if(iteration >= Maxiteration){
							converge += iNlocked;
							converge *= -1;
							break;
						}
						shijian.stop();
                		QDPIO::cout << "main loop step2" << ": total time = "
                            		<< shijian.getTimeInSeconds()
                            		<< " secs" << std::endl;

						//step3
						shijian.reset();
						shijian.start();
						iNtolock = 0;
						if(dim-Noeigen-iNbuffer-iNtolock >= MINFILTER) iNbuffer += iNtolock;
						else iNbuffer = dim-Noeigen-MINFILTER;
						Temp.resize(dim,dim);
						for(int i=0;i<dim;i++){
							for(int j=0;j<dim;j++){
								if(i==j) Temp(i,j) = 1;
								else Temp(i,j) = 0;
							}
						}
						for(int i=0;i<dim-Noeigen-iNbuffer;i++){
							for(int iCol=0;iCol<dim;iCol++) H(iCol,iCol) -= valtemp[Index[dim-1-i-iNlocked]];
							QR(true,dim);
							for(int iCol=0;iCol<dim;iCol++) H(iCol,iCol) += valtemp[Index[dim-1-i-iNlocked]];
						}
						iNlocked += iNtolock;
						vm_multiply(dim,Noeigen+1+iNbuffer);
						ctemp1.elem().elem().elem().real()= H(Noeigen+iNbuffer,Noeigen-1+iNbuffer).real();
						ctemp1.elem().elem().elem().imag()= H(Noeigen+iNbuffer,Noeigen-1+iNbuffer).imag();
						ctemp2.elem().elem().elem().real()= Temp(dim-1,Noeigen-1+iNbuffer).real();
                		ctemp2.elem().elem().elem().imag()= Temp(dim-1,Noeigen-1+iNbuffer).imag();
						b = ctemp1*Qv[Noeigen+iNbuffer] + ctemp2*b;
						shijian.stop();
                		QDPIO::cout << "main loop step3" << ": total time = "
                            		<< shijian.getTimeInSeconds()
                            		<< " secs" << std::endl;

						//step4
						shijian.reset();
						shijian.start();
						int num = arnoldi_factorization(Noeigen+iNbuffer,dim);
						if(num!=dim) {yn=true;break;}
						iteration += Noeigen-iNbuffer;
						shijian.stop();
                		QDPIO::cout << "main loop step4" << ": total time = "
                            		<< shijian.getTimeInSeconds()
                            		<< " secs" << std::endl;
						//main loop
					}
				}
				else yn = true;

				take_matrix(dim,dim,0,0);
				eigensystem();
				for(int i=0;i<dim;i++){
					valtemp[i]=Htemp(i,i);
					Index[i]=i;
				}
				eigsort(0,dim-1);

				shijian.reset();
				shijian.start();
				vm_multiply(dim,Noeigen);
				shijian.stop();
        		QDPIO::cout << "vector rotation" << ": total time = "
                    		<< shijian.getTimeInSeconds()
                    		<< " secs" << std::endl;
		
				for(int i=0;i<Noeigen;i++){
					es[i].vec=Qv[Index[i]];
				}

				if(Layout::primaryNode()) printf("\nconverge=%d \n",converge);
				if(yn) return 0;
				else return converge;
			}

			int arnoldi_factorization(int indx,int num){
				EigenOperator<T> &es=*this;
				double beta,nrm;
				int start;

				if(indx == 0){
					beta = toDouble(RealD(sqrt(norm2(b))));
					if(beta < TOLERANCE) return 0;
					Qv[0] = b/beta;
					chebyshev_op(b,Qv[0]);
					ctemp1 = innerProduct(Qv[0],b);
					H(0,0).real(toDouble(real(ctemp1)));
					H(0,0).imag(toDouble(imag(ctemp1)));
					b = b - ctemp1*Qv[0];
				}
				else indx--;

				beta = toDouble(RealD(sqrt(norm2(b))));
				for(int i=indx;i<num-1;i++){
					if(beta < TOLERANCE) return i+1;
					Qv[i+1] = b/beta;
					H(i+1,i) = beta;
					chebyshev_op(b,Qv[i+1]);
					nrm = toDouble(RealD(sqrt(norm2(b))));
					if(symflag) start=i;
					else start=0;
					for(int j=start;j<i+2;j++){
						ctemp1 = innerProduct(Qv[j],b);
						H(j,i+1).real(toDouble(real(ctemp1)));
						H(j,i+1).imag(toDouble(imag(ctemp1)));
						b = b - ctemp1*Qv[j];
					}
					beta = toDouble(RealD(sqrt(norm2(b))));

					if(beta < 0.717*nrm){
						for(int j=0;j<i+2;j++){
							ctemp1 = innerProduct(Qv[j],b);
							ctemp3.real(toDouble(real(ctemp1)));
							ctemp3.imag(toDouble(imag(ctemp1)));
							H(j,i+1) += ctemp3;
							b = b - ctemp1*Qv[j];
						}
						beta = toDouble(RealD(sqrt(norm2(b))));
					}
				}
				return num;
			}

			void QR(bool flag, int iDim){
				Eigen::MatrixXcd pcH;
				if(flag) pcH = H;
				else pcH = Htemp;
				Eigen::VectorXcd pcR11(iDim-1),pcR12(iDim-1),pcR21(iDim-1),pcR22(iDim-1);
				std::complex<double> cT11, cT12, cT21, cT22, cTemp1, cTemp2, cTemp3, cU1, cU2;
				double dV;

				for(int iRow = 0; iRow < iDim-1; iRow++){
    				if (abs(pcH(iRow+1,iRow)) < TOLERANCE){
      					pcR11(iRow) = 0;
						pcR12(iRow) = 0;
						pcR21(iRow) = 0;
						pcR22(iRow) = 0;
      					pcH(iRow+1,iRow) = 0;
      					continue;
    				}

    				dV = sqrt(std::norm(pcH(iRow,iRow)) + std::norm(pcH(iRow+1,iRow)));
    				cU1 = pcH(iRow,iRow);
    				dV = (cU1.real()>0) ? dV : -dV;
    				cU1 += dV;
    				cU2 = pcH(iRow+1,iRow);

    				cT11 = std::conj(cU1);
    				cT11 /= dV;
    				pcR11(iRow) = std::conj(cT11);

    				cT12 = std::conj(cU2);
    				cT12 /= dV;
    				pcR12(iRow) = std::conj(cT12);
    
    				cT21 = std::conj(cT12);
    				cTemp1 = std::conj(cU1);
    				cTemp1 /= cU1;
    				cT21 *= cTemp1;
    				pcR21(iRow) = std::conj(cT21);

    				cTemp1 = cU2 / cU1;
    				cT22 = cT12 * cTemp1;
    				pcR22(iRow) = std::conj(cT22);

    				cTemp1 = pcH(iRow,iRow);
    				cTemp2 = cT11 * cTemp1;
    				cTemp3 = cT12 * pcH(iRow+1,iRow);
    				cTemp2 += cTemp3;
    				pcH(iRow,iRow) -= cTemp2;
    				pcH(iRow+1,iRow) = 0;

    				for(int iCol=iRow+1; iCol < iDim; iCol++){
						cTemp1 = pcH(iRow,iCol);
						cTemp2 = cT11 * cTemp1;
						cTemp2 += cT12 * pcH(iRow+1,iCol);
						pcH(iRow,iCol) -= cTemp2;

						cTemp2 = cT21 * cTemp1;
						cTemp2 += cT22 * pcH(iRow+1,iCol);
						pcH(iRow+1,iCol) -= cTemp2;
    				}
				}

  				for(int iCol = 0; iCol < iDim - 1; iCol++){
  					if(abs(pcR11(iCol)) > TOLERANCE){
    					for(int iRow = 0; iRow < iCol+2; iRow++){
      						cTemp1 = pcH(iRow,iCol);
      						cTemp2 = pcR11(iCol) * cTemp1;
      						cTemp2 += pcR12(iCol) * pcH(iRow,iCol+1);
      						pcH(iRow,iCol) -= cTemp2;

      						cTemp2 = pcR21(iCol) * cTemp1;
      						cTemp2 += pcR22(iCol) * pcH(iRow,iCol+1);
      						pcH(iRow,iCol+1) -= cTemp2;
    					}

    					for(int iRow = 0; iRow < iDim; iRow++){
      						cTemp1 = Temp(iRow,iCol);
      						cTemp2 = pcR11(iCol) * cTemp1;
      						cTemp2 += pcR12(iCol) * Temp(iRow,iCol+1);
      						Temp(iRow,iCol) -= cTemp2;

      						cTemp2 = pcR21(iCol) * cTemp1;
      						cTemp2 += pcR22(iCol) * Temp(iRow,iCol+1);
      						Temp(iRow,iCol+1) -= cTemp2;
    					}
  					}
				}
				if(flag) H = pcH;
				else Htemp = pcH;
			}

			void eigensystem(){
				int idim = Htemp.row(0).size();
				std::complex<double> cTemp,cDiscr,cDist1,cDist2,cEigen;
				int iterations = 0;
				for(int iCol=idim-2;iCol>-1;iCol--){
					for(;iterations<Maxiteration;iterations++){
						if(abs(Htemp(iCol+1,iCol))<TOLERANCE){
							Htemp(iCol+1,iCol)=0;
							break;
						}
						cTemp = Htemp(iCol,iCol)-Htemp(iCol+1,iCol+1);
						cTemp *= cTemp;
						cTemp /= 4;

						cDiscr = Htemp(iCol+1,iCol)*Htemp(iCol,iCol+1);
						cDiscr += cTemp;
						cDiscr = sqrt(cDiscr);
						cTemp = Htemp(iCol,iCol) + Htemp(iCol+1,iCol+1);
						cTemp /= 2;
						cDist1 = cTemp + cDiscr;
						cDist1 = cDist1 - Htemp(iCol+1,iCol+1);
						cDist2 = cTemp - cDiscr;
						cDist2 = cDist2 - Htemp(iCol+1,iCol+1);
						if(abs(cDist1)<abs(cDist2)) cEigen = cDist1 + Htemp(iCol+1,iCol+1);
						else cEigen = cDist2 + Htemp(iCol+1,iCol+1);

						for(int iRow=0;iRow<idim;iRow++){
							Htemp(iRow,iRow) -= cEigen;
						}
						QR(false,idim);

						for(int iRow=0;iRow<idim;iRow++){
							Htemp(iRow,iRow) += cEigen;
						}	

					}
				}
				if(Layout::primaryNode()) printf("\nsmall eigensystem iterations=%d \n",iterations);
			}

			void take_matrix(int lx,int ly,int x0,int y0){
				Htemp.resize(lx,ly);
				Temp.resize(lx,ly);
				for(int i=x0;i<lx;i++){
					for(int j=y0;j<ly;j++){
						Htemp(i,j) = H(i,j);
						if(i==j) Temp(i,j) = 1;
						else Temp(i,j) = 0;
					}
				}
			}

			void swap(int i, int j){
				int temp;
				temp = Index[i];
				Index[i] = Index[j];
				Index[j] = temp;
			}

			void eigsort(int left, int right){
				int last;

				if(left>=right) return;

				swap(left,(left+right)/2);
				last=left;
				for(int i=left+1;i<=right;i++){
					if(abs(valtemp[Index[i]])>abs(valtemp[Index[left]])) swap(++last,i);
				}
	
				swap(left,last);
				eigsort(left,last-1);
				eigsort(last+1,right);
			}
	
			void vm_multiply(int idim,int icol){
				std::vector<T> ssr = Qv;
				for(int j=0;j<icol;j++){
					for(int i=0;i<idim;i++){
						ctemp1.elem().elem().elem().real() = Temp(i,Index[j]).real();
            			ctemp1.elem().elem().elem().imag() = Temp(i,Index[j]).imag();
						if(i==0) Qv[Index[j]] = ctemp1*ssr[i];
						else Qv[Index[j]] += ctemp1*ssr[i];
					}
				}
			}

			void reconstruct(){
				StopWatch shijian;
        		shijian.reset();
        		shijian.start();
				EigenOperator<T> &es=*this;
				LatticeFermion vecTemp,vecTemp2;
				double ZEROMODE = 5e-11;
				int offset = 0;
				//determine the number of zero modes
				do{
					es(vecTemp, es[offset].vec, PLUS);
					es[offset].val = innerProduct(es[offset].vec, vecTemp);
				}
				while(toDouble(real(es[offset++].val)) < ZEROMODE);
				offset--;
				//determine chirality
				vecTemp = Gamma(15)*es[0].vec;
				chirality = toDouble(real(innerProduct(es[0].vec, vecTemp))) > 0?1:-1;
				if(Layout::primaryNode()) printf("chirality is %d\n",chirality);
				//get full vec
				for(int i=0; i<Noeigen-offset; i++){
					// Compute |temp> = D|vec>
					es(vecTemp, es[i+offset].vec, PLUS);
					// Compute the partner of |vec> (opposite chirality partner) and normalize it
					vecTemp2 = (vecTemp-chirality*(Gamma(15)*vecTemp))/2;
					vecTemp2 /= RealD(sqrt(norm2(vecTemp2)));
					// Determine the multiplication factor
					ctemp1 = innerProduct(vecTemp2, vecTemp);
					ctemp3.real(toDouble(real(ctemp1)));
					ctemp3.imag(toDouble(imag(ctemp1)));
					ctemp3 = ctemp3/std::abs(ctemp3);
					ctemp2.elem().elem().elem().real() = ctemp3.imag();
					ctemp2.elem().elem().elem().imag() = ctemp3.real();
					// Determine the eigenvectors
					es[offset+i].vec = vecTemp2 + ctemp2*es[offset+i].vec;
				}
				//normalize
				for(int i=0; i<Noeigen-offset; i++){
					es[offset+i].vec /= RealD(sqrt(norm2(es[offset+i].vec)));
				}
				//Reorthogonalize
				for(int i=0;i<Noeigen;i++){
					for(int j=0;j<i;j++){
						ctemp1 = innerProduct(es[j].vec, es[i].vec);
						es[i].vec -= ctemp1*es[j].vec;
					}
					es[i].vec /= RealD(sqrt(norm2(es[i].vec)));
				}
				shijian.stop();
				QDPIO::cout << "reconstruct" << ": total time = "
                    		<< shijian.getTimeInSeconds()
                    		<< " secs" << std::endl;
			}
//eigen maker end

			const std::string &get_fermact()
			{
				return fermact_id;
			}

			Real Residual(EigenPair<LatticeFermion>& op){
				Real res=sqrt(norm2(op.vec));
				LatticeFermion tmp;
				(*this)(tmp,op.vec,PLUS);
				tmp=tmp-op.val*op.vec;
				res=sqrt(norm2(tmp));

				return res;
			}

			int check_residual(Real residual_goal)
			{
				LatticeFermion tmp;
				EigenOperator<T> &op=*this;
				for(int ind=0;ind<this->size();ind++)
				{
                    Real res=sqrt(norm2(op[ind].vec));
                    double bias=1.0-toDouble(res);
                    (*this)(tmp,op[ind].vec,PLUS);
                    tmp=tmp-op[ind].val*op[ind].vec;
                    res=sqrt(norm2(tmp));
					if(Layout::primaryNode()) printf("Norm evec[%4d] bias = %13.5e ", ind, bias);
					if(Layout::primaryNode()) printf(" eval[%4d] = %11.8f + I %10.8f, ||lambda vec - mat vec|| = %10.5e\n",
						ind,toDouble(real(op[ind].val)),
						toDouble(imag(op[ind].val)),
						toDouble(res));
				
				}
			}

			virtual void operator()(T& chi, const T& psi, enum PlusMinus isign) const
			{
				QDPIO::cerr << "The operator in this base class should not be used" <<std::endl;
				QDP_abort(1);
			}

			void save(const std::string filename, bool to_single=true)
			{
				//save as kentuchy format
				//save eigenvector and eigenvalue seperately

				//eigen=multi1d<int,EigenPair<LatticeFermion>>
				//Residues[]=    //the array of residuals

				//eigenvalue
				StopWatch shijian;
        		shijian.reset();
        		shijian.start();    
				EigenOperator<T> &eigen=*this;
//				if(Layout::primaryNode()) 
				{
					std::string file_eval=filename+".eigvals";
					fprintf(stderr,"Saving eigensysterm ...\n");
					FILE* fileEigval=fopen(file_eval.c_str(),"w");
					fprintf(fileEigval, "Eigenvalues and eigenvectors for overlap.\n");
					fprintf(fileEigval, "Each eigenvector is preceded by a line describing the eigenvalue.\n");
					fprintf(fileEigval, "The residue is defined as norm(mat.vec-lambda.vec).\n");
					fprintf(fileEigval, "The format is: a tag EIGV, the real and imaginary part of the eigenvalue and the residue.\n");
					for(int iIndex = 0; iIndex < eigen.size(); iIndex++){
						fprintf(fileEigval, "EIGV %+.15le\t%+.15le\t%.10le\n",
								toDouble(real(eigen[iIndex].val)),
								toDouble(imag(eigen[iIndex].val)),
								toDouble(eigen[iIndex].residual));
					}
					fclose(fileEigval);
				}
				

				//save eigenvector
				multi2d<LatticeComplexD> vec(Nc,Ns);
				LatticeColorVectorD cv;
				BinaryFileWriter bin(filename);
				for(int i=0; i<eigen.size(); i++){
					QDPIO::cout << i << " " ;
					LatticeFermion f=adj(U)*eigen[i].vec;
					for(int spin=0; spin<Ns; ++spin){
						cv=peekSpin(f,spin);
						for(int color=0; color<Nc;++color){
							vec(color,spin)=peekColor(cv,color);

						}
					}

					for(int spin=0; spin<Ns;++spin){
						for(int color=0; color<Nc; ++color){
							LatticeRealD tmp=real(vec(color,spin));
							write(bin,tmp);
						}
					}

					for(int spin=0; spin<Ns; ++spin){
						for(int color=0; color<Nc; ++color){
							LatticeRealD tmp=imag(vec(color,spin));
							write(bin,tmp);
						}
					}
				}//end i  				
				shijian.stop();
                QDPIO::cout << "save" << ": total time = "
                            << shijian.getTimeInSeconds()
                            << " secs" << std::endl;
			}




			void load(std::string filename, bool from_single=true)
			{
				StopWatch shijian;
        		shijian.reset();
        		shijian.start();
				//load eigenvalues
				std::string file_eval=filename+".eigvals";

				FILE* fileEigval = fopen(file_eval.c_str(),"r");
				EigenOperator<T> &eigen=*this;

				char sTemp[100];
				int iIndex=0;
				double Residual;
				double a,b;
				while(fscanf(fileEigval,"%s",sTemp)!=EOF){
					if(strcmp(sTemp,"EIGV")==0){
						fscanf(fileEigval,"%le %le %le",
								&a,&b,
								&Residual);
						eigen[iIndex].val.elem().elem().elem().real()=a;
						eigen[iIndex].val.elem().elem().elem().imag()=b;
						eigen[iIndex].residual.elem().elem().elem().elem()=Residual;
						iIndex++;
					}
				}//end while
				fclose(fileEigval); 

				//load eigenvectors
				BinaryFileReader bin(filename);
				LatticeFermionD f;
				multi2d<LatticeRealD> re(Nc,Ns);
				multi2d<LatticeRealD> im(Nc,Ns);
				LatticeColorVectorD cv;
				for(int i=0; i<eigen.size();i++){
					QDPIO::cout << i << " " ;
					for(int spin=0;spin<Ns;spin++)
						for(int color=0;color<Nc;color++){
							read(bin,re(color,spin));
						}
					for(int spin=0;spin<Ns;spin++){
						for(int color=0;color<Nc;color++){
							read(bin,im(color,spin));
						}
					}
					for(int spin=0; spin<Ns;++spin){
						for(int color=0; color<Nc; color++){
							pokeColor(cv,cmplx(re(color,spin),im(color,spin)),color);
						}
						pokeSpin(f,cv,spin);

					}
					eigen[i].vec=U*f;
					RealD res=sqrt(norm2(f));
				}//end i
				shijian.stop();
                QDPIO::cout << "load" << ": total time = "
                            << shijian.getTimeInSeconds()
                            << " secs" << std::endl;

			}//end load

			const FermBC<T,P,Q>& getFermBC() const {return *bc;}

			const FermState<T,P,Q>& getFermState() const {return *fs;}

			virtual int create_eigen(InlineEigenMakerParams &Params){};

		protected:


			void KYToDRMat()
			{
			
			  U= zero;
			  RealD     foo = RealD(1) / sqrt(RealD(2));
			  ComplexD  one = cmplx( foo,RealD(0));
			  ComplexD mone = cmplx(-foo,RealD(0));
			
			  pokeSpin(U,  one, 0, 1);
			  pokeSpin(U, mone, 0, 3);
			  pokeSpin(U, mone, 1, 0);
			  pokeSpin(U,  one, 1, 2);
			  pokeSpin(U,  one, 2, 1);
			  pokeSpin(U,  one, 2, 3);
			  pokeSpin(U, mone, 3, 0);
			  pokeSpin(U, mone, 3, 2);
			}

			SpinMatrixD U;

			Handle<FermState<T,P,Q> > fs;
			Handle<FermBC<T,P,Q> > bc;
			std::string fermact_id;

			//variables for eigenmaker
			double error,max_eigen,cut_eigen;
			int dim,Noeigen,order;
	
			double LOCK_THRESHOLD,TOLERANCE;
			int MINFILTER,Maxiteration,chirality;
			bool symflag;

			std::vector<LatticeFermion> Qv;
			LatticeFermion b;
			Eigen::MatrixXcd H;
			Eigen::MatrixXcd Htemp;
			Eigen::MatrixXcd Temp;
			std::vector<int> Index;
			std::vector<std::complex<double>> valtemp;
			Complex ctemp1,ctemp2;
			std::complex<double> ctemp3;
	};

	//-----------------------------------//
	//-----------------------------------//

	// Begin of class InlineEigenEraser

	struct InlineEigenEraserParams 
	{
		InlineEigenEraserParams() {};
		InlineEigenEraserParams(XMLReader& xml_in, const std::string& path)
		{
			XMLReader paramtop(xml_in, path);
			read(paramtop,"Frequency",frequency);
			XMLReader input(paramtop, "NamedObject");
			read(input, "eigen_id", named_obj.eigen_id);
			frequency=1;
		}
		void writeXML(XMLWriter& xml_out, const std::string& path)
		{
			push(xml_out, path);
			write(xml_out,"Frequency",frequency);
			push(xml_out, "NamedObject");
			write(xml_out, "eigen_id", named_obj.eigen_id);
			pop(xml_out);
			pop(xml_out);
		}       

		struct NamedObject_t
		{
			std::string  eigen_id;
		} named_obj;

		unsigned long frequency;
	};

	class InlineEigenEraser : public AbsInlineMeasurement 
	{
		public:
			~InlineEigenEraser() {}
			InlineEigenEraser(const InlineEigenEraserParams& p) : params(p) {}

			unsigned long getFrequency(void) const {return params.frequency;}

			//! Do the writing
			void operator()(const unsigned long update_no,
					XMLWriter& xml_out)
			{
				START_CODE();

				push(xml_out, "erase_eigen");
				write(xml_out, "update_no", update_no);

				QDPIO::cout << "EIGEN_ERASER: Eigen erase" << std::endl;

				// Erase the object
				QDPIO::cout << "Attempt to erase object name = "
					<< params.named_obj.eigen_id << std::endl;
				write(xml_out, "object_id", params.named_obj.eigen_id);
				if (TheNamedObjMap::Instance().check(params.named_obj.eigen_id)) {
					EigenOperator<LatticeFermion> *eigen=TheNamedObjMap::Instance().getData<EigenOperator<LatticeFermion>* >
						(params.named_obj.eigen_id);
					delete eigen;
					TheNamedObjMap::Instance().erase(params.named_obj.eigen_id);
				} else {
					QDPIO::cout << "Eigen system: " << params.named_obj.eigen_id
						<< " is not in the map. Cannot delete" << std::endl;
				}
				QDPIO::cout << "EIGEN_ERASER: ran successfully" << std::endl;

				pop(xml_out);  // erase_named_obj

				END_CODE();
			}

		private:
			InlineEigenEraserParams params;
	}; //end of class InlineEigenEraser   


}; // end of namespace Chroma

#endif

