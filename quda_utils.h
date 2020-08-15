#ifndef __QUDA_utils_h__
#define __QUDA_utils_h__

#include "actions/ferm/invert/quda_solvers/enum_quda_io.h"
#include <quda.h>

namespace Chroma 
{ 

  inline void set_BasicQudaParam(multi1d<LatticeColorMatrix> &u, QudaInvertParam &quda_inv_param,QudaPrecisionType cudaPrecision,bool tuneDslashP, bool verboseP, QudaReconsType cudaReconstruct)
  {
           QudaPrecision_s cpu_prec;
           QudaPrecision_s gpu_prec;

           int s = sizeof( WordType<Real>::Type_t );
           if (s == 4) {
             cpu_prec = QUDA_SINGLE_PRECISION;
           }
           else {
             cpu_prec = QUDA_DOUBLE_PRECISION;
           }

           switch( cudaPrecision ) {
           case HALF:
             gpu_prec = QUDA_HALF_PRECISION;
             break;
           case SINGLE:
             gpu_prec = QUDA_SINGLE_PRECISION;
             break;
           case DOUBLE:
             gpu_prec = QUDA_DOUBLE_PRECISION;
             break;
           default:
             gpu_prec = cpu_prec;
             break;
           }

           QudaGaugeParam q_gauge_param = newQudaGaugeParam();

           const multi1d<int>& latdims = Layout::subgridLattSize();
           q_gauge_param.X[0] = latdims[0];
           q_gauge_param.X[1] = latdims[1];
           q_gauge_param.X[2] = latdims[2];
           q_gauge_param.X[3] = latdims[3];
           q_gauge_param.type = QUDA_SU3_LINKS;
           q_gauge_param.gauge_order = QUDA_QDP_GAUGE_ORDER; 
           q_gauge_param.t_boundary = QUDA_PERIODIC_T;//smearing is independent to BC
           q_gauge_param.cpu_prec = cpu_prec;
           q_gauge_param.cuda_prec = gpu_prec;

           switch(cudaReconstruct ) {
           case RECONS_NONE:
             q_gauge_param.reconstruct = QUDA_RECONSTRUCT_NO;
             break;
           case RECONS_8:
             q_gauge_param.reconstruct = QUDA_RECONSTRUCT_8;
             break;
           case RECONS_12:
             q_gauge_param.reconstruct = QUDA_RECONSTRUCT_12;
             break;
           default:
             q_gauge_param.reconstruct = QUDA_RECONSTRUCT_12;
             break;
           };

           q_gauge_param.cuda_prec_sloppy = gpu_prec;
           q_gauge_param.reconstruct_sloppy = QUDA_RECONSTRUCT_12;

           q_gauge_param.gauge_fix = QUDA_GAUGE_FIXED_NO;
           q_gauge_param.anisotropy = 1.0;//smearing don't care this;

        // etup padding
           multi1d<int> face_size(4);
           face_size[0] = latdims[1]*latdims[2]*latdims[3]/2;
           face_size[1] = latdims[0]*latdims[2]*latdims[3]/2;
           face_size[2] = latdims[0]*latdims[1]*latdims[3]/2;
           face_size[3] = latdims[0]*latdims[1]*latdims[2]/2;

           int max_face = face_size[0];
           for(int i=1; i <=3; i++) {
             if ( face_size[i] > max_face ) {
               max_face = face_size[i];
             }
           }
           q_gauge_param.ga_pad = max_face;
           q_gauge_param.cuda_prec_precondition = gpu_prec;
           q_gauge_param.reconstruct_precondition = QUDA_RECONSTRUCT_12;

           void* gauge[4];
           for(int mu=0; mu < Nd; mu++) {
#ifndef BUILD_QUDA_DEVIFACE_GAUGE
                    gauge[mu] = (void *)&(u[mu].elem(all.start()).elem().elem(0,0).real());
#else
                    gauge[mu] = QDPCache::Instance().getDevicePtr( u[mu].getId() );
                    QDPIO::cout << "MDAGM CUDA gauge[" << mu << "] in = " << gauge[mu] << "\n";
#endif
           }
           QDPIO::cout << "before  loadGaugeQuda" << std::endl;
           loadGaugeQuda((void *)gauge, &q_gauge_param);

           quda_inv_param = newQudaInvertParam();
           quda_inv_param.cpu_prec = cpu_prec;
           quda_inv_param.cuda_prec = gpu_prec;
           quda_inv_param.cuda_prec_sloppy = gpu_prec;
           quda_inv_param.preserve_source = QUDA_PRESERVE_SOURCE_NO;
           quda_inv_param.dirac_order = QUDA_DIRAC_ORDER;
           quda_inv_param.input_location = QUDA_CPU_FIELD_LOCATION;
           quda_inv_param.output_location = QUDA_CPU_FIELD_LOCATION;

           quda_inv_param.dslash_type = QUDA_WILSON_DSLASH;
           quda_inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
           quda_inv_param.Ls = 1;
           quda_inv_param.solution_type = QUDA_MAT_SOLUTION ;

           if( tuneDslashP ) {
             QDPIO::cout << "Enabling Dslash Autotuning" << std::endl;

             quda_inv_param.tune = QUDA_TUNE_YES;
           }
           else {
             QDPIO::cout << "Disabling Dslash Autotuning" << std::endl;

             quda_inv_param.tune = QUDA_TUNE_NO;
           }
           quda_inv_param.sp_pad = 0;
           quda_inv_param.cl_pad = 0;
           if( verboseP ) {   
//             quda_inv_param.verbosity = QUDA_VERBOSE;
             quda_inv_param.verbosity = QUDA_DEBUG_VERBOSE;
           }
           else {
             quda_inv_param.verbosity = QUDA_SUMMARIZE;
           }
  }
  
}
#endif  
