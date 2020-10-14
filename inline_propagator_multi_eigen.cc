/*! \file
 * \brief Inline construction of propagator
 *
 * PropagatorMulti calculations
 */

#include "fermact.h"
#include "inline_propagator_multi_eigen.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "util/info/unique_id.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"
#include "inline_eigen_maker.h"
#include "overlap_eigenop.h"
#include "util/ferm/transf.h"
#include <complex.h>

#include "corr_contract.h"
namespace Chroma {
namespace InlinePropagatorMultiEnv {
namespace {
AbsInlineMeasurement *createMeasurement(XMLReader &xml_in,
                                        const std::string &path) {
  return new InlinePropagatorMulti(InlinePropagatorMultiParams(xml_in, path));
}

//! Local registration flag
bool registered = false;
}

const std::string name = "PROPAGATOR_MULTI";

//! Register all the factories
bool registerAll() {
  bool success = true;
  if (!registered) {
    success &= WilsonTypeFermActsEnv::registerAll();
    success &= TheInlineMeasurementFactory::Instance().registerObject(
        name, createMeasurement);
    registered = true;
  }
  return success;
}
} // end namespace

//! PropagatorMulti input
void read(XMLReader &xml, const std::string &path,
          InlinePropagatorMultiParams::NamedObject_t &input) {
  XMLReader inputtop(xml, path);

  read(inputtop, "gauge_id", input.gauge_id);
  read(inputtop, "source_id", input.source_id);
  read(inputtop, "prop_id", input.prop_id);
  read(inputtop, "eigen_id", input.eigen_id);
}

//! PropagatorMulti output
void write(XMLWriter &xml, const std::string &path,
           const InlinePropagatorMultiParams::NamedObject_t &input) {
  push(xml, path);

  write(xml, "gauge_id", input.gauge_id);
  write(xml, "source_id", input.source_id);
  write(xml, "prop_id", input.prop_id);
  write(xml, "eigen_id", input.eigen_id);

  pop(xml);
}

// Param stuff
InlinePropagatorMultiParams::InlinePropagatorMultiParams() { frequency = 0; }

InlinePropagatorMultiParams::InlinePropagatorMultiParams(
    XMLReader &xml_in, const std::string &path) {
  try {
    XMLReader paramtop(xml_in, path);

    if (paramtop.count("Frequency") == 1)
      read(paramtop, "Frequency", frequency);
    else
      frequency = 1;

    // Parameters for source construction
    read(paramtop, "Param", param);

    read(paramtop, "mass", mass);
    read(paramtop, "maxiter", maxiter);
    read(paramtop, "cg_error", cg_error);
    read(paramtop, "flag_dc", flag_dc);

    // Read in the output propagator/source configuration info
    read(paramtop, "NamedObject", named_obj);

    if (named_obj.prop_id.size() != mass.size()) {
      QDPIO::cerr << __func__
                  << "Size of propagators doesn't equal to that of the masses"
                  << std::endl;
      QDP_abort(1);
    }

    // Possible alternate XML file pattern
    if (paramtop.count("xml_file") != 0) {
      read(paramtop, "xml_file", xml_file);
    }
  }
  catch (const std::string &e) {
    QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e
                << std::endl;
    QDP_abort(1);
  }
}

void InlinePropagatorMultiParams::writeXML(XMLWriter &xml_out,
                                           const std::string &path) {
  push(xml_out, path);

  write(xml_out, "Param", param);
  write(xml_out, "mass", mass);
  write(xml_out, "NamedObject", named_obj);

  pop(xml_out);
}

// Function call
void InlinePropagatorMulti::operator()(unsigned long update_no,
                                       XMLWriter &xml_out) {
  // If xml file not empty, then use alternate
  if (params.xml_file != "") {
    std::string xml_file = makeXMLFileName(params.xml_file, update_no);

    push(xml_out, "propagator");
    write(xml_out, "update_no", update_no);
    write(xml_out, "xml_file", xml_file);
    pop(xml_out);

    XMLFileWriter xml(xml_file);
    func(update_no, xml);
  } else {
    func(update_no, xml_out);
  }
}

int get_poly_order(void) { return 1; }

void deflation_ov(multi1d<EigenPair<LatticeFermion> > &es,
                  multi1d<LatticeFermion> &source,
                  multi1d<LatticeFermion> &src_tmp) {
  StopWatch swatch;
  swatch.reset();
  swatch.start();

  LatticeFermion vec_tmp;
  for (int d = 0; d < Ns; d++)   // d<4
    for (int c = 0; c < Nc; c++) // c<3
    {
      int iIndex = c + Nc * d;
      src_tmp[iIndex] = source[iIndex];
      for (int i = 0; i < es.size(); i++) {
        DComplex alpha = innerProduct(es[i].vec, src_tmp[iIndex]);
        src_tmp[iIndex] = src_tmp[iIndex] - alpha * es[i].vec;
        vec_tmp = Gamma(Ns * Ns - 1) * es[i].vec;
        alpha = innerProduct(vec_tmp, src_tmp[iIndex]);
        src_tmp[iIndex] = src_tmp[iIndex] - alpha * vec_tmp;
      }
    }
  swatch.stop();
  double secs = swatch.getTimeInSeconds();
  QDPIO::cout << "Deflation time: " << secs << "seconds!" << std::endl;
}

template <typename TVector> struct matmult {
  virtual void operator()(TVector &src, TVector &dest) = 0;
  virtual ~matmult() {}
};
template <typename genvector> struct matmult_prec : matmult<genvector> {
  double max_prec;

  virtual void mult_prec(genvector &src, genvector &dest, double prec) = 0;

  void operator()(genvector &src, genvector &dest, double prec) {
    this->mult_prec(src, dest, prec);
  }

  void operator()(genvector &src, genvector &dest) {
    this->mult_prec(src, dest, max_prec);
  }

  matmult_prec(double m) : max_prec(m) {}
};

template <typename genvector> struct DD_ov_kernel : matmult_prec<genvector> {
  //OverlapEigenOperator &ov;
  EigenOperator<LatticeFermion> &ov;
  Real zero_shift;
  int chirality;

  DD_ov_kernel(EigenOperator<LatticeFermion> &ov, Real zero_shift, int c)
      : matmult_prec<genvector>(0), ov(ov), zero_shift(zero_shift),
        chirality(c) {}

  virtual void mult_prec(genvector &dsrc, genvector &dest, double prec) {
    genvector tmp;
    ov.general_dov(dest, dsrc, toDouble(zero_shift), 1.0, toDouble(chirality),
                   prec);
  }
};

template <typename genvector> int get_chirality(genvector &src) {
  int ichirality;
  double tol = 1e-8;
  genvector tmp;
  // gamma5_multiply (src, tmp);
  tmp = Gamma(15) * src;
  double dtmp =
      // cscalar_product (src, tmp).real / cscalar_product (src, src).real;
      innerProduct(src, tmp).elem().elem().elem().real() /
      innerProduct(src, src).elem().elem().elem().real();
  if (fabs(dtmp - 1) < tol) {
    ichirality = 1;
    printf("Chirality of source vector is positive\n");
  } else if (fabs(dtmp + 1) < tol) {
    ichirality = -1;
    printf("Chirality of source vector is negative\n");
  } else {
    ichirality = 0;
    printf("Source vector is not chiral\n");
  }

  return ichirality;
}

int cgm_adaptive(matmult_prec<LatticeFermion> &matMult, multi1d<Real> shifts,
                 LatticeFermion &src, multi1d<LatticeFermion> &sol, Real error,
                 int maxIters, int fudge = 100) {
  StopWatch swatch;
  swatch.reset();
  swatch.start();
  // shifts[0]==0, the zeroshift
  // shifts.size==masses.size()
  sol.resize(shifts.size());
  sol = zero;
  multi1d<LatticeFermion> p(shifts.size()), r(shifts.size());
  p = src;
  r = src;
  multi1d<Real> R(shifts.size());
  R = zero;
  multi1d<Real> c1(shifts.size()), c2(shifts.size()), c3(shifts.size());
  c1 = 1;
  c2 = 1;
  multi1d<Real> alpha(shifts.size());
  alpha = -1;
  Real alpha_old;
  alpha_old = -1;
  multi1d<Real> gamma(shifts.size());
  gamma = zero;
  Real gamma_old;
  gamma_old = zero;
  Real coef;
  LatticeFermion q;
  q = zero;
  Real zeta = 0;
  int iterCount = 1;
  gamma[0] = sqrt(real(innerProduct(r[0], r[0])));
  gamma_old = gamma[0];
  zeta = 1 / (gamma[0] * gamma[0]);
  error *= gamma[0];
  do {
    if (iterCount > maxIters) {
      QDPIO::cout << "cgm_adaptive: iterCount exceed maxIters!" << std::endl;
      swatch.stop();
      double secs = swatch.getTimeInSeconds();
      QDPIO::cout << "cgm_adaptive took " << secs << " seconds!" << std::endl;
      break;
    }
    if (toDouble(gamma[0]) < toDouble(error)) {
      swatch.stop();
      double secs = swatch.getTimeInSeconds();
      QDPIO::cout << "cgm_adaptive took " << secs << " seconds!" << std::endl;

      break;
    }
    for (int i = 0; i < shifts.size(); i++) {

      if (i == 0) {
        Real current_err = error * sqrt(zeta) / fudge;
        matMult(p[0], q, toDouble(current_err));
        alpha_old = alpha[0];
        alpha[0] = gamma[0] * gamma[0] / (real(innerProduct(q, p[0])));
        r[0] = r[0] - alpha[0] * q;
        sol[0] = sol[0] + alpha[0] * p[0];

        gamma_old = gamma[0];
        gamma[0] = sqrt(real(innerProduct(r[0], r[0])));
        zeta += 1.0 / (gamma[0] * gamma[0]);
        p[0] = r[0] + gamma[0] * gamma[0] / (gamma_old * gamma_old) * p[0];
      } else {
        c1[i] = alpha_old / (alpha[0] * coef * (1 - c1[i]) +
                             alpha_old * (1 + shifts[i] * alpha[0]));
        c2[i] = c2[i] * c1[i];
        alpha[i] = alpha[0] * c1[i];
        sol[i] = sol[i] + alpha[i] * p[i];
        R[i] = (gamma[0] * gamma[0] / (gamma_old * gamma_old)) * c1[i] *
               (alpha[i] / alpha[0]);
        r[i] = c2[i] * r[0];
        p[i] = r[i] + R[i] * p[i];
      }

    } // end i
    coef = (gamma[0] * gamma[0]) / (gamma_old * gamma_old);
    iterCount++;
  } while (true);

  swatch.stop();
  double secs = swatch.getTimeInSeconds();
  QDPIO::cout << "cgm_adaptive took " << secs << " seconds!" << std::endl;

  return 0;
}

template <typename T>
void check_residues(matmult<T> &M, T &src_pos, multi1d<T> &pos_Sol,
                    multi1d<Real> &shifts, int ishifts) {
  T tmp;
  M(pos_Sol[0], tmp);

  Real result = sqrt(real(innerProduct(tmp - src_pos, tmp - src_pos))) /
                sqrt(real(innerProduct(src_pos, src_pos)));
  QDPIO::cout << "The position solution residue for 0 solution=" << result
              << std::endl;
  for (int i = 1; i < ishifts; i++) {
    M(pos_Sol[i], tmp);
    result = sqrt(real(innerProduct(shifts[i] * pos_Sol[i] + tmp - src_pos,
                                    shifts[i] * pos_Sol[i] + tmp - src_pos))) /
             sqrt(real(innerProduct(src_pos, src_pos)));
    QDPIO::cout << "The solution residue for " << i << " solution= " << result
                << std::endl;
  }
}

//int overlap_inverter(OverlapEigenOperator &ov, LatticeFermion &src,
int overlap_inverter(EigenOperator<LatticeFermion> &ov, LatticeFermion &src,
                     multi1d<LatticeFermion> &pvecPropagators,
                     multi1d<Real> &masses, int max_iter, Real error,
                     int one_minus_halfD = 1) {
  StopWatch swatch;
  swatch.reset();
  swatch.start();

  // genvector --> LatticeFermion
  LatticeFermion tmp;
  tmp = zero;
  LatticeFermion tmp2;
  tmp2 = zero;
  LatticeFermion src_pos;
  src_pos = zero;
  LatticeFermion src_neg;
  src_neg = zero;
  multi1d<Real> shifts(masses.size());
  shifts = zero;
  Real zeroshift, rho;
  zeroshift = zero;
  rho = zero;
  int iters = 0;
  rho = ov.get_rho();
  for (int iIndex = 0; iIndex < masses.size(); iIndex++) {
    if (iIndex == 0) {
      shifts[0] = 2 * (rho * rho + masses[0] * masses[0] / 4.0) /
                  (rho * rho - masses[0] * masses[0] / 4.0);
    } else {
      shifts[iIndex] = 2 * (rho * rho + masses[iIndex] * masses[iIndex] / 4.0) /
                           (rho * rho - masses[iIndex] * masses[iIndex] / 4.0) -
                       shifts[0];
    }
  }
  // initializaing the DD^\dagger kernerl
  DD_ov_kernel<LatticeFermion> dd_pos(ov, shifts[0], +1);
  DD_ov_kernel<LatticeFermion> dd_neg(ov, shifts[0], -1);

  // check if src is chiral
  int ichirality = get_chirality(src);
  QDPIO::cout << "ichirality=" << ichirality << std::endl;
  // break src if chirality = 0
  if (ichirality == 0) {
    multi1d<LatticeFermion> pvecExtraSolutions;
    pvecExtraSolutions.resize(masses.size());
    pvecExtraSolutions = zero;

    tmp = Gamma(Ns * Ns - 1) * src;
    // positive chirality
    src_pos = 0.5 * (src + tmp);
    // negative chirality
    src_neg = 0.5 * (src - tmp);

    // positive solutions first then negative
    int pos_iter =
        cgm_adaptive(dd_pos, shifts, src_pos, pvecPropagators, error, max_iter);

    int neg_iter = cgm_adaptive(dd_neg, shifts, src_neg, pvecExtraSolutions,
                                error, max_iter);

    iters = pos_iter + neg_iter;
    check_residues(dd_pos, src_pos, pvecPropagators, shifts, shifts.size());
    check_residues(dd_neg, src_neg, pvecExtraSolutions, shifts, shifts.size());
    iters += shifts.size();
    // check residual
    // check_residues();

    for (int i = 0; i < masses.size(); i++)
      pvecPropagators[i] = pvecPropagators[i] + pvecExtraSolutions[i];

    // erase memory pvecExtraSolutions
  }

  if (ichirality == +1) {
    iters = cgm_adaptive(dd_pos, shifts, src, pvecPropagators, error, max_iter);
    check_residues(dd_pos, src, pvecPropagators, shifts, shifts.size());
  }
  if (ichirality == -1) {
    iters = cgm_adaptive(dd_neg, shifts, src, pvecPropagators, error, max_iter);
    check_residues(dd_neg, src, pvecPropagators, shifts, shifts.size());
  }
  // check residue

  // rescale solutions
  for (int i = 0; i < masses.size(); i++)
    pvecPropagators[i] =
        1 / (rho * rho - masses[i] * masses[i] / 4) * pvecPropagators[i];

  // multiply by D^{\dagger} to get full inversion
  for (int i = 0; i < masses.size(); i++) {
    // LatticeFermion tmp,tmp2;
    multi1d<LatticeFermion> tmp;
    tmp.resize(2);
    tmp(0) = Gamma(Ns * Ns - 1) * pvecPropagators[i];
    ov.general_dov(tmp(1), tmp(0), 0, 0, 1.0, 0);
    pvecPropagators[i] = (rho + masses[i] / 2) * pvecPropagators[i] +
                         (rho - masses[i] / 2) * tmp(1);
  }

  if (one_minus_halfD > 0)
    for (int i = 0; i < masses.size(); i++) {
      if (one_minus_halfD == 2) {
        multi1d<LatticeFermion> tmp(1);
        ov.general_dov(tmp[0], pvecPropagators[i], 1, 0, 0, 0);
        pvecPropagators[i] = 0.5 * (pvecPropagators[i] - tmp(0));
      } else {
        pvecPropagators[i] = (rho / (rho - masses[i] / 2)) *
                             (pvecPropagators[i] - src / (2 * rho));
      }
    }
  swatch.stop();
  QDPIO::cout << "Overlpa_inverter: " << swatch.getTimeInSeconds()
              << " seconds!" << std::endl;

  return iters;
} // end overlap_inverter

//int inverter_h(OverlapEigenOperator &ov, multi1d<LatticeFermion> &src,
int inverter_h(EigenOperator<LatticeFermion> &ov, multi1d<LatticeFermion> &src,
               multi1d<LatticeFermion> &prop, multi1d<Real> &masses,
               int max_iter, double cg_err, int one_minus_halfD) {
  // note that we assumed that prop has been initialized
  StopWatch swatch;
  swatch.reset();
  swatch.start();

  int totaliters = 0, iters;
  QDPIO::cout << "Polynomial order for sign function approximation is"
              << get_poly_order() << std::endl;
  multi1d<LatticeFermion> src_tmp(Nc * Ns);
  src_tmp = zero;

  int nmass = masses.size();
  multi1d<LatticeFermion> vec_prop(masses.size());
  vec_prop = zero;
  deflation_ov(ov, src, src_tmp);

  Real prec0 = cg_err;
  for (int d = 0; d < Ns; d++)
    for (int c = 0; c < Nc; c++) {
      int iIndex = c + Nc * d;
      for (int i = 0; i < nmass; i++)
        vec_prop[i] = prop[iIndex + i * Nc * Ns];
      iters = overlap_inverter(ov, src_tmp[iIndex], vec_prop, masses, max_iter,
                               prec0, one_minus_halfD);
      totaliters += iters;
      for (int i = 0; i < nmass; i++)
        prop[iIndex + i * Nc * Ns] = vec_prop[i];
    }
  printf("CGov: total iterations %d\n", totaliters);
  swatch.stop();
  QDPIO::cout << "Inverter_h: " << swatch.getTimeInSeconds() << " seconds!"
              << std::endl;

  return totaliters;
}

Complex inv(const Complex &lam, Real m, Real rho, int one_minus_halfD = 1) {
  if (one_minus_halfD > 0) {
    return (1.0 - lam / (2.0 * rho)) / (lam + m * (1.0 - lam / (2.0 * rho)));
  } else {
    return 1.0 / (lam + m * (1.0 - lam / (2.0 * rho)));
  }
}

Complex cscalar_product_chi(LatticeFermion &a, LatticeFermion &b, bool left) {
  int off = (left == true) ? 0 : 2;
  Complex res = 0.0;
  for (int d = off; d < off + 2; d++) {
    res += innerProduct(peekSpin(a, d), peekSpin(b, d));
  }

  return res;
}
// rho of inverter_l is not ov.rho
// TODO:why norm2(src[i])==0?
void inverter_l(multi1d<EigenPair<LatticeFermion> > &es,
                multi1d<LatticeFermion> &src, multi1d<LatticeFermion> &prop,
                const multi1d<Real> &mass, Real rho, int one_minus_halfD) {
  StopWatch swatch;
  swatch.reset();
  swatch.start();
QDPIO::cout<<"RH_CHECK!"<<std::endl;
  int nmass = mass.size();
  LatticeFermion src_tmp = zero;
  for (int d = 0; d < 4; d++)
    for (int c = 0; c < 3; c++) {
      int iIndex = c + 3 * d;
      for (int i = 0; i < es.size(); i++) {
        DComplex coef = innerProduct(es[i].vec, src[iIndex]);
        for (int im = 0; im < nmass; im++) {
          prop[iIndex + im * 12] =
              prop[iIndex + im * 12] +
              coef * inv(es[i].val, mass[im], rho, one_minus_halfD) * es[i].vec;
        }
        // TODO  es.prec?
        // This condition juged if it is zero mode
        // if it's nonzero sqrt(abs(eval.real))~0.707abs(eval.imag)
        // if it's zero mode sqrt(abs(eval.real))>100abs(eval.imag)

        if (sqrt(abs(es[i].val.elem().elem().elem().real())) >
            100 * abs(es[i].val.elem().elem().elem().imag()))
          continue;

        src_tmp = Gamma(Ns * Ns - 1) * es[i].vec;
        coef = innerProduct(src_tmp, src[iIndex]);
        for (int im = 0; im < nmass; im++) {
          prop[iIndex + im * 12] =
              prop[iIndex + im * 12] +
              coef * conj(inv(es[i].val, mass[im], rho, one_minus_halfD)) *
                  src_tmp;
        }
      }
    }
  swatch.stop();
  QDPIO::cout << "Inverter_l: " << swatch.getTimeInSeconds() << " seconds"
              << std::endl;
}

void quarkPropMult(multi1d<LatticePropagator> &psi, const multi1d<Real> &masses,
                   const LatticePropagator &chi, int j_decay,
                   const ChromaProp_t &Param_, int &ncg_had) {}

// Real work done here
void InlinePropagatorMulti::func(unsigned long update_no, XMLWriter &xml_out) {
  START_CODE();
  QDPIO::cout << InlinePropagatorMultiEnv::name << ": propagator calculation"
              << std::endl;
  StopWatch snoop;
  snoop.reset();
  snoop.start();
  // Test and grab a reference to the gauge field
  XMLBufferWriter gauge_xml;
  XMLBufferWriter eigen_xml;
  try {
    TheNamedObjMap::Instance().getData<multi1d<LatticeColorMatrix> >(
        params.named_obj.gauge_id);
    TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(
        gauge_xml);
  }
  catch (std::bad_cast) {
    QDPIO::cerr << InlinePropagatorMultiEnv::name
                << ": caught dynamic cast error" << std::endl;
    QDP_abort(1);
  }
  catch (const std::string &e) {
    QDPIO::cerr << InlinePropagatorMultiEnv::name
                << ": std::map call failed: " << e << std::endl;
    QDP_abort(1);
  }
  const multi1d<LatticeColorMatrix> &u =
      TheNamedObjMap::Instance().getData<multi1d<LatticeColorMatrix> >(
          params.named_obj.gauge_id);

/*
  OverlapEigenOperator *eigen =
      TheNamedObjMap::Instance().getData<OverlapEigenOperator *>(
          params.named_obj.op_id);
  OverlapEigenOperator &ov = *eigen;
*/
EigenOperator<LatticeFermion> *eigen=TheNamedObjMap::Instance().getData<EigenOperator<LatticeFermion> *>(
	params.named_obj.eigen_id);
	EigenOperator<LatticeFermion> &ov=*eigen;

  QDPIO::cout << "check residual from mult_prop!" << std::endl;
  ov.check_residual(fuzz);
  push(xml_out, "propagator");
  write(xml_out, "update_no", update_no);
  proginfo(xml_out); // Print out basic program info
  // Write out the input
  params.writeXML(xml_out, "Input");
  // Write out the config header
  write(xml_out, "Config_info", gauge_xml);
  push(xml_out, "Output_version");
  write(xml_out, "out_version", 1);
  pop(xml_out);
  // Calculate some gauge invariant observables just for info.
  MesPlq(xml_out, "Observables", u);
  //
  // Read in the source along with relevant information.
  //
  XMLReader source_file_xml, source_record_xml;
  // These pesky variables are needed in the quarkprop call - only chiral dudes
  // need this stuff, but it must be there for the bleeping virtual function
  // to live at the base class
  int t0;
  int j_decay;
  // Record the type of header
  bool make_sourceP = false;
  bool seqsourceP = false;
  QDPIO::cout << "Snarf the source from a named buffer" << std::endl;
  try {
    // Try the cast to see if this is a valid source
    LatticePropagator &source_tmp =
        TheNamedObjMap::Instance().getData<LatticePropagator>(
            params.named_obj.source_id);
    // Snarf the source info. This is will throw if the source_id is not there
    TheNamedObjMap::Instance().get(params.named_obj.source_id).getFileXML(
        source_file_xml);
    TheNamedObjMap::Instance().get(params.named_obj.source_id).getRecordXML(
        source_record_xml);
    // Try to invert this record XML into a source struct
    // First identify what kind of source might be here
    if (source_record_xml.count("/MakeSource") != 0) {
      make_sourceP = true;
      MakeSourceProp_t orig_header;
      read(source_record_xml, "/MakeSource", orig_header);
      j_decay = orig_header.source_header.j_decay;
      t0 = orig_header.source_header.t_source;
    } else if (source_record_xml.count("/SequentialSource") != 0) {
      seqsourceP = true;
      SequentialSource_t orig_header;
      read(source_record_xml, "/SequentialSource", orig_header);
      j_decay = orig_header.seqsource_header.j_decay;
      t0 = orig_header.seqsource_header.t_sink; // funky, but probably not
                                                // really needed
    } else {
      throw std::string("No appropriate header found");
    }

    // Write out the source header
    write(xml_out, "Source_file_info", source_file_xml);
    write(xml_out, "Source_record_info", source_record_xml);
  }
  catch (std::bad_cast) {
    QDPIO::cerr << InlinePropagatorMultiEnv::name
                << ": caught dynamic cast error" << std::endl;
    QDP_abort(1);
  }
  catch (const std::string &e) {
    QDPIO::cerr << InlinePropagatorMultiEnv::name
                << ": error extracting source_header: " << e << std::endl;
    QDP_abort(1);
  }

  // Should be a valid cast now
  LatticePropagator &quark_prop_source =
      TheNamedObjMap::Instance().getData<LatticePropagator>(
          params.named_obj.source_id);
  QDPIO::cout << "Source successfully read and parsed" << std::endl;
  // Sanity check - write out the norm2 of the source in the Nd-1 direction
  // Use this for any possible verification
  {
    // Initialize the slow Fourier transform phases
    SftMom phases(0, true, Nd - 1);
    multi1d<Double> source_corr =
        sumMulti(localNorm2(quark_prop_source), phases.getSet());
    push(xml_out, "Source_correlator");
    write(xml_out, "source_corr", source_corr);
    pop(xml_out);
  }

  //
  // Loop over the source color and spin, creating the source
  // and calling the relevant propagator routines. The QDP
  // terminology is that a propagator is a matrix in color
  // and spin space
  //
  multi1d<LatticePropagator> quark_propagator(0);
  try {
    quark_propagator.resize(params.mass.size());
  }
  catch (std::bad_cast) {
    QDPIO::cerr << InlinePropagatorMultiEnv::name
                << ": caught dynamic cast error" << std::endl;
    QDP_abort(1);
  }
  catch (const std::string &e) {
    QDPIO::cerr << InlinePropagatorMultiEnv::name
                << ": error creating prop: " << e << std::endl;
    QDP_abort(1);
  }

  int ncg_had = 0;
  //
  // Initialize fermion action
  //
  std::istringstream xml_s(params.param.fermact.xml);
  XMLReader fermacttop(xml_s);
  QDPIO::cout << "FermAct = " << params.param.fermact.id << std::endl;
  //
  // Try the factories
  //
  bool success = false;
  if (!success) {
    try {
      StopWatch swatch;
      swatch.reset();
      QDPIO::cout << "Try the various factories" << std::endl;
      QDPIO::cout << "Compute the multiple quark props" << std::endl;
      swatch.start();

      multi1d<LatticeFermion> prop(params.mass.size() * Nc * Ns), src(Nc * Ns);
      // prop.resize(params.mass.size() * Nc*Ns);
      for (int i = 0; i < Nc * Ns * (params.mass.size()); i++)
        prop[i] = zero;
      for (int d = 0; d < Ns; d++)
        for (int c = 0; c < Nc; c++) {
          PropToFerm(quark_prop_source, src[Nc * d + c], c, d);
        }

      quarkPropMult(quark_propagator, params.mass, quark_prop_source, j_decay,
                    params.param, ncg_had);
      // inverter_h(ov, src, prop, params.mass, 600, 1e-7, 1);
      inverter_h(ov, src, prop, params.mass, params.maxiter, params.cg_error,
                 params.flag_dc);
	QDPIO::cout<<"RH_CHECK!"<<std::endl;
      inverter_l(ov, src, prop, params.mass, ov.get_rho(), params.flag_dc);
      for (int i = 0; i < params.mass.size(); i++) {
        LatticePropagator sol = zero;
        for (int d = 0; d < Ns; d++)
          for (int c = 0; c < Nc; c++) {
            FermToProp(prop[d * Nc + c + Nc * Ns * i], sol, c, d);
          }
        multi1d<DComplex> corr = Corr::contract_meson(
            sol, sol, Gamma(Ns * Ns - 1), Gamma(Ns * Ns - 1));
        for (int t = 0; t < QDP::Layout::lattSize()[Nd - 1]; t++) {
          // for(int t=0;t<8;t++){

          QDPIO::cout << "CORR_ps:\t" << i << "\t" << corr[t] << std::endl;
        }
      }
      swatch.stop();
      QDPIO::cout << "PropagatorMulti computed: time= "
                  << swatch.getTimeInSeconds() << " secs" << std::endl;
      success = true;
    }
    catch (const std::string &e) {
      QDPIO::cout << InlinePropagatorMultiEnv::name
                  << ": caught exception around quarkprop: " << e << std::endl;
    }
  }

  if (!success) {
    QDPIO::cerr << "Error: no fermact found" << std::endl;
    QDP_abort(1);
  }

  push(xml_out, "Relaxation_Iterations");
  write(xml_out, "ncg_had", ncg_had);
  pop(xml_out);
  // Sanity check - write out the propagator (pion) correlator in the Nd-1
  // direction
  {
    // Initialize the slow Fourier transform phases
    SftMom phases(0, true, Nd - 1);
    for (int im = 0; im < params.mass.size(); im++) {
      multi1d<Double> prop_corr =
          sumMulti(localNorm2(quark_propagator[im]), phases.getSet());
      push(xml_out, "Prop_correlator");
      write(xml_out, "prop_corr", prop_corr);
      pop(xml_out);
    }
  }

  // Save the propagator info
  try {
    QDPIO::cout << "Start writing propagator info" << std::endl;
    XMLBufferWriter file_xml;
    push(file_xml, "propagator");
    write(file_xml, "id", uniqueId()); // NOTE: new ID form
    pop(file_xml);
    XMLBufferWriter record_xml;
    if (make_sourceP) {
      MakeSourceProp_t orig_header;
      read(source_record_xml, "/MakeSource", orig_header);
      Propagator_t new_header; // note, abandoning state_info
      new_header.prop_header = params.param;
      new_header.source_header = orig_header.source_header;
      new_header.gauge_header = orig_header.gauge_header;
      write(record_xml, "PropagatorMulti", new_header);
    } else if (seqsourceP) {
      SequentialSource_t orig_header;
      read(source_record_xml, "/SequentialSource", orig_header);
      SequentialProp_t new_header; // note, abandoning state_info
      new_header.seqprop_header = params.param;
      new_header.sink_header = orig_header.sink_header;
      new_header.seqsource_header = orig_header.seqsource_header;
      new_header.forward_props = orig_header.forward_props;
      new_header.gauge_header = orig_header.gauge_header;
      write(record_xml, "SequentialProp", new_header);
    }

    // Write the propagator xml info
    for (int im = 0; im < params.mass.size(); im++) {
      TheNamedObjMap::Instance().create<LatticePropagator>(
          params.named_obj.prop_id[im]);
      TheNamedObjMap::Instance().getData<LatticePropagator>(
          params.named_obj.prop_id[im]) = quark_propagator[im];
      TheNamedObjMap::Instance().get(params.named_obj.prop_id[im]).setFileXML(
          file_xml);
      TheNamedObjMap::Instance().get(params.named_obj.prop_id[im]).setRecordXML(
          record_xml);
    }

    QDPIO::cout << "PropagatorMulti successfully updated" << std::endl;
  }
  catch (std::bad_cast) {
    QDPIO::cerr << InlinePropagatorMultiEnv::name
                << ": caught dynamic cast error" << std::endl;
    QDP_abort(1);
  }
  catch (const std::string &e) {
    QDPIO::cerr << InlinePropagatorMultiEnv::name
                << ": error extracting prop_header: " << e << std::endl;
    QDP_abort(1);
  }

  pop(xml_out); // propagator
  snoop.stop();
  QDPIO::cout << InlinePropagatorMultiEnv::name
              << ": total time = " << snoop.getTimeInSeconds() << " secs"
              << std::endl;
  QDPIO::cout << InlinePropagatorMultiEnv::name << ": ran successfully"
              << std::endl;
  END_CODE();
}
}
