// -*- C++ -*-
// $Id: inline_eigen_maker.h, v0.9 2020-07-22$
/*! \file
 * \brief Inline construction of the H-Wilson and Overlap eigensystem, or just
 *load the existed one.
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

namespace Chroma {

/*! \ingroup inlinehadron */
namespace InlineEigenMakerEnv {
extern const std::string name;
bool registerAll();
}

//! Parameter structure
/*! \ingroup inlinehadron */
struct InlineEigenMakerParams {
  // Default constructor
  InlineEigenMakerParams(void);
  // Construct from XML
  InlineEigenMakerParams(XMLReader &xml_in, const std::string &path);
  // Write out the configuration of the parameters
  void write(XMLWriter &xml_out, const std::string &path);

  // How often should I be called during a gauge evolution?
  unsigned long frequency;

  // Holds the non-lattice parameters
  struct Param_t {
    GroupXML_t fermact;
    int Noeigen;   /*<! the amount of the eigenvectors*/
    int gpu_level; /* 0 for pure CPU, 1 for CPU(vector)+GPU(calc), 2 for pure
                      GPU*/
    bool check_residual;  /* flag to check the residual of the data*/
    std::string filename; /*!< bb output file name */

    int extra_space;      /*<! the extra space needed by this maker*/
    double chebyshev_cut; /*<! the upper band of the eigenvalues*/
    int chebyshev_order;  /*<! the maxinum length of t*/
    bool use_ckpoint;     /* flag to enable the checkpointing*/

    bool file_exist;
  } param;

  // Holds the names of the lattice objects
  struct NamedObject_t {
    std::string eigen_id; /*!< output eigen id */
    std::string gauge_id; /*!< input gauge id */
    std::string op_id;    /* id of Overlap/HwEigenOperator"*/
  } named_obj;

  std::string xml_file; // Alternate XML file pattern
};                      // end of struct InlineEigenMakerParams

//! Inline measurement of hadron spectrum
/*! \ingroup inlinehadron */
class InlineEigenMaker : public AbsInlineMeasurement {
public:
  // Default destructor
  ~InlineEigenMaker() {}
  // Constructor from param struct: copies param struct
  InlineEigenMaker(const InlineEigenMakerParams &p) : params(p) {}
  // Copy constructor
  InlineEigenMaker(const InlineEigenMaker &p) : params(p.params) {}

  // Getter for measurement frequency
  unsigned long getFrequency() const { return params.frequency; }

  //! Sets up the XML and invokes func, which does the acual work
  void operator()(const unsigned long update_no, XMLWriter &xml_out);

protected:
  //! Does the actual work
  void func(const unsigned long update_no, XMLWriter &xml_out);

private:
  //! The parameter structure; holds names of props, gauge field, XML, etc.
  InlineEigenMakerParams params;
}; // end of class InlineEigenMaker

//-----------------------------------//
//-----------------------------------//

template <typename T> struct EigenPair {
  T vec;
  Complex val;
  Real residual;
};

template <typename T>
class EigenOperator
    : public UnprecLinearOperator<T, multi1d<LatticeColorMatrix>,
                                  multi1d<LatticeColorMatrix> >,
      public multi1d<EigenPair<T> > {
public:
  typedef multi1d<LatticeColorMatrix> P;
  typedef multi1d<LatticeColorMatrix> Q;
  // multi1d< EigenPair<LatticeFermion> >  eigenpair;

  EigenOperator(GroupXML_t &fermact, multi1d<LatticeColorMatrix> &u, int neig)
      : multi1d<EigenPair<T> >(neig) {
    std::istringstream is(fermact.xml);
    XMLReader paramtop(is);
    read(paramtop, "FermAct", fermact_id);
    bc = new SimpleFermBC<T, P, Q>(SimpleFermBCParams(paramtop, "FermionBC"));
    fs = new SimpleFermState<T, P, Q>(bc, u);
    KYToDRMat();
  }

  const std::string &get_fermact() { return fermact_id; }

  Real Residual(EigenPair<LatticeFermion> &op) {
    Real res = sqrt(norm2(op.vec));
    LatticeFermion tmp;
    (*this)(tmp, op.vec, PLUS);
    tmp = tmp - op.val * op.vec;
    res = sqrt(norm2(tmp));

    return res;
  }

  int check_residual(Real residual_goal) {
    StopWatch swatch;
    swatch.reset();
    swatch.start();
    LatticeFermion tmp;
    EigenOperator<T> &op = *this;
    for (int ind = 0; ind < this->size(); ind++) {
      Real res = sqrt(norm2(op[ind].vec));
      double bias = 1.0 - res.elem().elem().elem().elem();
      if (Layout::primaryNode())
        (*this)(tmp, op[ind].vec, PLUS);
      tmp = tmp - op[ind].val * op[ind].vec;
      res = sqrt(norm2(tmp));
      if (Layout::primaryNode()) {
        printf("Norm evec[%4d] bias = %13.5e ", ind, bias);
        printf(" eval[%4d] = %11.8f + I %10.8f, ||lambda vec - mat vec|| = "
               "%10.5e\n",
               ind, toDouble(real(op[ind].val)), toDouble(imag(op[ind].val)),
               toDouble(res));
      }
    }
    swatch.stop();
    double secs = swatch.getTimeInSeconds();
    QDPIO::cout << "Checked residue of eigensystem in " << secs << " seconds!"
                << std::endl;
  }

  virtual void operator()(T &chi, const T &psi, enum PlusMinus isign) const {
    QDPIO::cerr << "The operator in this base class should not be used"
                << std::endl;
    QDP_abort(1);
  }

  void save(const std::string filename, bool to_single = true) {
    // save as kentuchy format
    // save eigenvector and eigenvalue seperately

    // eigen=multi1d<int,EigenPair<LatticeFermion>>
    // Residues[]=    //the array of residuals

    // eigenvalue
    StopWatch swatch;
    swatch.reset();
    swatch.start();
    EigenOperator<T> &eigen = *this;
    if (Layout::primaryNode()) {
      std::string file_eval = filename + ".eigvals";
      fprintf(stderr, "Saving eigensysterm ...\n");
      FILE *fileEigval = fopen(file_eval.c_str(), "w");
      fprintf(fileEigval, "Eigenvalues and eigenvectors for overlap.\n");
      fprintf(fileEigval, "Each eigenvector is preceded by a line describing "
                          "the eigenvalue.\n");
      fprintf(fileEigval,
              "The residue is defined as norm(mat.vec-lambda.vec).\n");
      fprintf(fileEigval, "The format is: a tag EIGV, the real and imaginary "
                          "part of the eigenvalue and the residue.\n");
      for (int iIndex = 0; iIndex < eigen.size(); iIndex++) {
        fprintf(fileEigval, "EIGV %+.15le\t%+.15le\t%.10le\n",
                toDouble(real(eigen[iIndex].val)),
                toDouble(imag(eigen[iIndex].val)),
                toDouble(eigen[iIndex].residual));
      }
      fclose(fileEigval);
    }

    // save eigenvector
    multi2d<LatticeComplexD> vec(Nc, Ns);
    LatticeColorVectorD cv;
    BinaryFileWriter bin(filename);
    for (int i = 0; i < eigen.size(); i++) {
      QDPIO::cout << i << " ";
      LatticeFermion f = adj(U) * eigen[i].vec;
      for (int spin = 0; spin < Ns; ++spin) {
        cv = peekSpin(f, spin);
        for (int color = 0; color < Nc; ++color) {
          vec(color, spin) = peekColor(cv, color);
        }
      }

      for (int spin = 0; spin < Ns; ++spin) {
        for (int color = 0; color < Nc; ++color) {
          LatticeRealD tmp = real(vec(color, spin));
          write(bin, tmp);
        }
      }

      for (int spin = 0; spin < Ns; ++spin) {
        for (int color = 0; color < Nc; ++color) {
          LatticeRealD tmp = imag(vec(color, spin));
          write(bin, tmp);
        }
      }
    } // end i
    QDPIO::cout << std::endl;
    swatch.stop();
    double secs = swatch.getTimeInSeconds();
    QDPIO::cout << "Saving eigensystem in " << secs << " seconds!" << std::endl;
  }

  void load(std::string filename, bool from_single = true) {
    // load eigenvalues
    StopWatch swatch;
    swatch.reset();
    swatch.start();
    std::string file_eval = filename + ".eigvals";

    FILE *fileEigval = fopen(file_eval.c_str(), "r");
    EigenOperator<T> &eigen = *this;

    char sTemp[100];
    int iIndex = 0;
    double Residual;
    double a, b;
    while (fscanf(fileEigval, "%s", sTemp) != EOF) {
      if (strcmp(sTemp, "EIGV") == 0) {
        fscanf(fileEigval, "%le %le %le", &a, &b, &Residual);
        eigen[iIndex].val.elem().elem().elem().real() = a;
        eigen[iIndex].val.elem().elem().elem().imag() = b;
        eigen[iIndex].residual.elem().elem().elem().elem() = Residual;
        iIndex++;
      }
    } // end while
    fclose(fileEigval);

    // load eigenvectors
    BinaryFileReader bin(filename);
    LatticeFermionD f;
    multi2d<LatticeRealD> re(Nc, Ns);
    multi2d<LatticeRealD> im(Nc, Ns);
    LatticeColorVectorD cv;
    for (int i = 0; i < eigen.size(); i++) {
      QDPIO::cout << i << " ";
      for (int spin = 0; spin < Ns; spin++)
        for (int color = 0; color < Nc; color++) {
          read(bin, re(color, spin));
        }
      for (int spin = 0; spin < Ns; spin++) {
        for (int color = 0; color < Nc; color++) {
          read(bin, im(color, spin));
        }
      }
      for (int spin = 0; spin < Ns; ++spin) {
        for (int color = 0; color < Nc; color++) {
          pokeColor(cv, cmplx(re(color, spin), im(color, spin)), color);
        }
        pokeSpin(f, cv, spin);
      }
      eigen[i].vec = U * f;
      RealD res = sqrt(norm2(f));
    } // end i
    QDPIO::cout << std::endl;
    swatch.stop();
    double secs = swatch.getTimeInSeconds();
    QDPIO::cout << "Loading eigensystem in " << secs << " seconds" << std::endl;

  } // end load

  const FermBC<T, P, Q> &getFermBC() const { return *bc; }

  const FermState<T, P, Q> &getFermState() const { return *fs; }

  virtual int create_eigen(InlineEigenMakerParams &Params) {};

protected:
  void KYToDRMat() {

    U = zero;
    RealD foo = RealD(1) / sqrt(RealD(2));
    ComplexD one = cmplx(foo, RealD(0));
    ComplexD mone = cmplx(-foo, RealD(0));

    pokeSpin(U, one, 0, 1);
    pokeSpin(U, mone, 0, 3);
    pokeSpin(U, mone, 1, 0);
    pokeSpin(U, one, 1, 2);
    pokeSpin(U, one, 2, 1);
    pokeSpin(U, one, 2, 3);
    pokeSpin(U, mone, 3, 0);
    pokeSpin(U, mone, 3, 2);
  }

  SpinMatrixD U;

  Handle<FermState<T, P, Q> > fs;
  Handle<FermBC<T, P, Q> > bc;
  std::string fermact_id;
};

//-----------------------------------//
//-----------------------------------//

// Begin of class InlineEigenEraser

struct InlineEigenEraserParams {
  InlineEigenEraserParams() {};
  InlineEigenEraserParams(XMLReader &xml_in, const std::string &path) {
    XMLReader paramtop(xml_in, path);
    read(paramtop, "Frequency", frequency);
    XMLReader input(paramtop, "NamedObject");
    read(input, "eigen_id", named_obj.eigen_id);
    frequency = 1;
  }
  void writeXML(XMLWriter &xml_out, const std::string &path) {
    push(xml_out, path);
    write(xml_out, "Frequency", frequency);
    push(xml_out, "NamedObject");
    write(xml_out, "eigen_id", named_obj.eigen_id);
    pop(xml_out);
    pop(xml_out);
  }

  struct NamedObject_t {
    std::string eigen_id;
  } named_obj;

  unsigned long frequency;
};

class InlineEigenEraser : public AbsInlineMeasurement {
public:
  ~InlineEigenEraser() {}
  InlineEigenEraser(const InlineEigenEraserParams &p) : params(p) {}

  unsigned long getFrequency(void) const { return params.frequency; }

  //! Do the writing
  void operator()(const unsigned long update_no, XMLWriter &xml_out) {
    START_CODE();

    push(xml_out, "erase_eigen");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << "EIGEN_ERASER: Eigen erase" << std::endl;

    // Erase the object
    QDPIO::cout << "Attempt to erase object name = "
                << params.named_obj.eigen_id << std::endl;
    write(xml_out, "object_id", params.named_obj.eigen_id);
    if (TheNamedObjMap::Instance().check(params.named_obj.eigen_id)) {
      EigenOperator<LatticeFermion> *eigen =
          TheNamedObjMap::Instance().getData<EigenOperator<LatticeFermion> *>(
              params.named_obj.eigen_id);
      delete eigen;
      TheNamedObjMap::Instance().erase(params.named_obj.eigen_id);
    } else {
      QDPIO::cout << "Eigen system: " << params.named_obj.eigen_id
                  << " is not in the map. Cannot delete" << std::endl;
    }
    QDPIO::cout << "EIGEN_ERASER: ran successfully" << std::endl;

    pop(xml_out); // erase_named_obj

    END_CODE();
  }

private:
  InlineEigenEraserParams params;
}; // end of class InlineEigenEraser

}; // end of namespace Chroma

#endif
