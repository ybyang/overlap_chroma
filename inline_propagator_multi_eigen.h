// -*- C++ -*-
/*! \file
 * \brief Inline construction of propagator
 *
 * Propagator calculations
 */

#ifndef __inline_propagator_multi_h__
#define __inline_propagator_multi_h__

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement.h"
#include "io/qprop_io.h"

namespace Chroma {
/*! \ingroup inlinehadron */
namespace InlinePropagatorMultiEnv {
extern const std::string name;
bool registerAll();
}

//! Parameter structure
/*! \ingroup inlinehadron */
struct InlinePropagatorMultiParams {
  InlinePropagatorMultiParams();
  InlinePropagatorMultiParams(XMLReader &xml_in, const std::string &path);
  void writeXML(XMLWriter &xml_out, const std::string &path);

  unsigned long frequency;

  ChromaProp_t param;

  multi1d<Real> mass;
  int maxiter;
  double cg_error;
  int flag_dc;


  struct NamedObject_t {
    std::string gauge_id;
    std::string source_id;
    std::string eigen_id;
    multi1d<std::string> prop_id;
  } named_obj;

  std::string xml_file; // Alternate XML file pattern
};

//! Inline propagator calculation
/*! \ingroup inlinehadron */
class InlinePropagatorMulti : public AbsInlineMeasurement {
public:
  ~InlinePropagatorMulti() {}
  InlinePropagatorMulti(const InlinePropagatorMultiParams &p) : params(p) {}
  InlinePropagatorMulti(const InlinePropagatorMulti &p) : params(p.params) {}

  unsigned long getFrequency(void) const { return params.frequency; }

  //! Do the measurement
  void operator()(const unsigned long update_no, XMLWriter &xml_out);

protected:
  //! Do the measurement
  void func(const unsigned long update_no, XMLWriter &xml_out);

private:
  InlinePropagatorMultiParams params;
};
}

#endif
