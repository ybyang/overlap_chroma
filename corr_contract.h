#include "subset.h"
namespace Corr {
multi1d<DComplex> contract_meson(LatticePropagator prop1,
                                 LatticePropagator prop2, Gamma g1, Gamma g2) {
  multi1d<DComplex> corr;
  corr.resize(QDP::Layout::lattSize()[Nd - 1]);

  LatticeComplex tmp = trace(
      g1 * adj(Gamma(Ns * Ns - 1) * prop1 * Gamma(Ns * Ns - 1)) * g2 * prop2);

  Set tslice;
  tslice.make(TimeSliceFunc(Nd - 1));
  corr = sumMulti(tmp, tslice);
  return corr;
}
}
