#include "CutsInterface.h"

FiveMomentum mult (const FiveMomentum& mom, double fact) {
  FiveMomentum result = mom;
  result.E *= fact;
  result.px *= fact;
  result.py *= fact;
  result.pz *= fact;
  result.invMass2 *= sqr(fact);
  result.initial = mom.initial;
  return result;
}

FiveMomentum add (const FiveMomentum& p, const FiveMomentum& q) {
  FiveMomentum result;
  result.E = p.E+q.E;
  result.px = p.px + q.px;
  result.py = p.py + q.py;
  result.pz = p.pz + q.pz;
  result.initial = p.initial || q.initial;
  return result;
}

void resetEvent_ () {
  jets.clear();
}

void registerFinalParton_ (double* E, double* px, double* py, double* pz) {
  double invMass2 = sqr(*E)-(sqr(*px)+sqr(*py)+sqr(*pz));
  FiveMomentum mom = { *E, *px, *py, *pz, invMass2, false };
  jets.push_back(mom);
}

void registerInitialParton_ (double* E, double* px, double* py, double* pz) {
  FiveMomentum mom = { *E, *px, *py, *pz, 0., true };
  jets.push_back(mom);
}
