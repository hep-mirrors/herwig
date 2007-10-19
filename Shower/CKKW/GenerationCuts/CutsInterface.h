// -*- C++ -*-
#ifndef HERWIG_CutsInterface_H
#define HERWIG_CutsInterface_H

#include <vector>

inline double sqr (double x) {
  return x*x;
}

inline double min (double a, double b) {
  return (a < b ? a : b);
}

/** A simple five momentum */
struct FiveMomentum {

  /** E */  
  double E;
  
  /** px */  
  double px;

  /** px */  
  double py;

  /** px */ 
  double pz;
  
  /** The invariant mass squared */ 
  double invMass2;

  /** Wether the momentum is incoming */
  bool initial;
  
};

inline double dotProduct (const FiveMomentum& p, const FiveMomentum& q) {
  return p.E*q.E - (p.px*q.px+p.py*q.py+p.pz*q.pz);
}

FiveMomentum mult (const FiveMomentum&, double);
FiveMomentum add (const FiveMomentum&, const FiveMomentum&);

std::vector<FiveMomentum> jets;

/**
 * Reset the event.
 */
extern "C" void resetEvent_ ();

/**
 * Register a final state parton's momentum.
 * The invariant mass is calculated the momentum
 * components. Units are assumed to be GeV^2.
 */
extern "C" void registerFinalParton_ (double* E, double* px, double* py, double* pz);

/**
 * Register a initial state parton's momentum.
 * It is assumed to be light-like.
 * Units are assumed to be GeV^2.
 */
extern "C" void registerInitialParton_ (double* E, double* px, double* py, double* pz);

#endif
