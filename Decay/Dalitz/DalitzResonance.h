// -*- C++ -*-
//
// DalitzResonance.h is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2018 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef HERWIG_DalitzResonance_H
#define HERWIG_DalitzResonance_H

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/EnumIO.h"

namespace Herwig {
using namespace ThePEG;
namespace ResonanceType {

/**
 *  Enum for the type of resonace
 */
enum Type {NonResonant=0,
	   Spin0=1,Spin1=3,Spin2=5,
	   Spin0E691=11,Spin1E691=13,Spin2E691=15,
	   BABARf0=21, Spin0Gauss=31,
	   SpecialSpin0=-1,SpecialSpin1=-3,SpecialSpin2=-5};

}
/**
 *  Struct to contain the properties of the intermediate 
 */
struct DalitzResonance {

  /**
   *  Default constructor
   */
  DalitzResonance() {}

  /**
   *  Constructor specifiying the parameters
   */
  DalitzResonance(long pid, ResonanceType::Type rtype, Energy m, Energy w,
		  unsigned int d1, unsigned int d2, unsigned int s,
		  double mag, double phi, InvEnergy rr)
    : id(pid), type(rtype), mass(m),width(w),
      daughter1(d1),daughter2(d2),spectator(s),
      amp(mag*exp(Complex(0,phi))), R(rr)
  {}

  /**
   *  PID of resonant particle
   */
  long id;

  /**
   *  Type of the resonance
   */
  ResonanceType::Type type;

  /**
   *  Mass of the resonance
   */
  Energy mass;

  /**
   *  Width of the resonance
   */
  Energy width;

  /**
   *  The children
   */
  unsigned int daughter1,daughter2;

  /**
   *   The spectactor
   */
  unsigned int spectator;

  /**
   *  The amplitude
   */
  Complex amp;

  /**
   *  Radius for the Ballt-Weisskopf formfactor
   */
  InvEnergy R;
};


/** 
 * Output operator to allow the structure to be persistently written
 * @param os The output stream
 * @param x The resonance
 */
inline PersistentOStream & operator<<(PersistentOStream & os, 
				      const DalitzResonance  & x) {
  os << x.id << oenum(x.type) << ounit(x.mass,GeV) << ounit(x.width,GeV)
     << x.daughter1 << x.daughter2 << x.spectator
     << x.amp << ounit(x.R,1./GeV);
  return os;
}

/** 
 * Input operator to allow the structure to be persistently read
 * @param is The input stream
 * @param x The resonance
 */
inline PersistentIStream & operator>>(PersistentIStream & is, 
				      DalitzResonance  & x) {
  is >> x.id >> ienum(x.type) >> iunit(x.mass,GeV) >> iunit(x.width,GeV)
     >> x.daughter1 >> x.daughter2 >> x.spectator
     >> x.amp >> iunit(x.R,1./GeV);
  return is;
}
}

#endif
