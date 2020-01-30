// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CLEOD0toKmPipPi0 class.
//

#include "CLEOD0toKmPipPi0.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

CLEOD0toKmPipPi0::CLEOD0toKmPipPi0() : WeakDalitzDecay(5./GeV,1.5/GeV,true)
{}

IBPtr CLEOD0toKmPipPi0::clone() const {
  return new_ptr(*this);
}

IBPtr CLEOD0toKmPipPi0::fullclone() const {
  return new_ptr(*this);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<CLEOD0toKmPipPi0,WeakDalitzDecay>
describeHerwigCLEOD0toKmPipPi0("Herwig::CLEOD0toKmPipPi0", "HwSMDecay.so");

void CLEOD0toKmPipPi0::Init() {

  static ClassDocumentation<CLEOD0toKmPipPi0> documentation
    ("The CLEOD0toKmPipPi0 class implements the model of CLEO for "
     "D0 -> K- pi+ pi0, Phys. Rev. D63 (2001) 092001.",
     "The CLEO fit of \\cite{Kopp:2000gv} was"
     " used for the decay $D^0\\to K^-\\pi^+\\pi^0$.",
     "\\bibitem{Kopp:2000gv} S.~Kopp {\\it et al.}  [CLEO Collaboration], "
     "Phys.\\ Rev.\\  D {\\bf 63} (2001) 092001 [arXiv:hep-ex/0011065].");

}

void CLEOD0toKmPipPi0::doinit() {
  WeakDalitzDecay::doinit();
  static const double degtorad = Constants::pi/180.;
  // create the resonances
  addResonance(DalitzResonance(getParticleData( 213)  , 0.770*GeV,0.1507*GeV,1,2,0,-1.  , 0.            ));
  addResonance(DalitzResonance(getParticleData(-323)  ,0.8915*GeV,0.050 *GeV,0,2,1,-0.44, 163  *degtorad));
  addResonance(DalitzResonance(getParticleData(-313)  ,0.8961*GeV,0.0505*GeV,0,1,2,-0.39,-  0.2*degtorad));
  addResonance(DalitzResonance(getParticleData(-10321),1.412 *GeV,0.294*GeV ,0,2,1, 0.77,  55.5*degtorad));
  addResonance(DalitzResonance(getParticleData(-10311),1.412 *GeV,0.294*GeV ,0,1,2, 0.85, 166. *degtorad));
  addResonance(DalitzResonance(getParticleData( 30213),1.717 *GeV,0.322*GeV ,1,2,0,-2.50, 171  *degtorad));
  addResonance(DalitzResonance(getParticleData(-30323),1.717 *GeV,0.322*GeV ,0,2,1,-2.50, 103  *degtorad));
  // D+ -> K- pi+ pi+
  createMode(getParticleData(ParticleID::D0),
	     {getParticleData(ParticleID::Kminus),
		 getParticleData(ParticleID::piplus),
		 getParticleData(ParticleID::pi0)});
}

void CLEOD0toKmPipPi0::doinitrun() {
  WeakDalitzDecay::doinitrun();
}

int CLEOD0toKmPipPi0::modeNumber(bool & cc,tcPDPtr parent,
				 const tPDVector & children) const {
  int id0(parent->id());
  // incoming particle must be D0 or D+
  if(abs(id0)!=ParticleID::D0) return -1;
  cc = id0<0;
  // must be three decay products
  if(children.size()!=3) return -1;
  tPDVector::const_iterator pit = children.begin();
  unsigned int npip(0),npim(0),nkm(0),npi0(0);
  int id;
  for( ;pit!=children.end();++pit) {
    id=(**pit).id();
    if(id          ==ParticleID::piplus)  ++npip;
    else if(id     ==ParticleID::pi0)     ++npi0;
    else if(id     ==ParticleID::piminus) ++npim;
    else if(abs(id)==ParticleID::Kplus)   ++nkm;
  }
  if(nkm==1&&(npip+npim)==1&&npi0==1) return  0;
  else                                return -1;
}

Complex CLEOD0toKmPipPi0::amplitude(int ichan) const {
  Complex output(0.);
  unsigned int imin=0, imax(resonances().size());
  if(ichan>=0) {
    imin=ichan;
    imax=imin+1;
  }
  for(unsigned int ix=imin;ix<imax;++ix) {
    output += resAmp(ix);
  }
  if(ichan<0) output += 1.75*Complex(cos(0.5445427266222308),sin(0.5445427266222308));
  return output;
}
