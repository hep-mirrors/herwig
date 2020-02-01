// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CLEOD0toK0PipPim class.
//

#include "CLEOD0toK0PipPim.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

CLEOD0toK0PipPim::CLEOD0toK0PipPim() : WeakDalitzDecay(5./GeV,1.5/GeV,true)
{}

IBPtr CLEOD0toK0PipPim::clone() const {
  return new_ptr(*this);
}

IBPtr CLEOD0toK0PipPim::fullclone() const {
  return new_ptr(*this);
}

void CLEOD0toK0PipPim::persistentOutput(PersistentOStream & ) const {
}

void CLEOD0toK0PipPim::persistentInput(PersistentIStream & , int) {
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<CLEOD0toK0PipPim,WeakDalitzDecay>
describeHerwigCLEOD0toK0PipPim("Herwig::CLEOD0toK0PipPim", "CLEOD0toK0PipPim.so");

void CLEOD0toK0PipPim::Init() {

  static ClassDocumentation<CLEOD0toK0PipPim> documentation
    ("The CLEOD0toK0PipPim class implements the model of CLEO for"
     " D0 -> Kbar0 pi+pi-",
     "The CLEO fit of \\cite{Muramatsu:2002jp} was used for the decay $D^0\\to\\bar{K}^0\\pi^+\\pi^-$",
     "\\bibitem{Muramatsu:2002jp} H.~Muramatsu {\\it et al.}  "
     "[CLEO Collaboration],Phys.\\ Rev.\\ Lett.\\  {\\bf 89} (2002) 251802"
     "[Erratum-ibid.\\  {\\bf 90} (2003) 059901] [arXiv:hep-ex/0207067].\n");

}

void CLEOD0toK0PipPim::doinit() {
  WeakDalitzDecay::doinit();
  static const double degtorad = Constants::pi/180.;
  // create the resonances
  addResonance(DalitzResonance(getParticleData(    323),0.89166*GeV,0.0508*GeV,0,1,2,-0.11  , 321*degtorad));
  addResonance(DalitzResonance(getParticleData(    113),0.7693 *GeV,0.1502*GeV,1,2,0,-1.    , 0           ));
  addResonance(DalitzResonance(getParticleData(    223),0.78257*GeV,  8.44*MeV,1,2,0,-0.037 , 114*degtorad));
  addResonance(DalitzResonance(getParticleData(   -323),0.89166*GeV,0.0508*GeV,0,2,1,-1.56  , 150*degtorad));
  addResonance(DalitzResonance(getParticleData(9010221),0.977  *GeV,0.050 *GeV,1,2,0, 0.34  , 188*degtorad));
  addResonance(DalitzResonance(getParticleData(    225),1.2754 *GeV,0.1851*GeV,1,2,0, 0.7   , 308*degtorad));
  addResonance(DalitzResonance(getParticleData(  10211),1.310  *GeV,0.2720*GeV,1,2,0, 1.8   ,  85*degtorad));
  addResonance(DalitzResonance(getParticleData( -10321),1.412  *GeV,0.294 *GeV,0,2,1, 2.0   ,   3*degtorad));
  addResonance(DalitzResonance(getParticleData(   -325),1.4256 *GeV,0.0985*GeV,0,2,1, 1.0   , 335*degtorad));
  addResonance(DalitzResonance(getParticleData( -30323),1.717  *GeV,0.322 *GeV,0,2,1,-5.6   , 174*degtorad));
  // D0 -> K- pi+ pi0
  createMode(getParticleData(ParticleID::D0),
	     {getParticleData(ParticleID::Kbar0),
		 getParticleData(ParticleID::piplus),
		 getParticleData(ParticleID::piminus)});
}

void CLEOD0toK0PipPim::doinitrun() {
  WeakDalitzDecay::doinitrun();
}

int CLEOD0toK0PipPim::modeNumber(bool & cc,tcPDPtr parent,
				 const tPDVector & children) const {
  int id0(parent->id());
  // incoming particle must be D0
  if(abs(id0)!=ParticleID::D0) return -1;
  cc = id0==ParticleID::Dbar0;
  // must be three decay products
  if(children.size()!=3) return -1;
  tPDVector::const_iterator pit = children.begin();
  unsigned int npip(0),npim(0),nkm(0),nk0(0),npi0(0);
  for( ;pit!=children.end();++pit) {
    id0=(**pit).id();
    if(id0          ==ParticleID::piplus)  ++npip;
    else if(id0     ==ParticleID::pi0)     ++npi0;
    else if(id0     ==ParticleID::piminus) ++npim;
    else if(abs(id0)==ParticleID::K0)      ++nk0;
    else if(id0     ==ParticleID::K_L0)    ++nk0;
    else if(id0     ==ParticleID::K_S0)    ++nk0;
    else if(abs(id0)==ParticleID::Kplus)   ++nkm;
  }
  if(npim==1&&npip==1&&nk0==1) return  0;
  else                        return -1;
}

Complex CLEOD0toK0PipPim::amplitude(int ichan) const {
  Complex output(0.);
  static const Complex ii(0.,1.);
  unsigned int imin=0, imax(resonances().size());
  if(ichan>=0) {
    imin=ichan;
    imax=imin+1;
  }
  for(unsigned int ix=imin;ix<imax;++ix) {
    // all resonances bar f0(980)
    if(resonances()[ix].resonance->id()!=9010221) {
      output += resAmp(ix);
    }
    // special treatment (Flatte) for f0(980)
    else {
      // output += resAmp(ix);
      static const double gpi=0.09, gK=0.02;
      const Energy & mAB = mInv(resonances()[ix].daughter1,resonances()[ix].daughter2);
      Energy Gamma_pi = gpi*sqrt(0.25*sqr(mAB)-sqr(mOut(resonances()[ix].daughter1)));
      Energy2 arg = 0.25*sqr(mAB)-sqr(mOut(resonances()[ix].spectator));
      complex<Energy> Gamma_K  = arg>=ZERO ? gK*sqrt(arg) : gK*ii*sqrt(-arg);
      output += resonances()[ix].amp*GeV2/
      	(sqr(resonances()[ix].mass)-sqr(mAB)-ii*resonances()[ix].mass*(Gamma_pi+Gamma_K));
    }
  }
  if(ichan<0) output += 1.1*Complex(cos(5.93411945678072),sin(5.93411945678072));
  return output;
}
