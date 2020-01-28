// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CLEODptoKmPipPipDalitz class.
//

#include "CLEODptoKmPipPipDalitz.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

CLEODptoKmPipPipDalitz::CLEODptoKmPipPipDalitz() : WeakDalitzDecay(5./GeV,1.5/GeV,false)
{}

IBPtr CLEODptoKmPipPipDalitz::clone() const {
  return new_ptr(*this);
}

IBPtr CLEODptoKmPipPipDalitz::fullclone() const {
  return new_ptr(*this);
}

void CLEODptoKmPipPipDalitz::doinit() {
  WeakDalitzDecay::doinit();
  static const double degtorad = Constants::pi/180.;
  // create the resonances
  addResonance(DalitzResonance(getParticleData(-313    ), 896  *MeV, 50.3*MeV,0,1,2,-1.   ,   0.          ));
  addResonance(DalitzResonance(getParticleData(-313    ), 896  *MeV, 50.3*MeV,0,2,1,-1.   ,   0.          ));
  addResonance(DalitzResonance(getParticleData(-10311  ),1463  *MeV,163.8*MeV,0,1,2,3.   ,  49.7*degtorad));
  addResonance(DalitzResonance(getParticleData(-10311  ),1463  *MeV,163.8*MeV,0,2,1,3.   ,  49.7*degtorad));
  addResonance(DalitzResonance(getParticleData(-315    ),1432.4*MeV,  109*MeV,0,1,2,0.962, -29.9*degtorad));
  addResonance(DalitzResonance(getParticleData(-315    ),1432.4*MeV,  109*MeV,0,2,1,0.962, -29.9*degtorad));
  addResonance(DalitzResonance(getParticleData(-30313  ),1717  *MeV,  322*MeV,0,1,2,-6.5  ,  29.0*degtorad));
  addResonance(DalitzResonance(getParticleData(-30313  ),1717  *MeV,  322*MeV,0,2,1,-6.5  ,  29.0*degtorad));
  addResonance(DalitzResonance(getParticleData(-9000311), 809  *MeV,  470*MeV,0,1,2,5.01 ,-163.7*degtorad));
  addResonance(DalitzResonance(getParticleData(-9000311), 809  *MeV,  470*MeV,0,2,1,5.01 ,-163.7*degtorad));
  // D+ -> K- pi+ pi+
  createMode(getParticleData(ParticleID::Dplus),
	     {getParticleData(ParticleID::Kminus),
		 getParticleData(ParticleID::piplus),
		 getParticleData(ParticleID::piplus)});
}

void CLEODptoKmPipPipDalitz::doinitrun() {
  WeakDalitzDecay::doinitrun();
}

void CLEODptoKmPipPipDalitz::persistentOutput(PersistentOStream & os) const {
}

void CLEODptoKmPipPipDalitz::persistentInput(PersistentIStream & is, int) {
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<CLEODptoKmPipPipDalitz,WeakDalitzDecay>
describeHerwigCLEODptoKmPipPipDalitz("Herwig::CLEODptoKmPipPipDalitz", "HwSMDecay.so");

void CLEODptoKmPipPipDalitz::Init() {

  static ClassDocumentation<CLEODptoKmPipPipDalitz> documentation
    ("There is no documentation for the CLEODptoKmPipPipDalitz class");

}


int CLEODptoKmPipPipDalitz::modeNumber(bool & cc,tcPDPtr parent,
				const tPDVector & children) const {
  int id0(parent->id());
  // incoming particle must be D+
  if(abs(id0)!=ParticleID::Dplus) return -1;
  cc = id0<0;
  // must be three decay products
  if(children.size()!=3) return -1;
  unsigned int npip(0),npim(0),nk(0);
  for(tPDPtr child : children) {
    long id= child->id();
    if(id          ==ParticleID::piplus)  ++npip;
    else if(id     ==ParticleID::piminus) ++npim;
    else if(abs(id)==ParticleID::Kplus)   ++nk;
  }
  if((id0==ParticleID::Dplus &&(nk==1&&npip==2))||
     (id0==ParticleID::Dminus&&(nk==1&&npim==2))) return 0;
  else return -1;
}

Complex CLEODptoKmPipPipDalitz::amplitude(int ichan) const {
  Complex output(0.);
  unsigned int imin=0, imax(nRes());
  if(ichan>=0) {
    imin=ichan;
    imax=imin+1;
  }
  for(unsigned int ix=imin;ix<imax;++ix) {
    if(ix==2 || ix==3 || ix==8 || ix==9) {
      output += resAmp(ix,true);
    }
    else
      output += resAmp(ix);
  }
  if(ichan<0) {
    output += 7.4*Complex(cos(-0.3211405823669566),sin(-0.3211405823669566));
  }
  return output;
}
