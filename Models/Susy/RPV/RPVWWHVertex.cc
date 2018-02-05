// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RPVWWHVertex class.
//

#include "RPVWWHVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "RPV.h"

using namespace Herwig;

RPVWWHVertex::RPVWWHVertex() : coupLast_(ZERO), q2Last_(ZERO) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::SINGLET);
}

void RPVWWHVertex::doinit() {
  // extract some models parameters to decide if we need sneutrinos
  tRPVPtr model = dynamic_ptr_cast<tRPVPtr>(generator()->standardModel());
  if( !model ) throw InitException() << "RPVWWHVertex::doinit() - "
				     << "The pointer to the RPV object is null!"
				     << Exception::abortnow;
  // get the Higgs mixing matrix
  MixingMatrixPtr mix = model->CPevenHiggsMix();
  // possible interactions
  vector<long> higgs(2);
  higgs[0] = 25; higgs[1] = 35;
  if(mix->size().first>2) {
    higgs.push_back(1000012);
    higgs.push_back(1000014);
    higgs.push_back(1000016);
  }
  for(unsigned int ix=0;ix<higgs.size();++ix) {
    addToList( 23, 23,higgs[ix]);
    addToList(-24, 24,higgs[ix]);
  }
  VVSVertex::doinit();
  // SM parameters
  Energy mw = getParticleData(ParticleID::Wplus)->mass();
  Energy mz = getParticleData(ParticleID::Z0)->mass();
  double sw = sqrt(sin2ThetaW());
  double cw = sqrt(1.-sin2ThetaW());
  vector<Energy> vnu = model->sneutrinoVEVs();
  Energy v = 2.*mw/electroMagneticCoupling(sqr(mw))*sw;
  double tanb = model->tanBeta();
  Energy vd = sqrt((sqr(v)-sqr(vnu[0])-sqr(vnu[1])-sqr(vnu[2]))/
		   (1.+sqr(tanb)));
  Energy vu = vd*tanb;
  for(unsigned int ix=0;ix<higgs.size();++ix) {
    complex<Energy> c = vd*(*mix)(ix,0)+vu*(*mix)(ix,1);
    for(size_t iy=2; iy<mix->size().second; ++iy) c += vnu[iy-2]*(*mix)(ix,iy);
    vector<complex<Energy> > coup(2);
    coup[0] = c/v*mw;
    coup[1] = c/v*mz/cw;
    couplings_.push_back(coup);
  }
}

IBPtr RPVWWHVertex::clone() const {
  return new_ptr(*this);
}

IBPtr RPVWWHVertex::fullclone() const {
  return new_ptr(*this);
}

void RPVWWHVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(couplings_,GeV);
}

void RPVWWHVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(couplings_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<RPVWWHVertex,Helicity::VVSVertex>
describeHerwigRPVWWHVertex("Herwig::RPVWWHVertex", "HwSusy.so HwRPV.so");

void RPVWWHVertex::Init() {

  static ClassDocumentation<RPVWWHVertex> documentation
    ("The RPVWWHVertex class implements the couplings of a pair of electroweak"
     " gauge bosons to the higgs boson in he R-parity violating MSSM.");

}

void RPVWWHVertex::setCoupling(Energy2 q2, tcPDPtr particle1, tcPDPtr,
			      tcPDPtr particle3) {
  long bosonID = abs(particle1->id());
  long higgsID =     particle3->id();
  assert( bosonID == ParticleID::Wplus || bosonID == ParticleID::Z0 );
  int ihiggs = higgsID>1000000 ? (higgsID-1000008)/2 : (higgsID-25)/10;
  assert(ihiggs>=0 && ihiggs<=4);
  complex<Energy> fact = bosonID==ParticleID::Wplus ? 
    couplings_[ihiggs][0] : couplings_[ihiggs][1];
  if( q2 != q2Last_ ) {
    q2Last_ = q2;
    coupLast_ = weakCoupling(q2);
  }
  norm(coupLast_*fact*UnitRemoval::InvE);
}
