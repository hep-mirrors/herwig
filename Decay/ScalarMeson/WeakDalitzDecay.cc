// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the WeakDalitzDecay class.
//

#include "WeakDalitzDecay.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

WeakDalitzDecay::WeakDalitzDecay(InvEnergy rP, InvEnergy rR, bool useResonanceMass)
  : rParent_(rP), rResonance_(rR), useResonanceMass_(useResonanceMass), maxWgt_(1.)
{}

void WeakDalitzDecay::persistentOutput(PersistentOStream & os) const {
  os << resonances_ << maxWgt_ << weights_
     << ounit(rParent_,1./GeV) << ounit(rResonance_,1./GeV);
}

void WeakDalitzDecay::persistentInput(PersistentIStream & is, int) {
  is >> resonances_ >> maxWgt_ >> weights_
     >> iunit(rParent_,1./GeV) >> iunit(rResonance_,1./GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<WeakDalitzDecay,DecayIntegrator>
describeHerwigWeakDalitzDecay("Herwig::WeakDalitzDecay", "HwSMDecay.so");

void WeakDalitzDecay::Init() {

  static ClassDocumentation<WeakDalitzDecay> documentation
    ("The WeakDalitzDecay class provides a base class for "
     "weak three-body decays of bottom and charm mesons");

  static Parameter<WeakDalitzDecay,InvEnergy> interfaceDRadius
    ("DRadius",
     "The radius parameter for the Blatt-Weisskopf form-factor for the D",
     &WeakDalitzDecay::rParent_, 1./GeV, 5./GeV, ZERO, 10./GeV,
     false, false, Interface::limited);

  static Parameter<WeakDalitzDecay,InvEnergy> interfaceResonanceRadius
    ("ResonanceRadius",
     "The radius parameter for the Blatt-Weisskopf form-factor for the"
     "intermediate resonances",
     &WeakDalitzDecay::rResonance_, 1./GeV, 1.5/GeV, ZERO, 10./GeV,
     false, false, Interface::limited);

  static Parameter<WeakDalitzDecay,double> interfaceMaximumWeight
    ("MaximumWeight",
     "The maximum weight for the phase-space sampling",
     &WeakDalitzDecay::maxWgt_, 1.0, 0.0, 1e10,
     false, false, Interface::limited);

  static ParVector<WeakDalitzDecay,double> interfaceWeights
    ("Weights",
     "The weights for the different channels for the phase-space integration",
     &WeakDalitzDecay::weights_, -1, 1.0, 0.0, 1.0,
     false, false, Interface::limited);
}

void WeakDalitzDecay::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // set up the spin information for the decay products
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  for(unsigned int ix=0;ix<3;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
}

void WeakDalitzDecay::doinit() {
  DecayIntegrator::doinit();
}

void WeakDalitzDecay::doinitrun() {
  DecayIntegrator::doinitrun();
  weights_.resize(mode(0)->channels().size());
  maxWgt_ = mode(0)->maxWeight();
  for(unsigned int iz=0;iz<mode(0)->channels().size();++iz) {
    weights_[iz]=mode(0)->channels()[iz].weight();
  }
}

double WeakDalitzDecay::me2(const int ichan, const Particle & part,
			 const tPDVector & ,
			 const vector<Lorentz5Momentum> & momenta,
			 MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin0,PDT::Spin0,PDT::Spin0)));
  useMe();
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&part),incoming);
  }
  // set the kinematics
  mD_ = part.mass();
  for(unsigned int ix=0;ix<momenta.size();++ix) {
    mOut_[ix]=momenta[ix].mass();
    for(unsigned int iy=ix+1;iy<momenta.size();++iy) {
      m2_[ix][iy]=(momenta[ix]+momenta[iy]).m();
      m2_[iy][ix]=m2_[ix][iy];
    }
  }
  // compute the amplitude (in inherited class)
  Complex amp = amplitude(ichan);
  // now compute the matrix element
  (*ME())(0,0,0,0)=amp;
  return norm(amp);
}

void WeakDalitzDecay::createMode(tPDPtr in, tPDVector out) {
  PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(in,out,maxWgt_));
  if(weights_.size()!=resonances_.size()) {
    weights_=vector<double>(resonances_.size(),1./double(resonances_.size()));
  }
  unsigned int ix=0;
  for(DalitzResonance res : resonances_) {
    mode->addChannel((PhaseSpaceChannel(mode),0,res.resonance,0,res.spectator+1,1,res.daughter1+1,1,res.daughter2+1));
    resetIntermediate(res.resonance,res.mass,res.width);
    ++ix;
  }
  addMode(mode);
}

Complex WeakDalitzDecay::resAmp(unsigned int i, bool gauss) const {
  Complex output = resonances_[i].amp;
  // locations of the outgoing particles
  const unsigned int &d1 = resonances_[i].daughter1;
  const unsigned int &d2 = resonances_[i].daughter2;
  const unsigned int &sp = resonances_[i].spectator;
  // mass and width of the resonance
  const Energy & mR = resonances_[i].mass ;
  const Energy & wR = resonances_[i].width;
  // momenta for the resonance decay
  // off-shell
  Energy pAB=sqrt(0.25*sqr(sqr(m2_[d1][d2]) -sqr(mOut_[d1])-sqr(mOut_[d2])) - sqr(mOut_[d1]*mOut_[d2]))/m2_[d1][d2];
  //  on-shell
  Energy  pR=sqrt(0.25*sqr(    mR*mR        -sqr(mOut_[d1])-sqr(mOut_[d2])) - sqr(mOut_[d1]*mOut_[d2]))/mR;
  // for the D decay
  Energy pD  = sqrt(max(ZERO,(0.25*sqr(sqr(mD_)-sqr(mR)-sqr(mOut_[sp])) - sqr(mR*mOut_[sp]))/sqr(mD_)));
  Energy pDAB= sqrt( 0.25*sqr(sqr(mD_)-sqr(m2_[d1][d2])-sqr(mOut_[sp])) - sqr(m2_[d1][d2]*mOut_[sp]))/mD_;
  // Blatt-Weisskopf factors
  double fR=1, fD=1;
  unsigned int power(1);
  double r1A(rResonance_*pR),r1B(rResonance_*pAB );
  double r2A(rParent_   *pD),r2B(rParent_   *pDAB);
  // mass for thre denominator
  Energy mDen = useResonanceMass_ ? mR : m2_[d1][d2];
  // Blatt-Weisskopf factors and spin piece
  switch (resonances_[i].resonance->iSpin()) {
  case PDT::Spin0:
    if(gauss) {
      fR = exp(-(r1B-r1A)/12.);
      fD = exp(-(r2B-r2A)/12.);
    }
    break;
  case PDT::Spin1:
    fR=sqrt( (1. + sqr(r1A)) / (1. + sqr(r1B)) );
    fD=sqrt( (1. + sqr(r2A)) / (1. + sqr(r2B)) );
    power=3;
    output *= fR*fD*(sqr(m2_[d2][sp])-sqr(m2_[d1][sp])
		  + (  sqr(mD_)-sqr(mOut_[sp]))*(sqr(mOut_[d1])-sqr(mOut_[d2]))/sqr(mDen) )/GeV2;
    break;
  case PDT::Spin2:
    fR = sqrt( (9. + sqr(r1A)*(3.+sqr(r1A))) / (9. + sqr(r1B)*(3.+sqr(r1B))));
    fD = sqrt( (9. + sqr(r2A)*(3.+sqr(r2A))) / (9. + sqr(r2B)*(3.+sqr(r2B))));
    power=5;
    output *= fR*fD/GeV2/GeV2*( sqr( sqr(m2_[d2][sp]) - sqr(m2_[d1][sp]) +(sqr(mD_)-sqr(mOut_[sp]))*(sqr(mOut_[d1])-sqr(mOut_[d2]))/(sqr(mDen))) -
				(sqr(m2_[d1][d2])-2*      sqr(mD_)-2*sqr(mOut_[sp]) + sqr((sqr(      mD_) - sqr(mOut_[sp]))/mDen))*
				(sqr(m2_[d1][d2])-2*sqr(mOut_[d1])-2*sqr(mOut_[d2]) + sqr((sqr(mOut_[d1]) - sqr(mOut_[d2]))/mDen))/3.);
    break;
  default :
    assert(false);
  }
  // multiply by Breit-Wigner piece and return
  Energy gam = wR*pow(pAB/pR,power)*(mR/m2_[d1][d2])*fR*fR;
  return output*GeV2/(sqr(mR)-sqr(m2_[d1][d2])-mR*gam*Complex(0.,1.));
}

void WeakDalitzDecay::dataBaseOutput(ofstream & output, bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  output << "newdef " << name() << ":DRadius " << rParent_*GeV << "\n";
  output << "newdef " << name() << ":ResonanceRadius " << rResonance_*GeV << "\n";
  output << "newdef " << name() << ":MaximumWeight " << maxWgt_ << "\n";
  for(unsigned int ix=0;ix<weights_.size();++ix) {
    output << "insert " << name() << ":Weights "
	   << ix << " " << weights_[ix] << "\n";
  }
  if(header) {
    output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
  }
}
