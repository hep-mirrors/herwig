// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ScalarTo3ScalarDalitz class.
//

#include "ScalarTo3ScalarDalitz.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;

IBPtr ScalarTo3ScalarDalitz::clone() const {
  return new_ptr(*this);
}

IBPtr ScalarTo3ScalarDalitz::fullclone() const {
  return new_ptr(*this);
}

void ScalarTo3ScalarDalitz::persistentOutput(PersistentOStream & os) const {
  os << f0gpi_ << f0gK_ << useResonanceMass_ ;
}

void ScalarTo3ScalarDalitz::persistentInput(PersistentIStream & is, int) {
  is >> f0gpi_ >> f0gK_ >> useResonanceMass_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ScalarTo3ScalarDalitz,DalitzBase>
describeHerwigScalarTo3ScalarDalitz("Herwig::ScalarTo3ScalarDalitz", "HwDalitzDecay.so");

void ScalarTo3ScalarDalitz::Init() {

  static ClassDocumentation<ScalarTo3ScalarDalitz> documentation
    ("The ScalarTo3ScalarDalitz class provides a base class for "
     "weak three-body decays of bottom and charm mesons");
  
  static Parameter<ScalarTo3ScalarDalitz,double> interfacegPi
    ("f0gPi",
     "The g_pi coupling for the f_0(980) width",
     &ScalarTo3ScalarDalitz::f0gpi_, 0.09, 0.0, 10.,
     false, false, Interface::limited);

  static Parameter<ScalarTo3ScalarDalitz,double> interfacegK
    ("f0gK",
     "The g_K coupling for the f_0(980) width",
     &ScalarTo3ScalarDalitz::f0gK_, 0.02, 0.0, 10.,
     false, false, Interface::limited);

  static Switch<ScalarTo3ScalarDalitz,bool> interfaceResonanceMass
    ("ResonanceMass",
     "Whether to use the kinematic mass or the resonance pole mass for the denominator in kinematic expressions",
     &ScalarTo3ScalarDalitz::useResonanceMass_, false, false, false);
  static SwitchOption interfaceResonanceMassYes
    (interfaceResonanceMass,
     "Yes",
     "Use the resonance mass, to be avoided only use if do in experimental fit",
     true);
  static SwitchOption interfaceResonanceMassNo
    (interfaceResonanceMass,
     "No",
     "Use the correct kinematic mass",
     false);

}

void ScalarTo3ScalarDalitz::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // set up the spin information for the decay products
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  for(unsigned int ix=0;ix<3;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
}

double ScalarTo3ScalarDalitz::me2(const int ichan, const Particle & part,
			    const tPDVector & ,
			    const vector<Lorentz5Momentum> & momenta,
			    MEOption meopt) const {
  static const Complex ii(0.,1.);
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin0,PDT::Spin0,PDT::Spin0)));
  useMe();
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&part),incoming);
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
  // now compute the matrix element
  Complex amp = amplitude(ichan);
  (*ME())(0,0,0,0) = amp;
  return norm(amp);
}

Complex ScalarTo3ScalarDalitz::resAmp(unsigned int i) const {
  Complex output = resonances()[i]->amp;
  if (resonances()[i]->type==ResonanceType::NonResonant) return output;
  // mass of the resonance
  const Energy & mR = resonances()[i]->mass ;
  // locations of the outgoing particles
  const unsigned int &d1 = resonances()[i]->daughter1;
  const unsigned int &d2 = resonances()[i]->daughter2;
  const unsigned int &sp = resonances()[i]->spectator;
  // compute the Breit-Wigner times resonance form factor piece
  output *= resonances()[i]->BreitWigner(m2_[d1][d2],mOut_[d1],mOut_[d2]);
  // spin zero is done now
  if(resonances()[i]->type%10==1 && resonances()[i]->type != ResonanceType::Spin0Gauss) {
    return output;
  }
  // spin piece x Blatt-Weisskopf for parent 
  else {
    double fD=1;
    // for the D decay
    Energy pD  = sqrt(max(ZERO,(0.25*sqr(sqr(mD_)-sqr(mR)-sqr(mOut_[sp])) - sqr(mR*mOut_[sp]))/sqr(mD_)));
    Energy pDAB= sqrt( 0.25*sqr(sqr(mD_)-sqr(m2_[d1][d2])-sqr(mOut_[sp])) - sqr(m2_[d1][d2]*mOut_[sp]))/mD_;
    double r2A(parentRadius()   *pD),r2B(parentRadius()   *pDAB);
    // mass for thre denominator
    Energy mDen = useResonanceMass_ ? resonances()[i]->mass : m2_[d1][d2];
    // denominator for the older form of the amplitude
    Energy2 denom = GeV2;
    if (resonances()[i]->type/10 == 1 ) { 
      Energy2 pa2 = 0.25*(sqr(m2_[d1][d2])-2.*(sqr(mOut_[d1])+sqr(mOut_[d2])) + sqr(sqr(mOut_[d1])-sqr(mOut_[d2]))/sqr(m2_[d1][d2]));
      Energy2 pc2 = 0.25*(sqr(m2_[d1][d2])-2.*(sqr(mD_      )+sqr(mOut_[sp])) + sqr(sqr(mD_      )-sqr(mOut_[sp]))/sqr(m2_[d1][d2]));
      denom = 4.*sqrt(pa2*pc2);
    }
    // Blatt-Weisskopf factors and spin piece
    switch (resonances()[i]->type) {
    case ResonanceType::Spin0Gauss:
      fD = exp(-(r2B-r2A)/12.);
      output *= fD;
      break;
    case ResonanceType::Spin1: case ResonanceType::Spin1E691 :
      fD=sqrt( (1. + sqr(r2A)) / (1. + sqr(r2B)) );
      output *= fD*(sqr(m2_[d2][sp])-sqr(m2_[d1][sp])
		    + (  sqr(mD_)-sqr(mOut_[sp]))*(sqr(mOut_[d1])-sqr(mOut_[d2]))/sqr(mDen) )/denom;
      break;
    case ResonanceType::Spin2: case ResonanceType::Spin2E691:
      fD = sqrt( (9. + sqr(r2A)*(3.+sqr(r2A))) / (9. + sqr(r2B)*(3.+sqr(r2B))));
      output *= fD/sqr(denom)*( sqr( sqr(m2_[d2][sp]) - sqr(m2_[d1][sp]) +(sqr(mD_)-sqr(mOut_[sp]))*(sqr(mOut_[d1])-sqr(mOut_[d2]))/(sqr(mDen))) -
				(sqr(m2_[d1][d2])-2*      sqr(mD_)-2*sqr(mOut_[sp]) + sqr((sqr(      mD_) - sqr(mOut_[sp]))/mDen))*
				(sqr(m2_[d1][d2])-2*sqr(mOut_[d1])-2*sqr(mOut_[d2]) + sqr((sqr(mOut_[d1]) - sqr(mOut_[d2]))/mDen))/3.);
      break;
    default :
      assert(false);
    }
    return output;
  }
}

void ScalarTo3ScalarDalitz::dataBaseOutput(ofstream & output, bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DalitzBase base class
  DalitzBase::dataBaseOutput(output,false);
  output << "newdef " << name() << ":ResonanceMass " << useResonanceMass_ << "\n";
  output << "newdef " << name() << ":f0gPi "      << f0gpi_   << "\n";
  output << "newdef " << name() << ":f0gK "       << f0gK_    << "\n";
  output << "\n";
  if(header) {
    output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
  }
}
