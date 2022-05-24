// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CLEOD0toKmPipPi0 class.
//

#include "CLEOD0toKmPipPi0.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

CLEOD0toKmPipPi0::CLEOD0toKmPipPi0() : ScalarTo3ScalarDalitz(5./GeV,1.5/GeV,true) {
  // masses and widths
  mrho_     =  770   *MeV; wrho_     = 150.7 *MeV;
  mK892_    =  891.5 *MeV; wK892_    =  50   *MeV;
  mK8920_   =  896.1 *MeV; wK8920_   =  50.5 *MeV;
  mK14300_  = 1412   *MeV; wK14300_  = 294   *MeV;
  mrho1700_ = 1700   *MeV; wrho1700_ = 240   *MeV;
  mK1680_   = 1717   *MeV; wK1680_   = 322   *MeV;
  // amplitudes and phases for D0 -> K-pi+pi0
  aNR_      = 1.75     ; phiNR_      =  31.2;
  arho_     = 1.00     ; phirho_     =   0. ;
  aKstarm_  = 0.44     ; phiKstarm_  = 163  ;
  aKstar0_  = 0.39     ; phiKstar0_  =  -0.2;
  aK1430m_  = 0.77*GeV2; phiK1430m_  =  55.5;
  aK14300_  = 0.85*GeV2; phiK14300_  = 166  ;
  arho1700_ = 2.50     ; phirho1700_ = 171  ;
  aK1680_   = 2.50     ; phiK1680_   = 103  ;  
  // intermediates
  generateIntermediates(true);
}

IBPtr CLEOD0toKmPipPi0::clone() const {
  return new_ptr(*this);
}

IBPtr CLEOD0toKmPipPi0::fullclone() const {
  return new_ptr(*this);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<CLEOD0toKmPipPi0,ScalarTo3ScalarDalitz>
describeHerwigCLEOD0toKmPipPi0("Herwig::CLEOD0toKmPipPi0", "HwDalitzDecay.so");

void CLEOD0toKmPipPi0::persistentOutput(PersistentOStream & os) const {
  os << ounit(mK14300_,GeV) << ounit(wK14300_,GeV) << ounit(mK1680_,GeV) << ounit(wK1680_,GeV)
     << ounit(mrho1700_,GeV) << ounit(wrho1700_,GeV) << ounit(mK8920_,GeV) << ounit(wK8920_,GeV)
     << ounit(mK892_,GeV) << ounit(wK892_,GeV) << ounit(mrho_,GeV) << ounit(wrho_,GeV)
     << aNR_ << phiNR_ << arho_ << phirho_ << aKstarm_ << phiKstarm_ << aKstar0_ << phiKstar0_
     <<  ounit(aK1430m_,GeV2) << phiK1430m_ << ounit(aK14300_,GeV2) << phiK14300_
     << arho1700_ << phirho1700_ << aK1680_ << phiK1680_ << cNR_;
}

void CLEOD0toKmPipPi0::persistentInput(PersistentIStream & is, int) {
  is >> iunit(mK14300_,GeV) >> iunit(wK14300_,GeV) >> iunit(mK1680_,GeV) >> iunit(wK1680_,GeV)
     >> iunit(mrho1700_,GeV) >> iunit(wrho1700_,GeV) >> iunit(mK8920_,GeV) >> iunit(wK8920_,GeV)
     >> iunit(mK892_,GeV) >> iunit(wK892_,GeV) >> iunit(mrho_,GeV) >> iunit(wrho_,GeV)
     >> aNR_ >> phiNR_ >> arho_ >> phirho_ >> aKstarm_ >> phiKstarm_ >> aKstar0_ >> phiKstar0_
     >>  iunit(aK1430m_,GeV2) >> phiK1430m_ >> iunit(aK14300_,GeV2) >> phiK14300_
     >> arho1700_ >> phirho1700_ >> aK1680_ >> phiK1680_ >> cNR_;
}

void CLEOD0toKmPipPi0::Init() {

  static ClassDocumentation<CLEOD0toKmPipPi0> documentation
    ("The CLEOD0toKmPipPi0 class implements the model of CLEO for "
     "D0 -> K- pi+ pi0, Phys. Rev. D63 (2001) 092001.",
     "The CLEO fit of \\cite{Kopp:2000gv} was"
     " used for the decay $D^0\\to K^-\\pi^+\\pi^0$.",
     "\\bibitem{Kopp:2000gv} S.~Kopp {\\it et al.}  [CLEO Collaboration], "
     "Phys.\\ Rev.\\  D {\\bf 63} (2001) 092001 [arXiv:hep-ex/0011065].");

  static Parameter<CLEOD0toKmPipPi0,Energy> interfaceK_01430Mass
    ("K_01430Mass",
     "The mass of the K_0(1430) meson",
     &CLEOD0toKmPipPi0::mK14300_, MeV, 1412   *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<CLEOD0toKmPipPi0,Energy> interfaceK_01430Width
    ("K_01430Width",
     "The width of the K_0(1430) meson",
     &CLEOD0toKmPipPi0::wK14300_, MeV, 294   *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);
  
  static Parameter<CLEOD0toKmPipPi0,Energy> interfaceKstar1680Mass
    ("Kstar1680Mass",
     "The mass of the K*(1680) meson",
     &CLEOD0toKmPipPi0::mK1680_, MeV, 1717   *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<CLEOD0toKmPipPi0,Energy> interfaceKstar1680Width
    ("Kstar1680Width",
     "The width of the K*(1680) meson",
     &CLEOD0toKmPipPi0::wK1680_, MeV, 322   *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);
  
  static Parameter<CLEOD0toKmPipPi0,Energy> interfacerho1700Mass
    ("rho1700Mass",
     "The mass of the rho(1700) meson",
     &CLEOD0toKmPipPi0::mrho1700_, MeV, 1700   *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);
  
  static Parameter<CLEOD0toKmPipPi0,Energy> interfacerho1700Width
    ("rho1700Width",
     "The width of the rho(1700) meson",
     &CLEOD0toKmPipPi0::wrho1700_, MeV, 240   *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<CLEOD0toKmPipPi0,Energy> interfaceKstar0892Mass
    ("Kstar0892Mass",
     "The mass of the K*0(892) meson",
     &CLEOD0toKmPipPi0::mK8920_, MeV, 896.1 *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<CLEOD0toKmPipPi0,Energy> interfaceKstar0892Width
    ("Kstar0892Width",
     "The width of the K*0(892) meson",
     &CLEOD0toKmPipPi0::wK8920_, MeV, 50.5 *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);
  
  static Parameter<CLEOD0toKmPipPi0,Energy> interfaceKstarPlus892Mass
    ("KstarPlus892Mass",
     "The mass of the K*+(892) meson",
     &CLEOD0toKmPipPi0::mK892_, MeV, 891.5 *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<CLEOD0toKmPipPi0,Energy> interfaceKstarPlus892Width
    ("KstarPlus892Width",
     "The width of the K*+(892) meson in D0 -> K-pi+pi0",
     &CLEOD0toKmPipPi0::wK892_, MeV,  50   *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);
  
  static Parameter<CLEOD0toKmPipPi0,Energy> interfacerhoMass
    ("RhoMass",
     "The mass of the rho+ meson",
     &CLEOD0toKmPipPi0::mrho_, MeV, 770   *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);
  
  static Parameter<CLEOD0toKmPipPi0,Energy> interfacerhoWidth
    ("RhoWidth",
     "The width of the rho+ meson",
     &CLEOD0toKmPipPi0::wrho_, MeV, 150.7 *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);
  
  static Parameter<CLEOD0toKmPipPi0,double> interfaceNonResonantAmplitude
    ("NonResonantAmplitude",
     "Amplitude for the non-resonant component for D0 -> K- pi+ pi0",
     &CLEOD0toKmPipPi0::aNR_, 1.75, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<CLEOD0toKmPipPi0,double> interfaceNonResonantPhase
    ("NonResonantPhase",
     "Phase for the non-resonant component for D0 -> K- pi+ pi0",
     &CLEOD0toKmPipPi0::phiNR_, 31.2, -180.0, 180.0,
     false, false, Interface::limited);

  static Parameter<CLEOD0toKmPipPi0,double> interfaceRhoAmplitude
    ("RhoAmplitude",
     "Amplitude for the rho+ component for D0 -> K- pi+ pi0",
     &CLEOD0toKmPipPi0::arho_, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<CLEOD0toKmPipPi0,double> interfaceRhoPhase
    ("RhoPhase",
     "Phase for the rho+ component for D0 -> K- pi+ pi0",
     &CLEOD0toKmPipPi0::phirho_, 0.0, -180.0, 180.0,
     false, false, Interface::limited);

  static Parameter<CLEOD0toKmPipPi0,double> interfaceKStarMinusAmplitude
    ("KStarMinusAmplitude",
     "Amplitude for the K*(892)- component for D0 -> K- pi+ pi0",
     &CLEOD0toKmPipPi0::aKstarm_, 0.44, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<CLEOD0toKmPipPi0,double> interfaceKStarMinusPhase
    ("KStarMinusPhase",
     "Phase for the K*(892)- component for D0 -> K- pi+ pi0",
     &CLEOD0toKmPipPi0::phiKstarm_, 163, -180.0, 180.0,
     false, false, Interface::limited);

  static Parameter<CLEOD0toKmPipPi0,double> interfaceKStar0Amplitude
    ("KStar0Amplitude",
     "Amplitude for the K*(892)0 component for D0 -> K- pi+ pi0",
     &CLEOD0toKmPipPi0::aKstar0_, 0.39, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<CLEOD0toKmPipPi0,double> interfaceKStar0Phase
    ("KStar0Phase",
     "Phase for the K*(892)0 component for D0 -> K- pi+ pi0",
     &CLEOD0toKmPipPi0::phiKstar0_, -0.2, -180.0, 180.0,
     false, false, Interface::limited);

  static Parameter<CLEOD0toKmPipPi0,Energy2> interfaceK_0MinusAmplitude
    ("K_0MinusAmplitude",
     "Amplitude for the K_0(1430)- component for D0 -> K- pi+ pi0",
     &CLEOD0toKmPipPi0::aK1430m_, GeV2, 0.77*GeV2, ZERO, 10.0*GeV2,
     false, false, Interface::limited);

  static Parameter<CLEOD0toKmPipPi0,double> interfaceK_0MinusPhase
    ("K_0MinusPhase",
     "Phase for the K_0(1430)- component for D0 -> K- pi+ pi0",
     &CLEOD0toKmPipPi0::phiK1430m_, 55.5, -180.0, 180.0,
     false, false, Interface::limited);

  static Parameter<CLEOD0toKmPipPi0,Energy2> interfaceK_00Amplitude
    ("K_00Amplitude",
     "Amplitude for the K_0(1430)0 component for D0 -> K- pi+ pi0",
     &CLEOD0toKmPipPi0::aK14300_, GeV2, 0.85*GeV2, ZERO, 10.0*GeV2,
     false, false, Interface::limited);

  static Parameter<CLEOD0toKmPipPi0,double> interfaceK_00Phase
    ("K_00Phase",
     "Phase for the K_0(1430)0 component for D0 -> K- pi+ pi0",
     &CLEOD0toKmPipPi0::phiK14300_, 166, -180.0, 180.0,
     false, false, Interface::limited);

  static Parameter<CLEOD0toKmPipPi0,double> interfaceRho1700Amplitude
    ("Rho1700Amplitude",
     "Amplitude for the rho1700+ component for D0 -> K- pi+ pi0",
     &CLEOD0toKmPipPi0::arho1700_, 2.5, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<CLEOD0toKmPipPi0,double> interfaceRho1700Phase
    ("Rho1700Phase",
     "Phase for the rho1700+ component for D0 -> K- pi+ pi0",
     &CLEOD0toKmPipPi0::phirho1700_, 171., -180.0, 180.0,
     false, false, Interface::limited);

  static Parameter<CLEOD0toKmPipPi0,double> interfaceK1680MinusAmplitude
    ("K1680MinusAmplitude",
     "Amplitude for the K*(1680)- component for D0 -> K- pi+ pi0",
     &CLEOD0toKmPipPi0::aK1680_, 2.5, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<CLEOD0toKmPipPi0,double> interfaceK1680MinusPhase
    ("K1680MinusPhase",
     "Phase for the K*(1680)- component for D0 -> K- pi+ pi0",
     &CLEOD0toKmPipPi0::phiK1680_, 103, -180.0, 180.0,
     false, false, Interface::limited);
}

void CLEOD0toKmPipPi0::doinit() {
  ScalarTo3ScalarDalitz::doinit();
  static const double degtorad = Constants::pi/180.;  
  // non-resonant amplitude
  cNR_ = aNR_*Complex(cos(phiNR_*degtorad),sin(phiNR_*degtorad));
  // create the resonances
  addResonance(DalitzResonance( 213  ,ResonanceType::Spin1, mrho_   , wrho_    ,1,2,0,-arho_        , phirho_    *degtorad));
  addResonance(DalitzResonance(-323  ,ResonanceType::Spin1, mK892_  , wK892_   ,0,2,1,-aKstarm_     , phiKstarm_ *degtorad));
  addResonance(DalitzResonance(-313  ,ResonanceType::Spin1,mK8920_  , wK8920_  ,0,1,2,-aKstar0_     , phiKstar0_ *degtorad));
  addResonance(DalitzResonance(-10321,ResonanceType::Spin0,mK14300_ , wK14300_ ,0,2,1, aK1430m_/GeV2, phiK1430m_ *degtorad));
  addResonance(DalitzResonance(-10311,ResonanceType::Spin0,mK14300_ , wK14300_ ,0,1,2, aK14300_/GeV2, phiK14300_ *degtorad));
  addResonance(DalitzResonance( 30213,ResonanceType::Spin1,mrho1700_, wrho1700_,1,2,0,-arho1700_    , phirho1700_*degtorad));
  addResonance(DalitzResonance(-30323,ResonanceType::Spin1,mK1680_  , wK1680_  ,0,2,1,-aK1680_      , phiK1680_  *degtorad));
  // D+ -> K- pi+ pi+
  createMode(getParticleData(ParticleID::D0),
	     {getParticleData(ParticleID::Kminus),
		 getParticleData(ParticleID::piplus),
		 getParticleData(ParticleID::pi0)});
}

void CLEOD0toKmPipPi0::doinitrun() {
  ScalarTo3ScalarDalitz::doinitrun();
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

void CLEOD0toKmPipPi0::dataBaseOutput(ofstream & output, bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  ScalarTo3ScalarDalitz::dataBaseOutput(output,false);
  output << "newdef " << name() << ":K_01430Mass "        << mK14300_/MeV  << "\n";
  output << "newdef " << name() << ":K_01430Width "       << wK14300_/MeV  << "\n";
  output << "newdef " << name() << ":Kstar1680Mass "      << mK1680_/MeV   << "\n";
  output << "newdef " << name() << ":Kstar1680Width "     << wK1680_/MeV   << "\n";
  output << "newdef " << name() << ":rho1700Mass "        << mrho1700_/MeV << "\n";
  output << "newdef " << name() << ":rho1700Width "       << wrho1700_/MeV << "\n";
  output << "newdef " << name() << ":Kstar0892Mass "      << mK8920_/MeV   << "\n";
  output << "newdef " << name() << ":Kstar0892Width "     << wK8920_/MeV   << "\n";
  output << "newdef " << name() << ":KstarPlus892Mass "   << mK892_/MeV    << "\n";
  output << "newdef " << name() << ":KstarPlus892Width "  << wK892_/MeV    << "\n";
  output << "newdef " << name() << ":RhoMass "            << mrho_/MeV     << "\n";
  output << "newdef " << name() << ":RhoWidth "           << wrho_/MeV     << "\n";
  output << "newdef " << name() << ":NonResonantAmplitude " << aNR_          << "\n";
  output << "newdef " << name() << ":NonResonantPhase "     << phiNR_        << "\n";
  output << "newdef " << name() << ":RhoAmplitude "         << arho_         << "\n";
  output << "newdef " << name() << ":RhoPhase "             << phirho_       << "\n";
  output << "newdef " << name() << ":KStarMinusAmplitude "  << aKstarm_      << "\n";
  output << "newdef " << name() << ":KStarMinusPhase "      << phiKstarm_    << "\n";
  output << "newdef " << name() << ":KStar0Amplitude "      << aKstar0_      << "\n";
  output << "newdef " << name() << ":KStar0Phase "          << phiKstar0_    << "\n";
  output << "newdef " << name() << ":K_0MinusAmplitude "    << aK1430m_/GeV2 << "\n";
  output << "newdef " << name() << ":K_0MinusPhase "        << phiK1430m_    << "\n";
  output << "newdef " << name() << ":K_00Amplitude "        << aK14300_/GeV2 << "\n";
  output << "newdef " << name() << ":K_00Phase "            << phiK14300_    << "\n";
  output << "newdef " << name() << ":Rho1700Amplitude "     << arho1700_     << "\n";
  output << "newdef " << name() << ":Rho1700Phase "         << phirho1700_   << "\n";
  output << "newdef " << name() << ":K1680MinusAmplitude "  << aK1680_       << "\n";
  output << "newdef " << name() << ":K1680MinusPhase "      << phiK1680_     << "\n";
  if(header) {
    output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
  }
}
