// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CLEOD0toK0PipPim class.
//

#include "CLEOD0toK0PipPim.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

CLEOD0toK0PipPim::CLEOD0toK0PipPim() : WeakDalitzDecay(5./GeV,1.5/GeV,true) {
  // amplitudes and phases for D0 -> K0pi+pi-
  aKstarp_   = 0.11     ; phiKstarp_ = 321;
  arho_      = 1.00     ; phirho_    =   0;
  aomega_    = 0.037    ; phiomega_  = 114;
  aKstarm_   = 1.56     ; phiKstarm_ = 150;
  af980_     = 0.34*GeV2; phif980_   = 188;
  af2_       = 0.7/GeV2 ; phif2_     = 308;
  af1370_    = 1.8*GeV2 ; phif1370_  =  85;
  aK14300_   = 2.0*GeV2 ; phiK14300_ =   3;
  aK14302_   = 1.0/GeV2 ; phiK14302_ = 335;
  aK1680_    = 5.6      ; phiK1680_  = 174;
  aNR_       = 1.1      ; phiNR_     = 340;
  // masses and widths
  f0opt_ = true;
  gpi_   = 0.09;
  gK_    = 0.02;
  mK892_   =  891.66*MeV; wK892_   =  50.8 *MeV;
  mrho_    =  769.3 *MeV; wrho_    = 150.2 *MeV;
  momega_  =  782.57*MeV; womega_  =   8.44*MeV;   
  mf980_   =  977.00*MeV; wf980_   =  50.  *MeV;
  mf2_     = 1275.4 *MeV; wf2_     = 185.1 *MeV;
  mf1370_  = 1310   *MeV; wf1370_  = 272.0 *MeV;
  mK14300_ = 1412   *MeV; wK14300_ = 294   *MeV;
  mK14302_ = 1425.6 *MeV; wK14302_ =  98.5 *MeV;
  mK1680_  = 1717   *MeV; wK1680_  = 322   *MeV;
  // intermediates
  generateIntermediates(true);
}

IBPtr CLEOD0toK0PipPim::clone() const {
  return new_ptr(*this);
}

IBPtr CLEOD0toK0PipPim::fullclone() const {
  return new_ptr(*this);
}

void CLEOD0toK0PipPim::persistentOutput(PersistentOStream & os) const {
  os << ounit(mrho_,GeV)  << ounit(wrho_,GeV)  << ounit(momega_,GeV) << ounit(womega_,GeV)
     << ounit(mf980_,GeV) << ounit(wf980_,GeV) << gpi_ << gK_ << f0opt_
     << ounit(mf2_,GeV) << ounit(wf2_,GeV) << ounit(mf1370_,GeV) << ounit(wf1370_,GeV)
     << ounit(mK892_,GeV) << ounit(wK892_,GeV) << ounit(mK14300_,GeV) << ounit(wK14300_,GeV)
     << ounit(mK14302_,GeV) << ounit(wK14302_,GeV) << ounit(mK1680_,GeV) << ounit(wK1680_,GeV)
     << aKstarp_ << phiKstarp_ << arho_ << phirho_ << aomega_ << phiomega_
     << aKstarm_ <<  phiKstarm_ << ounit(af980_,GeV2) << phif980_ << ounit(af2_,1./GeV2) << phif2_
     << ounit(af1370_,GeV2) << phif1370_ << ounit(aK14300_,GeV2) << phiK14300_
     << ounit(aK14302_,1./GeV2) << phiK14302_
     << aK1680_ <<  phiK1680_ << aNR_ << phiNR_ << cNR_;
}

void CLEOD0toK0PipPim::persistentInput(PersistentIStream & is, int) {
  is >> iunit(mrho_,GeV)  >> iunit(wrho_,GeV)  >> iunit(momega_,GeV) >> iunit(womega_,GeV)
     >> iunit(mf980_,GeV) >> iunit(wf980_,GeV) >> gpi_ >> gK_ >> f0opt_
     >> iunit(mf2_,GeV) >> iunit(wf2_,GeV) >> iunit(mf1370_,GeV) >> iunit(wf1370_,GeV)
     >> iunit(mK892_,GeV) >> iunit(wK892_,GeV) >> iunit(mK14300_,GeV) >> iunit(wK14300_,GeV)
     >> iunit(mK14302_,GeV) >> iunit(wK14302_,GeV) >> iunit(mK1680_,GeV) >> iunit(wK1680_,GeV)
     >> aKstarp_ >> phiKstarp_ >> arho_ >> phirho_ >> aomega_ >> phiomega_
     >> aKstarm_ >>  phiKstarm_ >> iunit(af980_,GeV2) >> phif980_ >> iunit(af2_,1./GeV2) >> phif2_
     >> iunit(af1370_,GeV2) >> phif1370_ >> iunit(aK14300_,GeV2) >> phiK14300_
     >> iunit(aK14302_,1./GeV2) >> phiK14302_
     >> aK1680_ >>  phiK1680_ >> aNR_ >> phiNR_ >> cNR_;
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

  static Parameter<CLEOD0toK0PipPim,Energy> interfacerhoMass
    ("RhoMass",
     "The mass of the rho0 meson",
     &CLEOD0toK0PipPim::mrho_, MeV, 769.3 *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,Energy> interfacerhoWidth
    ("RhoWidth",
     "The width of the rho0 meson in D",
     &CLEOD0toK0PipPim::wrho_, MeV, 150.2 *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,Energy> interfaceOmegaMass
    ("OmegaMass",
     "The mass of the omega meson",
     &CLEOD0toK0PipPim::momega_, MeV, 782.57*MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,Energy> interfaceOmegaWidth
    ("OmegaWidth",
     "The width of the omega meson",
     &CLEOD0toK0PipPim::womega_, MeV, 8.44*MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,Energy> interfacef980Mass
    ("f980Mass",
     "The mass of the f_0(980) meson",
     &CLEOD0toK0PipPim::mf980_, MeV, 977.00*MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,Energy> interfacef980Width
    ("f980Width",
     "The width of the f_0(980) meson",
     &CLEOD0toK0PipPim::wf980_, MeV,  50.  *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,double> interfacegPi
    ("gPi",
     "The g_pi coupling for the f_0(980) width",
     &CLEOD0toK0PipPim::gpi_, 0.09, 0.0, 1.,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,double> interfacegK
    ("gK",
     "The g_K coupling for the f_0(980) width",
     &CLEOD0toK0PipPim::gK_, 0.02, 0.0, 1.,
     false, false, Interface::limited);

  static Switch<CLEOD0toK0PipPim,bool> interfacef0Option
    ("f0Option",
     "Option for the treatment of the f_0(980) width",
     &CLEOD0toK0PipPim::f0opt_, true, false, false);
  static SwitchOption interfacef0OptionCoupled
    (interfacef0Option,
     "Coupled",
     "Use the coupling pion and kaon channels",
     true);
  static SwitchOption interfacef0OptionSWave
    (interfacef0Option,
     "SWave",
     "Use a simple s-wave running width",
     false);

  static Parameter<CLEOD0toK0PipPim,Energy> interfacef_2Mass
    ("f_2Mass",
     "The mass of the f_2 meson",
     &CLEOD0toK0PipPim::mf2_, MeV, 1275.4 *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,Energy> interfacef_2Width
    ("f_2Width",
     "The width of the f_2 meson",
     &CLEOD0toK0PipPim::wf2_, MeV, 185.1 *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,Energy> interfacef1370Mass
    ("f1370Mass",
     "The mass of the f_0(1370) meson",
     &CLEOD0toK0PipPim::mf1370_, MeV, 1310   *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,Energy> interfacef1370Width
    ("f1370Width",
     "The width of the f_0(1370) meson",
     &CLEOD0toK0PipPim::wf1370_, MeV, 272.0 *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,Energy> interfaceKstarMass
    ("KstarMass",
     "The mass of the K*+(892) meson",
     &CLEOD0toK0PipPim::mK892_, MeV, 891.66*MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,Energy> interfaceKstarWidth
    ("KstarWidth",
     "The width of the K*+(892) meson",
     &CLEOD0toK0PipPim::wK892_, MeV, 50.8 *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);
  
  static Parameter<CLEOD0toK0PipPim,Energy> interfaceK_01430Mass
    ("K_01430Mass",
     "The mass of the K_0(1430) meson",
     &CLEOD0toK0PipPim::mK14300_, MeV, 1412   *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,Energy> interfaceK_01430Width
    ("K_01430Width",
     "The width of the K_0(1430) meson",
     &CLEOD0toK0PipPim::wK14300_, MeV, 294   *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,Energy> interfaceK_21430Mass
    ("K_21430Mass",
     "The mass of the K_2(1430) meson",
     &CLEOD0toK0PipPim::mK14302_, MeV, 1425.6 *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);
  static Parameter<CLEOD0toK0PipPim,Energy> interfaceK_21430Width
    ("K_21430Width",
     "The width of the K_2(1430) meson",
     &CLEOD0toK0PipPim::wK14302_, MeV,  98.5 *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,Energy> interfaceKstar1680Mass
    ("Kstar1680Mass",
     "The mass of the K*(1680) meson",
     &CLEOD0toK0PipPim::mK1680_, MeV, 1717   *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,Energy> interfaceKstar1680Width
    ("Kstar1680Width",
     "The width of the K*(1680) meson",
     &CLEOD0toK0PipPim::wK1680_, MeV, 322   *MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,double> interfaceKStarPlusAmplitude
    ("KStarPlusAmplitude",
     "Amplitude for the K*(892)+ component for D0 -> Kbar0 pi+ pi-",
     &CLEOD0toK0PipPim::aKstarp_, 0.11, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,double> interfaceKStarPlusPhase
    ("KStarPlusPhase",
     "Phase for the K*(892)+ component for D0 -> Kbar0 pi+ pi-",
     &CLEOD0toK0PipPim::phiKstarp_, 321., 0.0, 360.0,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,double> interfaceRhoAmplitude
    ("RhoAmplitude",
     "Amplitude for the rho0 component for D0 -> Kbar0 pi+ pi-",
     &CLEOD0toK0PipPim::arho_, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,double> interfaceRhoPhase
    ("RhoPhase",
     "Phase for the rho0 component for D0 -> Kbar0 pi+ pi-",
     &CLEOD0toK0PipPim::phirho_, 0.0, 0.0, 360.0,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,double> interfaceOmegaAmplitude
    ("OmegaAmplitude",
     "Amplitude for the omega component for D0 -> Kbar0 pi+ pi-",
     &CLEOD0toK0PipPim::aomega_, 0.037, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,double> interfaceOmegaPhase
    ("OmegaPhase",
     "Phase for the omega component for D0 -> Kbar0 pi+ pi-",
     &CLEOD0toK0PipPim::phiomega_, 114., 0.0, 360.0,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,double> interfaceKStarMinusAmplitude
    ("KStarMinusAmplitude",
     "Amplitude for the K*(892)- component for D0 -> Kbar0 pi+ pi-",
     &CLEOD0toK0PipPim::aKstarm_, 1.56, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,double> interfaceKStarMinusPhase
    ("KStarMinusPhase",
     "Phase for the K*(892)- component for D0 -> Kbar0 pi+ pi-",
     &CLEOD0toK0PipPim::phiKstarm_, 150., 0.0, 360.0,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,Energy2> interfacef980Amplitude
    ("f980Amplitude",
     "Amplitude for the f_0(980) component for D0 -> Kbar0 pi+ pi-",
     &CLEOD0toK0PipPim::af980_, GeV2, 0.34*GeV2, ZERO, 10.0*GeV2,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,double> interfacef980Phase
    ("f980Phase",
     "Phase for the f_0(980) component for D0 -> Kbar0 pi+ pi-",
     &CLEOD0toK0PipPim::phif980_, 188., 0.0, 360.0,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,InvEnergy2> interfacef2Amplitude
    ("f2Amplitude",
     "Amplitude for the f_2 component for D0 -> Kbar0 pi+ pi-",
     &CLEOD0toK0PipPim::af2_, 1./GeV2, 0.7/GeV2, ZERO, 10.0/GeV2,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,double> interfacef2Phase
    ("f2Phase",
     "Phase for the f_0(2) component for D0 -> Kbar0 pi+ pi-",
     &CLEOD0toK0PipPim::phif2_, 308., 0.0, 360.0,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,Energy2> interfacef1370Amplitude
    ("f1370Amplitude",
     "Amplitude for the f_0(1370) component for D0 -> Kbar0 pi+ pi-",
     &CLEOD0toK0PipPim::af1370_, GeV2, 1.8*GeV2, ZERO, 10.0*GeV2,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,double> interfacef1370Phase
    ("f1370Phase",
     "Phase for the f_0(1370) component for D0 -> Kbar0 pi+ pi-",
     &CLEOD0toK0PipPim::phif1370_, 85., 0.0, 360.0,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,Energy2> interfaceKK_0MinusAmplitude
    ("KK_0MinusAmplitude",
     "Amplitude for the K_0(1430)- component for D0 -> Kbar0 pi+ pi-",
     &CLEOD0toK0PipPim::aK14300_, GeV2, 2.0*GeV2, ZERO, 10.0*GeV2,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,double> interfaceKK_0MinusPhase
    ("KK_0MinusPhase",
     "Phase for the K_0(1430)- component for D0 -> Kbar0 pi+ pi-",
     &CLEOD0toK0PipPim::phiK14300_, 3, 0.0, 360.0,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,InvEnergy2> interfaceKK_2MinusAmplitude
    ("KK_2MinusAmplitude",
     "Amplitude for the K_2(1430)- component for D0 -> Kbar0 pi+ pi-",
     &CLEOD0toK0PipPim::aK14302_, 1./GeV2, 1.0/GeV2, ZERO, 10.0/GeV2,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,double> interfaceKK_2MinusPhase
    ("KK_2MinusPhase",
     "Phase for the K_2(1430)- component for D0 -> Kbar0 pi+ pi-",
     &CLEOD0toK0PipPim::phiK14302_, 335, 0.0, 360.0,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,double> interfaceK1680MinusAmplitude
    ("K1680MinusAmplitude",
     "Amplitude for the K*(892)- component for D0 -> Kbar0 pi+ pi-",
     &CLEOD0toK0PipPim::aK1680_, 5.6, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,double> interfaceK1680MinusPhase
    ("K1680MinusPhase",
     "Phase for the K*(892)- component for D0 -> Kbar0 pi+ pi-",
     &CLEOD0toK0PipPim::phiK1680_, 174, 0.0, 360.0,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,double> interfaceNonResonantAmplitude
    ("NonResonantAmplitude",
     "Amplitude for the non-resonant component for D0 -> Kbar0 pi+ pi-",
     &CLEOD0toK0PipPim::aNR_, 1.1, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<CLEOD0toK0PipPim,double> interfaceNonResonantPhase
    ("NonResonantPhase",
     "Phase for the non-resonant component for D0 -> Kbar0 pi+ pi-",
     &CLEOD0toK0PipPim::phiNR_, 340, 0.0, 360.0,
     false, false, Interface::limited);
}

void CLEOD0toK0PipPim::doinit() {
  WeakDalitzDecay::doinit();
  static const double degtorad = Constants::pi/180.;
  // non-resonant amplitude
  cNR_ = aNR_*Complex(cos(phiNR_*degtorad),sin(phiNR_*degtorad));
  // resonances
  addResonance(DalitzResonance(getParticleData(    323),mK892_  , wK892_  ,0,1,2,-aKstarp_     , phiKstarp_*degtorad));
  addResonance(DalitzResonance(getParticleData(    113),mrho_   , wrho_   ,1,2,0,-arho_        , phirho_   *degtorad));
  addResonance(DalitzResonance(getParticleData(    223),momega_ , womega_ ,1,2,0,-aomega_      , phiomega_ *degtorad));
  addResonance(DalitzResonance(getParticleData(   -323),mK892_  , wK892_  ,0,2,1,-aKstarm_     , phiKstarm_*degtorad));
  addResonance(DalitzResonance(getParticleData(9010221),mf980_  , wf980_  ,1,2,0, af980_/GeV2  , phif980_  *degtorad));
  addResonance(DalitzResonance(getParticleData(    225),mf2_    , wf2_    ,1,2,0, af2_*GeV2    , phif2_    *degtorad));
  addResonance(DalitzResonance(getParticleData(  10221),mf1370_ , wf1370_ ,1,2,0, af1370_/GeV2 , phif1370_ *degtorad));
  addResonance(DalitzResonance(getParticleData( -10321),mK14300_, wK14300_,0,2,1, aK14300_/GeV2, phiK14300_*degtorad));
  addResonance(DalitzResonance(getParticleData(   -325),mK14302_, wK14302_,0,2,1, aK14302_*GeV2, phiK14302_*degtorad));
  addResonance(DalitzResonance(getParticleData( -30323),mK1680_ , wK1680_ ,0,2,1,-aK1680_      , phiK1680_ *degtorad));
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
      if(!f0opt_) {
	output += resAmp(ix);
      }
      else {
	const Energy & mAB = mInv(resonances()[ix].daughter1,resonances()[ix].daughter2);
	Energy Gamma_pi = gpi_*sqrt(0.25*sqr(mAB)-sqr(mOut(resonances()[ix].daughter1)));
	Energy2 arg = 0.25*sqr(mAB)-sqr(mOut(resonances()[ix].spectator));
	complex<Energy> Gamma_K  = arg>=ZERO ? gK_*sqrt(arg) : gK_*ii*sqrt(-arg);
	output += resonances()[ix].amp*GeV2/
	  (sqr(resonances()[ix].mass)-sqr(mAB)-ii*resonances()[ix].mass*(Gamma_pi+Gamma_K));
      }
    }
  }
  if(ichan<0) output += cNR_;
  return output;
}

void CLEOD0toK0PipPim::dataBaseOutput(ofstream & output, bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  output << "newdef " << name() << ":Rho0Mass "           << mrho_/MeV    << "\n";
  output << "newdef " << name() << ":Rho0Width "          << wrho_/MeV    << "\n";
  output << "newdef " << name() << ":OmegaMass "          << momega_/MeV  << "\n";
  output << "newdef " << name() << ":OmegaWidth "         << womega_/MeV  << "\n";
  output << "newdef " << name() << ":f980Mass "           << mf980_/MeV   << "\n";
  output << "newdef " << name() << ":f980Width "          << wf980_/MeV   << "\n";
  output << "newdef " << name() << ":gPi "      << gpi_   << "\n";
  output << "newdef " << name() << ":gK "       << gK_    << "\n";
  output << "newdef " << name() << ":f0Option " << f0opt_ << "\n";
  output << "newdef " << name() << ":f_2Mass "     << mf2_/MeV     << "\n";
  output << "newdef " << name() << ":f_2Width "    << wf2_/MeV     << "\n";
  output << "newdef " << name() << ":f1370Mass "   << mf1370_/MeV  << "\n";
  output << "newdef " << name() << ":f1370Width "  << wf1370_/MeV  << "\n";
  output << "newdef " << name() << ":KstarMass "   << mK892_/MeV   << "\n";
  output << "newdef " << name() << ":KstarWidth "  << wK892_/MeV   << "\n";
  output << "newdef " << name() << ":K_01430Mass "        << mK14300_/MeV  << "\n";
  output << "newdef " << name() << ":K_01430Width "       << wK14300_/MeV  << "\n";
  output << "newdef " << name() << ":K_21430Mass "        << mK14302_/MeV  << "\n";
  output << "newdef " << name() << ":K_21430Width "       << wK14302_/MeV  << "\n";
  output << "newdef " << name() << ":Kstar1680Mass "      << mK1680_/MeV   << "\n";
  output << "newdef " << name() << ":Kstar1680Width "     << wK1680_/MeV   << "\n";
  output << "newdef " << name() << ":KStarPlusAmplitude "   << aKstarp_      << "\n";
  output << "newdef " << name() << ":KStarPlusPhase "       << phiKstarp_    << "\n";
  output << "newdef " << name() << ":RhoAmplitude "         << arho_         << "\n";
  output << "newdef " << name() << ":RhoPhase "             << phirho_       << "\n";
  output << "newdef " << name() << ":OmegaAmplitude "       << aomega_       << "\n";
  output << "newdef " << name() << ":OmegaPhase "           << phiomega_     << "\n";
  output << "newdef " << name() << ":KStarMinusAmplitude "  << aKstarm_      << "\n";
  output << "newdef " << name() << ":KStarMinusPhase "      << phiKstarm_    << "\n";
  output << "newdef " << name() << ":f980Amplitude "        << af980_/GeV2   << "\n";
  output << "newdef " << name() << ":f980Phase "            << phif980_      << "\n";
  output << "newdef " << name() << ":f2Amplitude "          << af2_*GeV2     << "\n";
  output << "newdef " << name() << ":f2Phase "              << phif2_        << "\n";
  output << "newdef " << name() << ":f1370Amplitude "       << af1370_/GeV2  << "\n";
  output << "newdef " << name() << ":f1370Phase "           << phif1370_     << "\n";
  output << "newdef " << name() << ":KK_0MinusAmplitude "   << aK14300_/GeV2 << "\n";
  output << "newdef " << name() << ":KK_0MinusPhase "       << phiK14300_    << "\n";
  output << "newdef " << name() << ":KK_2MinusAmplitude "   << aK14302_*GeV2 << "\n";
  output << "newdef " << name() << ":KK_2MinusPhase "       << phiK14302_    << "\n";
  output << "newdef " << name() << ":K1680MinusAmplitude "  << aK1680_       << "\n";
  output << "newdef " << name() << ":K1680MinusPhase "      << phiK1680_     << "\n";
  output << "newdef " << name() << ":NonResonantAmplitude " << aNR_          << "\n";
  output << "newdef " << name() << ":NonResonantPhase "     << phiNR_        << "\n";
  if(header) {
    output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
  }
}
