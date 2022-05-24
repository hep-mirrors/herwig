// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BABARDstoKpKmPip class.
//

#include "BABARDstoKpKmPip.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Utilities/Kinematics.h"

using namespace Herwig;

BABARDstoKpKmPip::BABARDstoKpKmPip() : ScalarTo3ScalarDalitz(3./GeV,1.5/GeV,false),
				       mKStar_(895.6*MeV), wKStar_(45.1*MeV),
				       mPhi_(1019.455*MeV), wPhi_(4.26*MeV),
				       mf0_980_(0.922*GeV),wf0_980_(0.24*GeV),
				       mK0_(1425*MeV), wK0_(270*MeV),
				       mf0_1710_(1720*MeV),wf0_1710_(135*MeV),
				       mf0_1370_(1.22*GeV),wf0_1370_(0.21*GeV),
				       aKStar_(1.),phiKStar_(0.),
				       aPhi_(1.15),phiPhi_(2.89),
				       af0_980_(2.67),phif0_980_(1.56),
				       aK0_(1.14), phiK0_(2.55),
				       af0_1710_(0.65), phif0_1710_(1.36),
				       af0_1370_(0.46), phif0_1370_(-0.45) {
  // intermediates
  generateIntermediates(true);
}

IBPtr BABARDstoKpKmPip::clone() const {
  return new_ptr(*this);
}

IBPtr BABARDstoKpKmPip::fullclone() const {
  return new_ptr(*this);
}

void BABARDstoKpKmPip::persistentOutput(PersistentOStream & os) const {
  os << ounit(mKStar_,GeV) << ounit(wKStar_,GeV) << ounit(mPhi_,GeV) << ounit(wPhi_,GeV)
     << ounit(mf0_980_,GeV) << ounit(wf0_980_,GeV) << ounit(mK0_,GeV) << ounit(wK0_,GeV)
     << ounit(mf0_1710_,GeV) << ounit(wf0_1710_,GeV) << ounit(mf0_1370_,GeV) << ounit(wf0_1370_,GeV) 
     << aKStar_ << phiKStar_ << aPhi_ << phiPhi_
     << af0_980_ << phif0_980_ << aK0_ << phiK0_
     << af0_1710_ << phif0_1710_ << af0_1370_ << phif0_1370_;
}

void BABARDstoKpKmPip::persistentInput(PersistentIStream & is, int) {
  is >> iunit(mKStar_,GeV) >> iunit(wKStar_,GeV) >> iunit(mPhi_,GeV) >> iunit(wPhi_,GeV)
     >> iunit(mf0_980_,GeV) >> iunit(wf0_980_,GeV) >> iunit(mK0_,GeV) >> iunit(wK0_,GeV)
     >> iunit(mf0_1710_,GeV) >> iunit(wf0_1710_,GeV) >> iunit(mf0_1370_,GeV) >> iunit(wf0_1370_,GeV) 
     >> aKStar_ >> phiKStar_ >> aPhi_ >> phiPhi_
     >> af0_980_ >> phif0_980_ >> aK0_ >> phiK0_
     >> af0_1710_ >> phif0_1710_ >> af0_1370_ >> phif0_1370_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<BABARDstoKpKmPip,ScalarTo3ScalarDalitz>
describeHerwigBABARDstoKpKmPip("Herwig::BABARDstoKpKmPip", "HwDalitzDecay.so");

void BABARDstoKpKmPip::Init() {

  static ClassDocumentation<BABARDstoKpKmPip> documentation
    ("The BABARDstoKpKmPip class implements the BaBar Dalitz plot analysis"
     " for D_s+ -> K+ K- pi+",
     "The BABARDstoKpKmPip class implementing the BaBar Dalitz plot analysis"
     " from \\cite{BaBar:2010wqe} for D_s+ -> K+ K- pi+",
     "\\bibitem{BaBar:2010wqe} P.~del Amo Sanchez \\textit{et al.} [BaBar],\n"
     "%``Dalitz plot analysis of $D_s^+ \\to K^+ K^- \\pi^+$,''\n"
     "Phys. Rev. D \\textbf{83} (2011), 052001\n"
     "doi:10.1103/PhysRevD.83.052001\n"
     "[arXiv:1011.4190 [hep-ex]].\n"
     "%70 citations counted in INSPIRE as of 20 May 2022\n");

  static Parameter<BABARDstoKpKmPip,Energy> interfaceKStarMass
    ("KStarMass",
    "The mass of the K* resonance",
    &BABARDstoKpKmPip::mKStar_, GeV, 895.6*MeV, 0.0*GeV, 10.0*GeV,
    false, false, Interface::limited);

  static Parameter<BABARDstoKpKmPip,Energy> interfaceKStarWidth
    ("KStarWidth",
     "The mass of the K* resonance",
     &BABARDstoKpKmPip::wKStar_, GeV, 45.1*MeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<BABARDstoKpKmPip,Energy> interfacePhiMass
    ("PhiMass",
    "The mass of the phi resonance",
    &BABARDstoKpKmPip::mPhi_, GeV, 1019.455*MeV, 0.0*GeV, 10.0*GeV,
    false, false, Interface::limited);

  static Parameter<BABARDstoKpKmPip,Energy> interfacePhiWidth
    ("PhiWidth",
    "The mass of the phi resonance",
    &BABARDstoKpKmPip::wPhi_, GeV, 4.26*MeV, 0.0*GeV, 10.0*GeV,
    false, false, Interface::limited);

  static Parameter<BABARDstoKpKmPip,Energy> interfacef0_980_Mass
    ("f0_980_Mass",
    "The mass of the f0(980) resonance",
    &BABARDstoKpKmPip::mf0_980_, GeV, 0.922*GeV, 0.0*GeV, 10.0*GeV,
    false, false, Interface::limited);

  static Parameter<BABARDstoKpKmPip,Energy> interfacef0_980_Width
    ("f0_980_Width",
    "The mass of the f0(980) resonance",
    &BABARDstoKpKmPip::wf0_980_, GeV, 0.24*GeV, 0.0*GeV, 10.0*GeV,
    false, false, Interface::limited);

  static Parameter<BABARDstoKpKmPip,Energy> interfaceK0Mass
    ("K0Mass",
    "The mass of the K_0* resonance",
    &BABARDstoKpKmPip::mK0_, GeV, 1425*MeV, 0.0*GeV, 10.0*GeV,
    false, false, Interface::limited);

  static Parameter<BABARDstoKpKmPip,Energy> interfaceK0Width
    ("K0Width",
     "The mass of the K_0* resonance",
     &BABARDstoKpKmPip::wK0_, GeV, 270*MeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<BABARDstoKpKmPip,Energy> interfacef0_1710_Mass
    ("f0_1710_Mass",
    "The mass of the f0(1710) resonance",
    &BABARDstoKpKmPip::mf0_1710_, GeV, 1720*MeV, 0.0*GeV, 10.0*GeV,
    false, false, Interface::limited);

  static Parameter<BABARDstoKpKmPip,Energy> interfacef0_1710_Width
    ("f0_1710_Width",
    "The mass of the f0(1710) resonance",
    &BABARDstoKpKmPip::wf0_1710_, GeV, 135*MeV, 0.0*GeV, 10.0*GeV,
    false, false, Interface::limited);
  
  static Parameter<BABARDstoKpKmPip,Energy> interfacef0_1370_Mass
    ("f0_1370_Mass",
    "The mass of the f0(1370) resonance",
    &BABARDstoKpKmPip::mf0_1370_, GeV, 1.22*GeV, 0.0*GeV, 10.0*GeV,
    false, false, Interface::limited);
  
  static Parameter<BABARDstoKpKmPip,Energy> interfacef0_1370_Width
    ("f0_1370_Width",
    "The mass of the f0(1370) resonance",
    &BABARDstoKpKmPip::wf0_1370_, GeV, 0.21*GeV, 0.0*GeV, 10.0*GeV,
    false, false, Interface::limited);

    static Parameter<BABARDstoKpKmPip,double> interfaceKStarAmplitude
    ("KStarAmplitude",
     "The amplitude for the K* channel",
     &BABARDstoKpKmPip::aKStar_, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

    static Parameter<BABARDstoKpKmPip,double> interfaceKStarPhase
    ("KStarPhase",
     "The phase for the K* channel in radians",
    &BABARDstoKpKmPip::phiKStar_, 0.0, -Constants::pi, Constants::pi,
     false, false, Interface::limited);

    static Parameter<BABARDstoKpKmPip,double> interfacePhiAmplitude
    ("PhiAmplitude",
     "The amplitude for the phi channel",
     &BABARDstoKpKmPip::aPhi_, 1.15, 0.0, 10.0,
     false, false, Interface::limited);

    static Parameter<BABARDstoKpKmPip,double> interfacePhiPhase
    ("PhiPhase",
     "The phase for the phi channel in radians",
    &BABARDstoKpKmPip::phiPhi_, 2.89, -Constants::pi, Constants::pi,
     false, false, Interface::limited);

    static Parameter<BABARDstoKpKmPip,double> interfacef0_980_Amplitude
    ("f0_980_Amplitude",
     "The amplitude for the phi channel",
     &BABARDstoKpKmPip::af0_980_, 2.67, 0.0, 10.0,
     false, false, Interface::limited);

    static Parameter<BABARDstoKpKmPip,double> interfacef0_980_Phase
    ("f0_980_Phase",
     "The phase for the phi channel in radians",
    &BABARDstoKpKmPip::phif0_980_, 1.56, -Constants::pi, Constants::pi,
     false, false, Interface::limited);

    static Parameter<BABARDstoKpKmPip,double> interfaceK0Amplitude
    ("K0Amplitude",
     "The amplitude for the K_0* channel",
     &BABARDstoKpKmPip::aK0_, 1.14, 0.0, 10.0,
     false, false, Interface::limited);

    static Parameter<BABARDstoKpKmPip,double> interfaceK0Phase
    ("K0Phase",
     "The phase for the K_0* channel in radians",
    &BABARDstoKpKmPip::phiK0_, 2.55, -Constants::pi, Constants::pi,
     false, false, Interface::limited);

    static Parameter<BABARDstoKpKmPip,double> interfacef0_1710_Amplitude
    ("f0_1710_Amplitude",
     "The amplitude for the phi channel",
     &BABARDstoKpKmPip::af0_1710_, 0.65, 0.0, 10.0,
     false, false, Interface::limited);

    static Parameter<BABARDstoKpKmPip,double> interfacef0_1710_Phase
    ("f0_1710_Phase",
     "The phase for the phi channel in radians",
    &BABARDstoKpKmPip::phif0_1710_, 1.36, -Constants::pi, Constants::pi,
     false, false, Interface::limited);

    static Parameter<BABARDstoKpKmPip,double> interfacef0_1370_Amplitude
    ("f0_1370_Amplitude",
     "The amplitude for the phi channel",
     &BABARDstoKpKmPip::af0_1370_, 0.46, 0.0, 10.0,
     false, false, Interface::limited);

    static Parameter<BABARDstoKpKmPip,double> interfacef0_1370_Phase
    ("f0_1370_Phase",
     "The phase for the phi channel in radians",
    &BABARDstoKpKmPip::phif0_1370_, -0.45, -Constants::pi, Constants::pi,
     false, false, Interface::limited);
}

void BABARDstoKpKmPip::doinit() {
  ScalarTo3ScalarDalitz::doinit();
  // resonances
  addResonance(DalitzResonance(   -313,ResonanceType::Spin1  ,mKStar_   , wKStar_   ,1,2,0, aKStar_   , phiKStar_ ));
  // N.B. - sign w.r.t. paper to get right interference sign
  addResonance(DalitzResonance(    333,ResonanceType::Spin1  ,mPhi_     , wPhi_     ,0,1,2,-aPhi_     , phiPhi_   ));
  addResonance(DalitzResonance(9010221,ResonanceType::BABARf0,mf0_980_  , wf0_980_  ,0,1,2, af0_980_  , phif0_980_));
  addResonance(DalitzResonance( -10311,ResonanceType::Spin0  ,mK0_      , wK0_      ,1,2,0, aK0_      , phiK0_    ));
  addResonance(DalitzResonance(  10331,ResonanceType::Spin0  ,mf0_1710_ , wf0_1710_ ,0,1,2, af0_1710_ , phif0_1710_));
  addResonance(DalitzResonance(  10221,ResonanceType::Spin0  ,mf0_1370_ , wf0_1370_ ,0,1,2, af0_1370_ , phif0_1370_));
  // Ds+ -> K+ K- pi+
  createMode(getParticleData(ParticleID::D_splus),
  	     {getParticleData(ParticleID::Kplus),
	      getParticleData(ParticleID::Kminus),
	      getParticleData(ParticleID::piplus)});
}

int BABARDstoKpKmPip::modeNumber(bool & cc,tcPDPtr parent,
				 const tPDVector & children) const {
  int id0(parent->id());
  // incoming particle must be Ds
  if(abs(id0)!=ParticleID::D_splus) return -1;
  cc = id0==ParticleID::D_sminus;
  // must be three decay products
  if(children.size()!=3) return -1;
  tPDVector::const_iterator pit = children.begin();
  unsigned int npip(0),npim(0),nKp(0),nKm(0);
  for( ;pit!=children.end();++pit) {
    id0=(**pit).id();
    if(     id0==ParticleID::piplus)  ++npip;
    else if(id0==ParticleID::piminus) ++npim;
    else if(id0==ParticleID::Kplus)  ++nKp;
    else if(id0==ParticleID::Kminus) ++nKm;
  }
  if (nKp == 1 && nKm == 1) {
    if ( (!cc && npip==1) ||
	 ( cc && npim==1)) return 0;
    else
      return -1;
  }
  else
    return -1;
}

void BABARDstoKpKmPip::dataBaseOutput(ofstream & output, bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  ScalarTo3ScalarDalitz::dataBaseOutput(output,false);
  output << "newdef " << name() << ":KStarMass"  <<  mKStar_/GeV << "\n";
  output << "newdef " << name() << ":KStarWidth"  <<  wKStar_/GeV << "\n";
  output << "newdef " << name() << ":PhiMass"  <<  mPhi_/GeV << "\n";
  output << "newdef " << name() << ":PhiWidth"  <<  wPhi_/GeV << "\n";
  output << "newdef " << name() << ":f0_980_Mass"  <<  mf0_980_/GeV << "\n";
  output << "newdef " << name() << ":f0_980_Width"  <<  wf0_980_/GeV << "\n";
  output << "newdef " << name() << ":K0Mass"  <<  mK0_/GeV << "\n";
  output << "newdef " << name() << ":K0Width"  <<  wK0_/GeV << "\n";
  output << "newdef " << name() << ":f0_1710_Mass"  <<  mf0_1710_/GeV << "\n";
  output << "newdef " << name() << ":f0_1710_Width"  <<  wf0_1710_/GeV << "\n";
  output << "newdef " << name() << ":f0_1370_Mass"  <<  mf0_1370_/GeV << "\n";
  output << "newdef " << name() << ":f0_1370_Width"  <<  wf0_1370_/GeV << "\n";
  output << "newdef " << name() << ":KStarAmplitude"  <<  aKStar_ << "\n";
  output << "newdef " << name() << ":KStarPhase"  <<  phiKStar_ << "\n";
  output << "newdef " << name() << ":PhiAmplitude"  <<  aPhi_ << "\n";
  output << "newdef " << name() << ":PhiPhase"  <<  phiPhi_ << "\n";
  output << "newdef " << name() << ":f0_980_Amplitude"  <<  af0_980_ << "\n";
  output << "newdef " << name() << ":f0_980_Phase"  <<  phif0_980_ << "\n";
  output << "newdef " << name() << ":K0Amplitude"  <<  aK0_ << "\n";
  output << "newdef " << name() << ":K0Phase"  <<  phiK0_ << "\n";
  output << "newdef " << name() << ":f0_1710_Amplitude"  <<  af0_1710_ << "\n";
  output << "newdef " << name() << ":f0_1710_Phase"  <<  phif0_1710_ << "\n";
  output << "newdef " << name() << ":f0_1370_Amplitude"  <<  af0_1370_ << "\n";
  output << "newdef " << name() << ":f0_1370_Phase"  <<  phif0_1370_ << "\n";
  if(header) {
    output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
  }
}
