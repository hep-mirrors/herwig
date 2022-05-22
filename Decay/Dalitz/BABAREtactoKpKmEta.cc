// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BABAREtactoKpKmEta class.
//

#include "BABAREtactoKpKmEta.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

BABAREtactoKpKmEta::BABAREtactoKpKmEta() :
  WeakDalitzDecay(1.5/GeV,1.5/GeV,false),
  mf0_1500_(1505*MeV), wf0_1500_(109*MeV),
  mf0_1710_(1720*MeV), wf0_1710_(135*MeV),
  mK0_1430_(1438*MeV), wK0_1430_(210*MeV),
  mf0_2200_(2189*MeV), wf0_2200_(238*MeV),
  mK0_1950_(1945*MeV), wK0_1950_(201*MeV),
  mf2_1525_(1525*MeV), wf2_1525_(73*MeV),
  mf0_1370_(1350*MeV),wf0_1370_(265*MeV),
  mf0_980_(0.922*GeV),wf0_980_(0.24*GeV),
  af0_1500_(1.), phif0_1500_( 0. ),
  af0_1710_(1.), phif0_1710_( 2.2),
  aK0_1430_(1.), phiK0_1430_( 2.3),
  af0_2200_(1.), phif0_2200_( 2.1),
  aK0_1950_(1.), phiK0_1950_(-0.2),
  af2_1525_(1.), phif2_1525_( 1.0),
  af0_1370_(1.), phif0_1370_( 0.9),
  af0_980_ (1.), phif0_980_ (-0.3),
  aNR_     (1.), phiNR_     (-1.2),
  channel1_(-1), channel2_(-1)  {
  // intermediates
  generateIntermediates(true);
}

IBPtr BABAREtactoKpKmEta::clone() const {
  return new_ptr(*this);
}

IBPtr BABAREtactoKpKmEta::fullclone() const {
  return new_ptr(*this);
}

void BABAREtactoKpKmEta::persistentOutput(PersistentOStream & os) const {
  os << ounit(mf0_1500_,GeV) << ounit(wf0_1500_,GeV) << ounit(mf0_1710_,GeV) << ounit(wf0_1710_,GeV)
     << ounit(mK0_1430_,GeV) << ounit(wK0_1430_,GeV) << ounit(mf0_2200_,GeV) << ounit(wf0_2200_,GeV)
     << ounit(mK0_1950_,GeV) << ounit(wK0_1950_,GeV) << ounit(mf2_1525_,GeV) << ounit(wf2_1525_,GeV)
     << ounit(mf0_1370_,GeV) << ounit(wf0_1370_,GeV) << ounit(mf0_980_ ,GeV) << ounit(wf0_980_ ,GeV)
     << af0_1500_ << phif0_1500_ << af0_1710_ << phif0_1710_
     << aK0_1430_ << phiK0_1430_ << af0_2200_ << phif0_2200_
     << aK0_1950_ << phiK0_1950_ << af2_1525_ << phif2_1525_
     << af0_1370_ << phif0_1370_ << af0_980_  << phif0_980_
     << aNR_ << phiNR_ << channel1_ << channel2_;
}

void BABAREtactoKpKmEta::persistentInput(PersistentIStream & is, int) {
  is >> iunit(mf0_1500_,GeV) >> iunit(wf0_1500_,GeV) >> iunit(mf0_1710_,GeV) >> iunit(wf0_1710_,GeV)
     >> iunit(mK0_1430_,GeV) >> iunit(wK0_1430_,GeV) >> iunit(mf0_2200_,GeV) >> iunit(wf0_2200_,GeV)
     >> iunit(mK0_1950_,GeV) >> iunit(wK0_1950_,GeV) >> iunit(mf2_1525_,GeV) >> iunit(wf2_1525_,GeV)
     >> iunit(mf0_1370_,GeV) >> iunit(wf0_1370_,GeV) >> iunit(mf0_980_ ,GeV) >> iunit(wf0_980_ ,GeV)
     >> af0_1500_ >> phif0_1500_ >> af0_1710_ >> phif0_1710_
     >> aK0_1430_ >> phiK0_1430_ >> af0_2200_ >> phif0_2200_
     >> aK0_1950_ >> phiK0_1950_ >> af2_1525_ >> phif2_1525_
     >> af0_1370_ >> phif0_1370_ >> af0_980_  >> phif0_980_
     >> aNR_ >> phiNR_ >> channel1_ >> channel2_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<BABAREtactoKpKmEta,WeakDalitzDecay>
describeHerwigBABAREtactoKpKmEta("Herwig::BABAREtactoKpKmEta", "HwDalitzDecay.so");

void BABAREtactoKpKmEta::Init() {

  static ClassDocumentation<BABAREtactoKpKmEta> documentation
    ("There is no documentation for the BABAREtactoKpKmEta class");

  // /**
  //  *  Masses and widths of the resonances
  //  */
  // //@{
  // /**
  //  *  Mass of the $f_0(1500)^0$
  //  */
  // Energy mf0_1500_;
  
  // /**
  //  *  Width of the $f_0(1500)^0$
  //  */
  // Energy wf0_1500_;
  
  // /**
  //  *  Mass of the $f_0(1710)^0$
  //  */
  // Energy mf0_1710_;
  
  // /**
  //  *  Width of the $f_0(1710)^0$
  //  */
  // Energy wf0_1710_;
  
  // /**
  //  *  Mass of the $K^*_0(1430)^+$
  //  */
  // Energy mK0_1430_;
  
  // /**
  //  *  Width of the $K^*_0(1430)^+$
  //  */
  // Energy wK0_1430_;
  
  // /**
  //  *  Mass of the $f_0(2200)^0$
  //  */
  // Energy mf0_2200_;
  
  // /**
  //  *  Width of the $f_0(2200)^0$
  //  */
  // Energy wf0_2200_;
  
  // /**
  //  *  Mass of the $K^*_0(1950)^+$
  //  */
  // Energy mK0_1950_;
  
  // /**
  //  *  Width of the $K^*_0(1950)^+$
  //  */
  // Energy wK0_1950_;
  
  // /**
  //  *  Mass of the $f_2(1525)^0$
  //  */
  // Energy mf2_1525_;
  
  // /**
  //  *  Width of the $f_2(1525)^0$
  //  */
  // Energy wf2_1525_;
  
  // /**
  //  *  Mass of the $f_0(1370)^0$
  //  */
  // Energy mf0_1370_;
  
  // /**
  //  *  Width of the $f_0(1370)^0$
  //  */
  // Energy wf0_1370_;

  // /**
  //  *  Mass of the \f$f_0(980\f$
  //  */
  // Energy mf0_980_;

  // /**
  //  *  Width of the \f$f_0(980\f$
  //  */
  // Energy wf0_980_;
  // //@}

  // /**
  //  *  Magnitudes and phases of the amplitudes
  //  */
  // //@{
  // /**
  //  *  Amplitude of the $f_0(1500)^0$
  //  */
  // double af0_1500_;
  
  // /**
  //  *  Phase of the $f_0(1500)^0$
  //  */
  // double phif0_1500_;
  
  // /**
  //  *  Amplitude of the $f_0(1710)^0$
  //  */
  // double af0_1710_;
  
  // /**
  //  *  Phase of the $f_0(1710)^0$
  //  */
  // double phif0_1710_;
  
  // /**
  //  *  Amplitude of the $K^*_0(1430)^+$
  //  */
  // double aK0_1430_;
  
  // /**
  //  *  Phase of the $K^*_0(1430)^+$
  //  */
  // double phiK0_1430_;
  
  // /**
  //  *  Amplitude of the $f_0(2200)^0$
  //  */
  // double af0_2200_;
  
  // /**
  //  *  Phase of the $f_0(2200)^0$
  //  */
  // double phif0_2200_;
  
  // /**
  //  *  Amplitude of the $K^*_0(1950)^+$
  //  */
  // double aK0_1950_;
  
  // /**
  //  *  Phase of the $K^*_0(1950)^+$
  //  */
  // double phiK0_1950_;
  
  // /**
  //  *  Amplitude of the $f_2(1525)^0$
  //  */
  // double af2_1525_;
  
  // /**
  //  *  Phase of the $f_2(1525)^0$
  //  */
  // double phif2_1525_;
  
  // /**
  //  *  Amplitude of the $f_0(1370)^0$
  //  */
  // double af0_1370_;
  
  // /**
  //  *  Phase of the $f_0(1370)^0$
  //  */
  // double phif0_1370_;

  // /**
  //  *  Amplitude of the \f$f_0(980\f$
  //  */
  // double af0_980_;

  // /**
  //  *  Phase of the \f$f_0(980\f$
  //  */
  // double phif0_980_;

  // /**
  //  *  Non-resonant amplitude
  //  */
  // double aNR_;

  // /**
  //  *  Non-resonant phase
  //  */
  // double phiNR_;
  static Parameter<BABAREtactoKpKmEta,int> interfaceChannel1
    ("Channel1",
     "The first allowed channel, for debugging/calculation of fit fractions only",
     &BABAREtactoKpKmEta::channel1_, -1, -1, 5,
     false, false, Interface::limited);
  
  static Parameter<BABAREtactoKpKmEta,int> interfaceChannel2
    ("Channel2",
     "The first allowed channel, for debugging/calculation of fit fractions only",
     &BABAREtactoKpKmEta::channel2_, -1, -1, 5,
     false, false, Interface::limited);
}

void BABAREtactoKpKmEta::doinit() {
  WeakDalitzDecay::doinit();
  // resonances
  addResonance(DalitzResonance(9030221,ResonanceType::Spin0,mf0_1500_ , wf0_1500_ ,0,1,2, af0_1500_ , phif0_1500_));
  addResonance(DalitzResonance(  10331,ResonanceType::Spin0,mf0_1710_ , wf0_1710_ ,0,1,2, af0_1710_ , phif0_1710_));
  addResonance(DalitzResonance(  10321,ResonanceType::Spin0,mK0_1430_ , wK0_1430_ ,0,2,1, aK0_1430_ , phiK0_1430_));
  addResonance(DalitzResonance( -10321,ResonanceType::Spin0,mK0_1430_ , wK0_1430_ ,1,2,0, aK0_1430_ , phiK0_1430_));
  addResonance(DalitzResonance(    335,ResonanceType::Spin2,mf2_1525_ , wf2_1525_ ,0,1,2, af2_1525_ , phif2_1525_));
  addResonance(DalitzResonance(  10221,ResonanceType::Spin0,mf0_1370_ , wf0_1370_ ,0,1,2, af0_1370_ , phif0_1370_));
  addResonance(DalitzResonance(9010221,ResonanceType::Spin0,mf0_980_  , wf0_980_  ,0,1,2, af0_980_  , phif0_980_ ));
  // eta_c -> K+ K- eta
  createMode(getParticleData(ParticleID::eta_c),
  	     {getParticleData(ParticleID::Kplus),
  	      getParticleData(ParticleID::Kminus),
  	      getParticleData(ParticleID::eta)});
}

int BABAREtactoKpKmEta::modeNumber(bool & cc,tcPDPtr parent,
				 const tPDVector & children) const {
  int id0(parent->id());
  cc = false;
  // incoming particle must be Ds
  if(id0!=ParticleID::eta_c) return -1;
  // must be three decay products
  if(children.size()!=3) return -1;
  tPDVector::const_iterator pit = children.begin();
  unsigned int nEta(0),nKp(0),nKm(0);
  for( ;pit!=children.end();++pit) {
    id0=(**pit).id();
    if(     id0==ParticleID::eta)  ++nEta;
    else if(id0==ParticleID::Kplus)  ++nKp;
    else if(id0==ParticleID::Kminus) ++nKm;
  }
  if (nKp == 1 && nKm == 1 && nEta==1) return  0;
  else                                 return -1;
}

Complex BABAREtactoKpKmEta::amplitude(int ichan) const {
  Complex output(0.);
  static const Complex ii(0.,1.);
  int imin=0, imax(resonances().size());
  if(ichan>=0) {
    imin=ichan;
    imax=imin+1;
  }
  for(int ix=imin;ix<imax;++ix) {
    if(channel1_>=0) {
      if(ix!=channel1_ && ix!=channel2_) continue;
    }
    // all resonances bar f0(980)
    if(ix!=6) {
      output += resAmp(ix);
    }
    // special treatment for f0(980)
    else {
      const Energy & mAB = mInv(resonances()[ix].daughter1,resonances()[ix].daughter2);
      const Energy & mK  = mOut(resonances()[ix].daughter1);
      double rho = 2.*Kinematics::pstarTwoBodyDecay(mAB,mK,mK)/mAB;
      output += resonances()[ix].amp*GeV2/
  	(sqr(resonances()[ix].mass)-sqr(mAB)-ii*resonances()[ix].mass*resonances()[ix].width*rho);
    }
  }
  // return if in the phase allocation
  if(ichan>=0) return output;
  // now for the channels where we don't have the resonances



  if(channel1_<0 || channel1_==10 || channel2_==10)
    output += aNR_*exp(ii*phiNR_);
  // return the outout
  return output;
}

void BABAREtactoKpKmEta::dataBaseOutput(ofstream & output, bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  WeakDalitzDecay::dataBaseOutput(output,false);
  // output << "newdef " << name() << ":KStarMass"  <<  mKStar_/GeV << "\n";
  // output << "newdef " << name() << ":KStarWidth"  <<  wKStar_/GeV << "\n";
  // output << "newdef " << name() << ":PhiMass"  <<  mPhi_/GeV << "\n";
  // output << "newdef " << name() << ":PhiWidth"  <<  wPhi_/GeV << "\n";
  // output << "newdef " << name() << ":f0_980_Mass"  <<  mf0_980_/GeV << "\n";
  // output << "newdef " << name() << ":f0_980_Width"  <<  wf0_980_/GeV << "\n";
  // output << "newdef " << name() << ":K0Mass"  <<  mK0_/GeV << "\n";
  // output << "newdef " << name() << ":K0Width"  <<  wK0_/GeV << "\n";
  // output << "newdef " << name() << ":f0_1710_Mass"  <<  mf0_1710_/GeV << "\n";
  // output << "newdef " << name() << ":f0_1710_Width"  <<  wf0_1710_/GeV << "\n";
  // output << "newdef " << name() << ":f0_1370_Mass"  <<  mf0_1370_/GeV << "\n";
  // output << "newdef " << name() << ":f0_1370_Width"  <<  wf0_1370_/GeV << "\n";
  // output << "newdef " << name() << ":KStarAmplitude"  <<  aKStar_ << "\n";
  // output << "newdef " << name() << ":KStarPhase"  <<  phiKStar_ << "\n";
  // output << "newdef " << name() << ":PhiAmplitude"  <<  aPhi_ << "\n";
  // output << "newdef " << name() << ":PhiPhase"  <<  phiPhi_ << "\n";
  // output << "newdef " << name() << ":f0_980_Amplitude"  <<  af0_980_ << "\n";
  // output << "newdef " << name() << ":f0_980_Phase"  <<  phif0_980_ << "\n";
  // output << "newdef " << name() << ":K0Amplitude"  <<  aK0_ << "\n";
  // output << "newdef " << name() << ":K0Phase"  <<  phiK0_ << "\n";
  // output << "newdef " << name() << ":f0_1710_Amplitude"  <<  af0_1710_ << "\n";
  // output << "newdef " << name() << ":f0_1710_Phase"  <<  phif0_1710_ << "\n";
  // output << "newdef " << name() << ":f0_1370_Amplitude"  <<  af0_1370_ << "\n";
  // output << "newdef " << name() << ":f0_1370_Phase"  <<  phif0_1370_ << "\n";
  if(header) {
    output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
  }
}
