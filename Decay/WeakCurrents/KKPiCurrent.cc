// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the KKPiCurrent class.
//

#include "KKPiCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/epsilon.h"

using namespace Herwig;

KKPiCurrent::KKPiCurrent() {
  // masses for the isoscalar component
  isoScalarMasses_ = {1019.461*MeV,1630*MeV,1960*MeV};
  isoScalarWidths_ = {  4.249*MeV, 218*MeV, 267*MeV};
  // masses for the isovector component
  isoVectorMasses_ = {775.26*MeV,1465*MeV,1720*MeV};
  isoVectorWidths_ = {149.1 *MeV, 400*MeV, 250*MeV};
  // amplitude and phases for the isoscalar
  isoScalarKStarAmp_  ={ZERO, 0.233/GeV, 0.0405/GeV};
  isoScalarKStarPhase_={ 0.,  1.1E-07,  5.19};
  // amplitudes and phase for the isovector component
  isoVectorKStarAmp_  ={-2.34/GeV, 0.594/GeV, -0.0179/GeV};
  isoVectorKStarPhase_={0.,0.317, 2.57};
  // Coupling for the K* to Kpi
  gKStar_ = 5.37392360229;
  // mstar masses
  mKStarP_ = 895.6*MeV;
  mKStar0_ = 895.6*MeV;
  wKStarP_ = 47.0*MeV;
  wKStar0_ = 47.0*MeV;
  // modes
  addDecayMode(3,-3);
  addDecayMode(3,-3);
  addDecayMode(3,-3);
  addDecayMode(3,-3);
  addDecayMode(3,-3);
  addDecayMode(3,-3);
}

IBPtr KKPiCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr KKPiCurrent::fullclone() const {
  return new_ptr(*this);
}

void KKPiCurrent::doinit() {
  WeakCurrent::doinit();
  static const Complex ii(0.,1.);
  assert(isoScalarKStarAmp_.size()==isoScalarKStarPhase_.size());
  for(unsigned int ix=0;ix<isoScalarKStarAmp_.size();++ix) {
    isoScalarKStarCoup_.push_back(isoScalarKStarAmp_[ix]*(cos(isoScalarKStarPhase_[ix])
							  +ii*sin(isoScalarKStarPhase_[ix])));
  }
  assert(isoVectorKStarAmp_.size()==isoVectorKStarPhase_.size());
  for(unsigned int ix=0;ix<isoVectorKStarAmp_.size();++ix)
    isoVectorKStarCoup_.push_back(isoVectorKStarAmp_[ix]*(cos(isoVectorKStarPhase_[ix])
						+ii*sin(isoVectorKStarPhase_[ix])));
}

void KKPiCurrent::persistentOutput(PersistentOStream & os) const {
  os << ounit(isoScalarMasses_,GeV) << ounit(isoScalarWidths_,GeV)
     << ounit(isoVectorMasses_,GeV) << ounit(isoVectorWidths_,GeV)
     << ounit(isoScalarKStarAmp_,1./GeV) << ounit(isoVectorKStarAmp_,1./GeV)
     << isoScalarKStarPhase_ << isoVectorKStarPhase_
     << ounit(isoScalarKStarCoup_,1./GeV) << ounit(isoVectorKStarCoup_,1./GeV)
     << gKStar_
     << ounit(mKStarP_,GeV) <<  ounit(mKStar0_,GeV)
     << ounit(wKStarP_,GeV) << ounit(wKStar0_,GeV);
}

void KKPiCurrent::persistentInput(PersistentIStream & is, int) {
  is >> iunit(isoScalarMasses_,GeV) >> iunit(isoScalarWidths_,GeV)
     >> iunit(isoVectorMasses_,GeV) >> iunit(isoVectorWidths_,GeV)
     >> iunit(isoScalarKStarAmp_,1./GeV) >> iunit(isoVectorKStarAmp_,1./GeV)
     >> isoScalarKStarPhase_ >> isoVectorKStarPhase_
     >> iunit(isoScalarKStarCoup_,1./GeV) >> iunit(isoVectorKStarCoup_,1./GeV)
     >> gKStar_
     >> iunit(mKStarP_,GeV) >>  iunit(mKStar0_,GeV)
     >> iunit(wKStarP_,GeV) >> iunit(wKStar0_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<KKPiCurrent,WeakCurrent>
describeHerwigKKPiCurrent("Herwig::KKPiCurrent", "HwWeakCurrents.so");

void KKPiCurrent::Init() {

  static ClassDocumentation<KKPiCurrent> documentation
    ("The KKPiCurrent class implements a simple model for the production"
     " of K K pi in e+e- collisions via rho and phi resonances with a"
     " K*K intermediate state.");

  static ParVector<KKPiCurrent,Energy> interfaceIsoScalarMasses
    ("IsoScalarMasses",
     "The masses for I=0 part of the current",
     &KKPiCurrent::isoScalarMasses_, GeV, -1, 1.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);
  
  static ParVector<KKPiCurrent,Energy> interfaceIsoVectorMasses
    ("IsoVectorMasses",
     "The masses for I=1 part of the current",
     &KKPiCurrent::isoVectorMasses_, GeV, -1, 1.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<KKPiCurrent,Energy> interfaceIsoScalarWidths
    ("IsoScalarWidths",
     "The widths for I=0 part of the current",
     &KKPiCurrent::isoScalarWidths_, GeV, -1, 1.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);
  
  static ParVector<KKPiCurrent,Energy> interfaceIsoVectorWidths
    ("IsoVectorWidths",
     "The widths for I=1 part of the current",
     &KKPiCurrent::isoVectorWidths_, GeV, -1, 1.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<KKPiCurrent,InvEnergy> interfaceIsoScalarKStarAmp
    ("IsoScalarKStarAmp",
     "Amplitude for the I=0 K* component",
     &KKPiCurrent::isoScalarKStarAmp_, 1./GeV, -1, 0.0/GeV, -1000.0/GeV, 1000.0/GeV,
     false, false, Interface::limited);
  
  static ParVector<KKPiCurrent,InvEnergy> interfaceIsoVectorKStarAmp
    ("IsoVectorKStarAmp",
     "Amplitude for the I=1 K* component",
     &KKPiCurrent::isoVectorKStarAmp_, 1./GeV, -1, 0.0/GeV, -1000.0/GeV, 1000.0/GeV,
     false, false, Interface::limited);
  
  static ParVector<KKPiCurrent,double> interfaceIsoScalarKStarPhase
    ("IsoScalarKStarPhase",
     "The phase of the couplings for the I=0 part of the current",
     &KKPiCurrent::isoScalarKStarPhase_, -1, 0., -Constants::pi, Constants::pi,
     false, false, Interface::limited);
  
  static ParVector<KKPiCurrent,double> interfaceIsoVectorKStarPhase
    ("IsoVectorKStarPhase",
     "The phase of the couplings for the I=1 part of the current",
     &KKPiCurrent::isoVectorKStarPhase_, -1, 0., -Constants::pi, Constants::pi,
     false, false, Interface::limited);

  static Parameter<KKPiCurrent,Energy> interfacemKStarPlus
    ("mKStarPlus",
     "The mass of the charged K* resonace",
     &KKPiCurrent::mKStarP_, GeV, 0.8956*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);
  
  static Parameter<KKPiCurrent,Energy> interfacemKStar0
    ("mKStar0",
     "The mass of the neutral K* resonace",
     &KKPiCurrent::mKStar0_, GeV, 0.8956*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<KKPiCurrent,Energy> interfacewKStarPlus
    ("wKStarPlus",
     "The width of the charged K* resonance",
     &KKPiCurrent::wKStarP_, GeV, 0.047*GeV, 0.0*GeV, 1.0*GeV,
     false, false, Interface::limited);
  
  static Parameter<KKPiCurrent,Energy> interfacewKStar0
    ("wKStar0",
     "The width of the neutral K* resonance",
     &KKPiCurrent::wKStar0_, GeV, 0.047*GeV, 0.0*GeV, 1.0*GeV,
     false, false, Interface::limited);
  
  static Parameter<KKPiCurrent,double> interfacegKStar
    ("gKStar",
     "The coupling of K* K pi",
     &KKPiCurrent::gKStar_, 5.37392360229, 0.0, 10.0,
     false, false, Interface::limited);

}


// complete the construction of the decay mode for integration
bool KKPiCurrent::createMode(int icharge, tcPDPtr resonance,
			     FlavourInfo flavour,
			     unsigned int imode,PhaseSpaceModePtr mode,
			     unsigned int iloc,int ires,
			     PhaseSpaceChannel phase, Energy upp ) {
  // check the charge
  if(icharge!=0) return false;
  if(imode>5) return false;
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I==IsoSpin::IZero) {
      if(flavour.I3!=IsoSpin::I3Zero) return false;
    }
    else if(flavour.I==IsoSpin::IOne) {
      if(flavour.I3!=IsoSpin::I3Zero) return false;
    }
    else
      return false;
  }
  if(flavour.strange != Strangeness::Unknown)
     if(flavour.strange != Strangeness::Zero and flavour.strange != Strangeness::ssbar) return false;
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero       ) return false;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  != Beauty::Zero      ) return false;
  // get the external particles
  tPDVector out = particles(0,imode,0,0);
  // check the kinematics
  Energy mout(ZERO);
  for(unsigned int ix=0;ix<out.size();++ix)
    mout += out[ix]->mass();
  if(mout>upp) return false;
  // resonances we need
  vector<tPDPtr> resI1 = {getParticleData(   113),getParticleData(100113),getParticleData( 30113)};
  vector<tPDPtr> resI0 = {getParticleData(   333),getParticleData(100333)};
  tPDPtr res[2];
  if(imode==0) {
    res[0] = getParticleData(ParticleID::Kstar0);
    res[1] = getParticleData(ParticleID::Kstarbar0);
  }
  else if(imode==1) {
    res[0] = getParticleData(ParticleID::Kstarplus);
    res[1] = getParticleData(ParticleID::Kstarminus);
  }
  else if(imode==2||imode==4) {
    res[0] = getParticleData(ParticleID::Kstarplus);
    res[1] = getParticleData(ParticleID::Kstarbar0);
  }
  else if(imode==3||imode==5) {
    res[0] = getParticleData(ParticleID::Kstarminus);
    res[1] = getParticleData(ParticleID::Kstar0);
  }
  for(unsigned int ix=0;ix<resI0.size();++ix) {
    if(resonance && resonance != resI0[ix]) continue;
    mode->addChannel((PhaseSpaceChannel(phase),ires,resI0[ix],ires+1,res[0],ires+1,iloc+2,
		      ires+2,iloc+1,ires+2,iloc+3));
    mode->addChannel((PhaseSpaceChannel(phase),ires,resI0[ix],ires+1,res[1],ires+1,iloc+1,
		      ires+2,iloc+2,ires+2,iloc+3));
  }
  for(unsigned int ix=0;ix<resI1.size();++ix) {
    if(resonance && resonance != resI1[ix]) continue;
    mode->addChannel((PhaseSpaceChannel(phase),ires,resI1[ix],ires+1,res[0],ires+1,iloc+2,
		      ires+2,iloc+1,ires+2,iloc+3));
    mode->addChannel((PhaseSpaceChannel(phase),ires,resI1[ix],ires+1,res[1],ires+1,iloc+1,
		      ires+2,iloc+2,ires+2,iloc+3));
  }
  return true;
}

// the particles produced by the current
tPDVector KKPiCurrent::particles(int icharge, unsigned int imode,
				 int,int) {
  assert(icharge==0);
  if(imode==0)
    return {getParticleData(ParticleID::K_S0 ),getParticleData(ParticleID::K_L0  ),getParticleData(ParticleID::pi0)};
  else if(imode==1) 
    return {getParticleData(ParticleID::Kplus),getParticleData(ParticleID::Kminus),getParticleData(ParticleID::pi0)};
  else if(imode==2)
    return {getParticleData(ParticleID::K_S0 ),getParticleData(ParticleID::Kminus),getParticleData(ParticleID::piplus)};
  else if(imode==3)
    return {getParticleData(ParticleID::K_S0 ),getParticleData(ParticleID::Kplus ),getParticleData(ParticleID::piminus)};
  else if(imode==4)
    return {getParticleData(ParticleID::K_L0 ),getParticleData(ParticleID::Kminus),getParticleData(ParticleID::piplus)};
  else if(imode==5)
    return {getParticleData(ParticleID::K_L0 ),getParticleData(ParticleID::Kplus ),getParticleData(ParticleID::piminus)};
  else
    assert(false);
}


// hadronic current   
vector<LorentzPolarizationVectorE> 
KKPiCurrent::current(tcPDPtr resonance,
		     FlavourInfo flavour,
		     const int imode, const int ichan, Energy & scale, 
		     const tPDVector & ,
		     const vector<Lorentz5Momentum> & momenta,
		     DecayIntegrator::MEOption) const {
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I==IsoSpin::IZero) {
      if(flavour.I3!=IsoSpin::I3Zero) return vector<LorentzPolarizationVectorE>();
    }
    else if(flavour.I==IsoSpin::IOne) {
      if(flavour.I3!=IsoSpin::I3Zero) return vector<LorentzPolarizationVectorE>();
    }
    else
      return vector<LorentzPolarizationVectorE>();
  }
  int ssbar=0;
  if(flavour.strange != Strangeness::Unknown) {
    if(flavour.strange == Strangeness::Zero) ssbar=1;
    else if (flavour.strange == Strangeness::ssbar) ssbar=2;
    else assert(false);
  }
  useMe();
  // calculate q2,s1,s2
  Lorentz5Momentum q;
  for(unsigned int ix=0;ix<momenta.size();++ix) q+=momenta[ix];
  q.rescaleMass();
  scale=q.mass();
  Energy2 q2=q.mass2();
  Energy2 s1 = (momenta[0]+momenta[2]).m2();
  Energy2 s2 = (momenta[1]+momenta[2]).m2();
  // I=0 coefficient
  complex<InvEnergy> A0(ZERO);
  int ires=-1;
  if(ichan>=0) ires=ichan/2;
  if(resonance) {
    int ires2=-1;
    switch(abs(resonance->id())) {
    case 113:
      ires2=2;
      break;
    case 100113:
      ires2=3;
      break;
    case 30113:
      ires2=4;
      break;
    case 333:
      ires2=0;
      break;
    case 100333:
      ires2=1;
      break;
    };
    if(ires>=0 && ires!=ires2) {
      return vector<LorentzPolarizationVectorE>();
    }
    ires=ires2;
  }
  if((flavour.I==IsoSpin::IUnknown || flavour.I==IsoSpin::IZero) && ssbar!=1) {
    if(ires>=0) {
      if(ires<int(isoScalarMasses_.size()))
	A0 = isoScalarKStarCoup_[ires]*Resonance::BreitWignerFW(q2,isoScalarMasses_[ires],isoScalarWidths_[ires]);
    }
    else {
      for(unsigned int ix=0;ix<isoScalarMasses_.size();++ix) {
	A0 += isoScalarKStarCoup_[ix]*Resonance::BreitWignerFW(q2,isoScalarMasses_[ix],isoScalarWidths_[ix]);
      }
    }
  }
  ires-=2;
  // I=1 coefficient
  complex<InvEnergy> A1(ZERO);
  if((flavour.I==IsoSpin::IUnknown || flavour.I==IsoSpin::IOne) && ssbar!=2) {
    if(ires>=0) {
      if(ires<int(isoVectorMasses_.size()))
	A1  = isoVectorKStarCoup_[ires]*Resonance::BreitWignerFW(q2,isoVectorMasses_[ires],isoVectorWidths_[ires]);
    }
    else {
      for(unsigned int ix=0;ix<isoVectorMasses_.size();++ix) {
	A1 += isoVectorKStarCoup_[ix]*Resonance::BreitWignerFW(q2,isoVectorMasses_[ix],isoVectorWidths_[ix]);
      }
    }
  }
  complex<InvEnergy3> amp(ZERO);
  ires = -1;
  if(ichan>=0) ires = ichan%2;
  if(imode==0) {
    complex<InvEnergy2> r1 = (ires<0||ires==0) ?
     Resonance::BreitWignerPWave(s1,mKStar0_,wKStar0_,momenta[0].mass(),momenta[2].mass())/sqr(mKStar0_) : InvEnergy2(); 
    complex<InvEnergy2> r2 = (ires<0||ires==1) ?
     Resonance::BreitWignerPWave(s2,mKStar0_,wKStar0_,momenta[1].mass(),momenta[2].mass())/sqr(mKStar0_) : InvEnergy2();
    amp = sqrt(1./6.)*(A0+A1)*(r1+r2);
  }
  else if(imode==1) {
    complex<InvEnergy2> r1 = (ires<0||ires==0) ?
     Resonance::BreitWignerPWave(s1,mKStarP_,wKStarP_,momenta[0].mass(),momenta[2].mass())/sqr(mKStarP_) : InvEnergy2(); 
    complex<InvEnergy2> r2 = (ires<0||ires==1) ?
     Resonance::BreitWignerPWave(s2,mKStarP_,wKStarP_,momenta[1].mass(),momenta[2].mass())/sqr(mKStarP_) : InvEnergy2();
    amp = sqrt(1./6.)*(A0-A1)*(r1+r2);
  }
  else {
    complex<InvEnergy2> r1 = (ires<0||ires==0) ?
     Resonance::BreitWignerPWave(s1,mKStarP_,wKStarP_,momenta[0].mass(),momenta[2].mass())/sqr(mKStarP_) : InvEnergy2(); 
    complex<InvEnergy2> r2 = (ires<0||ires==1) ?
     Resonance::BreitWignerPWave(s2,mKStar0_,wKStar0_,momenta[1].mass(),momenta[2].mass())/sqr(mKStar0_) : InvEnergy2();
    amp = sqrt(1./6.)*((A0+A1)*r1+(A0-A1)*r2);
  }
  amp *= 2.*gKStar_;
  // the current
  LorentzPolarizationVector vect = amp*Helicity::epsilon(momenta[0],momenta[1],momenta[2]);
  // factor to get dimensions correct
  return vector<LorentzPolarizationVectorE>(1,scale*vect);
}
   
bool KKPiCurrent::accept(vector<int> id) {
  if(id.size()!=3) return false;
  int npip(0),npim(0),nkp(0),nkm(0),
    npi0(0),nks(0),nkl(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)       ++npip;
    else if(id[ix]==ParticleID::piminus) ++npim;
    else if(id[ix]==ParticleID::Kplus)   ++nkp;
    else if(id[ix]==ParticleID::Kminus)  ++nkm;
    else if(id[ix]==ParticleID::pi0)     ++npi0;
    else if(id[ix]==ParticleID::K_S0)    ++nks;
    else if(id[ix]==ParticleID::K_L0)    ++nkl;
  }
  if ( (npi0==1 && (( nks==1&&nkl==1 ) ||
		    ( nkp==1&&nkm==1 )) ) ||
       ( (nkl==1||nks==1) &&
	 ( (nkm==1&&npip==1) || (nkp==1&&npim==1) ) ) ) return true;
  return false;
}

// the decay mode
unsigned int KKPiCurrent::decayMode(vector<int> id) {
  assert(id.size()==3);
  int npip(0),npim(0),nkp(0),nkm(0),
    npi0(0),nks(0),nkl(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)       ++npip;
    else if(id[ix]==ParticleID::piminus) ++npim;
    else if(id[ix]==ParticleID::Kplus)   ++nkp;
    else if(id[ix]==ParticleID::Kminus)  ++nkm;
    else if(id[ix]==ParticleID::pi0)     ++npi0;
    else if(id[ix]==ParticleID::K_S0)    ++nks;
    else if(id[ix]==ParticleID::K_L0)    ++nkl;
  }
  if     ( nks==1&&nkl==1&&npi0==1 ) return 0;
  else if( nkp==1&&nkm==1&&npi0==1 ) return 1;
  else if( nks==1&&nkm==1&&npip==1 ) return 2;
  else if( nks==1&&nkp==1&&npim==1 ) return 3;
  else if( nkl==1&&nkm==1&&npip==1 ) return 4;
  else if( nkl==1&&nkp==1&&npim==1 ) return 5;
  assert(false);
}

// output the information for the database
void KKPiCurrent::dataBaseOutput(ofstream & output,bool header,
					bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::KKPiCurrent " 
  		    << name() << " HwWeakCurrents.so\n";
  for(unsigned int ix=0;ix<isoScalarMasses_.size();++ix) {
    if(ix<2) output << "newdef ";
    else     output << "insert ";
    output << name() << "IsoScalarMasses " << ix << " " << isoScalarMasses_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<isoVectorMasses_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << "IsoVectorMasses " << ix << " " << isoVectorMasses_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<isoScalarWidths_.size();++ix) {
    if(ix<2) output << "newdef ";
    else     output << "insert ";
    output << name() << "IsoScalarWidths " << ix << " " << isoScalarWidths_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<isoVectorWidths_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << "IsoVectorWidths " << ix << " " << isoVectorWidths_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<isoScalarKStarAmp_.size();++ix) {
    if(ix<2) output << "newdef ";
    else     output << "insert ";
    output << name() << "IsoScalarKStarAmp " << ix << " " << isoScalarKStarAmp_[ix]*GeV << "\n";
  }
  for(unsigned int ix=0;ix<isoVectorKStarAmp_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << "IsoVectorKStarAmp " << ix << " " << isoVectorKStarAmp_[ix]*GeV << "\n";
  }
  for(unsigned int ix=0;ix<isoScalarKStarPhase_.size();++ix) {
    if(ix<2) output << "newdef ";
    else     output << "insert ";
    output << name() << "IsoScalarKStarPhase " << ix << " " << isoScalarKStarPhase_[ix] << "\n";
  }
  for(unsigned int ix=0;ix<isoVectorKStarPhase_.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << "IsoVectorKStarPhase " << ix << " " << isoVectorKStarPhase_[ix] << "\n";
  }
  output << "newdef " << name() << ":mKStarPlus " << mKStarP_/GeV  << "\n";
  output << "newdef " << name() << ":mKStar0 "    << mKStar0_/GeV  << "\n";
  output << "newdef " << name() << ":wKStarPlus " << wKStarP_/GeV  << "\n";
  output << "newdef " << name() << ":wKStar0 "    << wKStar0_/GeV  << "\n";
  output << "newdef " << name() << ":gKStar "     << gKStar_ << "\n";
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
  		    << fullName() << "\";" << endl;
}
