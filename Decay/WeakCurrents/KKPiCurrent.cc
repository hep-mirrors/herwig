// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the KKPiCurrent class.
//

#include "KKPiCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
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
  //  isoScalarMasses_ = {782.65*MeV,1019.461*MeV,1425*MeV,1680*MeV,1625*MeV,2188*MeV};
  //  isoScalarWidths_ = {  8.49*MeV,   4.249*MeV, 215*MeV, 150*MeV, 315*MeV,  83*MeV};
  isoScalarMasses_ = {1019.461*MeV,1630*MeV,1960*MeV};
  isoScalarWidths_ = {  4.249*MeV, 218*MeV, 267*MeV};
  // masses for the isovector component
  isoVectorMasses_ = {775.26*MeV,1465*MeV,1720*MeV};
  isoVectorWidths_ = {149.1 *MeV, 400*MeV, 250*MeV};
  // amplitude and phases for the isoscalar
  //  isoScalarKStarAmp_  ={ZERO,ZERO,ZERO,0.096/GeV,ZERO,ZERO};
  //  isoScalarKStarPhase_={  0.,  0.,  0.,     0.,  0.,  0.};
  isoScalarKStarAmp_  ={ZERO, 0.233/GeV, 0.0405/GeV};
  isoScalarKStarPhase_={ 0.,  1.1E-07,  5.19};
  // amplitudes and phase for the isovector component
  //  isoVectorKStarAmp_  ={ZERO,ZERO,0.04/GeV};
  //  isoVectorKStarPhase_={0.,0.,Constants::pi};
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
    ("There is no documentation for the KKPiCurrent class");


}


// complete the construction of the decay mode for integration
bool KKPiCurrent::createMode(int icharge, tcPDPtr resonance,
			     IsoSpin::IsoSpin Itotal, IsoSpin::I3 i3, Strangeness::Strange S,
			     unsigned int imode,PhaseSpaceModePtr mode,
			     unsigned int iloc,int ires,
			     PhaseSpaceChannel phase, Energy upp ) {
  // check the charge
  if(icharge!=0) return false;
  if(imode>5) return false;
  // check the total isospin
  //   if(Itotal!=IsoSpin::IUnknown) {
  //     if(Itotal==IsoSpin::IZero) {
  //       if(i3!=IsoSpin::I3Unknown) return false;
  //     }
  //     else if(Itotal==IsoSpin::IOne) {
  //       if(i3!=IsoSpin::I3Unknown&&
  // 	 i3!=IsoSpin::I3One) return false;
  //     }
  //     else
  //       return false;
  //   }
  
  
  
  // get the external particles
  tPDVector out = particles(0,imode,0,0);
  // check the kinematics
  Energy mout(ZERO);
  for(unsigned int ix=0;ix<out.size();++ix)
    mout += out[ix]->mass();
  if(mout>upp) return false;
  // resonances we need
  tPDPtr resI1[3] = {getParticleData(   113),getParticleData(100113),getParticleData( 30113)};
  tPDPtr resI0[6] = {getParticleData(   223),getParticleData(   333),
		     getParticleData(100223),getParticleData(100333),
		     getParticleData( 30223)};
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
  for(unsigned int ix=0;ix<5;++ix) {
    mode->addChannel((PhaseSpaceChannel(phase),ires,resI0[ix],ires+1,res[0],ires+1,iloc+2,
		      ires+2,iloc+1,ires+2,iloc+3));
    mode->addChannel((PhaseSpaceChannel(phase),ires,resI0[ix],ires+1,res[1],ires+1,iloc+1,
		      ires+2,iloc+2,ires+2,iloc+3));
  }
  for(unsigned int ix=0;ix<3;++ix) {
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
		     IsoSpin::IsoSpin Itotal, IsoSpin::I3 i3, Strangeness::Strange S,
		     const int imode, const int ichan, Energy & scale, 
		     const tPDVector & ,
		     const vector<Lorentz5Momentum> & momenta,
		     DecayIntegrator::MEOption) const {
  // check the total isospin
  if(Itotal!=IsoSpin::IUnknown) {
    if(Itotal==IsoSpin::IZero) {
      if(i3!=IsoSpin::I3Unknown) return vector<LorentzPolarizationVectorE>();
    }
    else if(Itotal==IsoSpin::IOne) {
      if(i3!=IsoSpin::I3Unknown&&
	 i3!=IsoSpin::I3Zero) return vector<LorentzPolarizationVectorE>();
    }
    else
      return vector<LorentzPolarizationVectorE>();
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
  if(Itotal==IsoSpin::IUnknown ||
     Itotal==IsoSpin::IZero) {
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
  ires-=5;
  // I=1 coefficient
  complex<InvEnergy> A1(ZERO);
  if(Itotal==IsoSpin::IUnknown ||
     Itotal==IsoSpin::IOne) {
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
//   for(unsigned int ix=0;ix<rhoMasses_.size();++ix) {
//     if(ix<3) output << "newdef ";
//     else     output << "insert ";
//     output << name() << ":RhoMassesI0 " << ix << " " << rhoMasses_[ix]/GeV << "\n";
//   }
//   for(unsigned int ix=0;ix<rhoWidths_.size();++ix) {
//     if(ix<3) output << "newdef ";
//     else     output << "insert ";
//     output << name() << ":RhoWidthsI0 " << ix << " " << rhoWidths_[ix]/GeV << "\n";
//   }
//   for(unsigned int ix=0;ix<omegaMasses_.size();++ix) {
//     if(ix<3) output << "newdef ";
//     else     output << "insert ";
//     output << name() << ":OmegaMassesI0 " << ix << " " << omegaMasses_[ix]/GeV << "\n";
//   }
//   for(unsigned int ix=0;ix<omegaWidths_.size();++ix) {
//     if(ix<3) output << "newdef ";
//     else     output << "insert ";
//     output << name() << ":OmegaWidthsI0 " << ix << " " << omegaWidths_[ix]/GeV << "\n";
//   }
//   output << "newdef " << name() << ":PhiMass "  << phiMass_/GeV  << "\n";
//   output << "newdef " << name() << ":PhiWidth " << phiWidth_/GeV << "\n";
//   for(unsigned int ix=0;ix<coup_I0_.size();++ix) {
//     if(ix<6) output << "newdef ";
//     else     output << "insert ";
//     output << name() << ":CouplingsI0 " << ix << " " << coup_I0_[ix]*GeV*GeV2 << "\n";
//   }
  
//   for(unsigned int ix=0;ix<rhoMasses_I1_.size();++ix) {
//     if(ix<3) output << "newdef ";
//     else     output << "insert ";
//     output << name() << ":RhoMassesI1 " << ix << " " << rhoMasses_I1_[ix]/GeV << "\n";
//   }
//   for(unsigned int ix=0;ix<rhoWidths_I1_.size();++ix) {
//     if(ix<3) output << "newdef ";
//     else     output << "insert ";
//     output << name() << ":RhoWidthsI1 " << ix << " " << rhoWidths_I1_[ix]/GeV << "\n";
//   }
//   output << "newdef " << name() << ":OmegaMass "  << omegaMass_I1_/GeV  << "\n";
//   output << "newdef " << name() << ":OmegaWidth " << omegaWidth_I1_/GeV << "\n";
//   output << "newdef " << name() << ":sigma "      << sigma_     << "\n";  
//   output << "newdef " << name() << ":GWPrefactor "      << GW_pre_*GeV     << "\n";  
//   output << "newdef " << name() << ":g_omega_pipi "      << g_omega_pi_pi_ << "\n";
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
  		    << fullName() << "\";" << endl;
}
