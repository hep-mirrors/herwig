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

using namespace Herwig;

KKPiCurrent::KKPiCurrent() {
  // masses for the isoscalar component
  isoScalarMasses_ = {782.65*MeV,1019.461*MeV,1425*MeV,1680*MeV,1625*MeV,2188*MeV};
  isoScalarWidths_ = {  8.49*MeV,   4.249*MeV, 215*MeV, 150*MeV, 315*MeV,  83*MeV};
  // masses for the isovector component
  isoVectorMasses_ = {775.26*MeV,1465*MeV,1720*MeV};
  isoVectorWidths_ = {149.1 *MeV, 400*MeV, 250*MeV};
  // amplitude and phases for the isoscalar
  isoScalarKStarAmp_  ={ZERO,ZERO,ZERO,11./GeV,ZERO,ZERO};
  isoScalarKStarPhase_={  0.,  0.,  0.,     0.,  0.,  0.};
  // amplitudes and phase for the isovector component
  isoVectorKStarAmp_  ={ZERO,ZERO,4.5/GeV};
  isoVectorKStarPhase_={0.,0.,Constants::pi};
  // Coupling for the K* to Kpi
  gKStar_ = 5.38;
  // mstar masses
  mKStarP_ = 895.5 *MeV;
  mKStar0_ = 895.55*MeV;
  wKStarP_ = 46.2  *MeV;
  wKStar0_ = 47.3  *MeV;
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
  for(unsigned int ix=0;ix<isoScalarKStarAmp_.size();++ix)
    isoScalarKStarCoup_.push_back(isoScalarKStarAmp_[ix]*(cos(isoScalarKStarPhase_[ix])
						+ii*sin(isoScalarKStarPhase_[ix])));
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
			     IsoSpin::IsoSpin Itotal, IsoSpin::I3 i3,
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
  tPDPtr resI0[3] = {getParticleData(   113),getParticleData(100113),getParticleData( 30113)};
  tPDPtr resI1[6] = {getParticleData(   223),getParticleData(   333),
		     getParticleData(100223),getParticleData(100333),
		     getParticleData( 30223),getParticleData( 30333)};
  tPDPtr res[2];
  if(imode==0) {
    res[0] = getParticleData(ParticleID::Kstar0);
    res[0] = getParticleData(ParticleID::Kstarbar0);
  }
  else if(imode==1) {
    res[0] = getParticleData(ParticleID::Kstarplus);
    res[0] = getParticleData(ParticleID::Kstarminus);
  }
  else if(imode==2||imode==4) {
    res[0] = getParticleData(ParticleID::Kstarplus);
    res[0] = getParticleData(ParticleID::Kstarbar0);
  }
  else if(imode==3||imode==5) {
    res[0] = getParticleData(ParticleID::Kstarminus);
    res[0] = getParticleData(ParticleID::Kstar0);
  }
  for(unsigned int ix=0;ix<3;++ix) {
    mode->addChannel((PhaseSpaceChannel(phase),ires,resI0[ix],ires+1,res[0],ires+1,iloc+2,
		      ires+2,iloc+1,ires+2,iloc+3));
    mode->addChannel((PhaseSpaceChannel(phase),ires,resI0[ix],ires+1,res[1],ires+1,iloc+1,
		      ires+2,iloc+2,ires+2,iloc+3));
  }
  for(unsigned int ix=0;ix<6;++ix) {
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
		     IsoSpin::IsoSpin Itotal, IsoSpin::I3 i3,
		     const int, const int ichan, Energy & scale, 
		     const tPDVector & ,
		     const vector<Lorentz5Momentum> & momenta,
		     DecayIntegrator::MEOption) const {
//   // check the total isospin
//   if(Itotal!=IsoSpin::IUnknown) {
//     if(Itotal==IsoSpin::IZero) {
//       if(i3!=IsoSpin::I3Unknown) return vector<LorentzPolarizationVectorE>();
//     }
//     else if(Itotal==IsoSpin::IOne) {
//       if(i3!=IsoSpin::I3Unknown&&
// 	 i3!=IsoSpin::I3One) return vector<LorentzPolarizationVectorE>();
//     }
//     else
//       return vector<LorentzPolarizationVectorE>();
//   }
//   useMe();
//   // calculate q2,s1,s2,s3
//   Lorentz5Momentum q;
//   for(unsigned int ix=0;ix<momenta.size();++ix) q+=momenta[ix];
//   q.rescaleMass();
//   scale=q.mass();
//   Energy2 q2=q.mass2();
//   Energy2 sm = (momenta[1]+momenta[2]).m2();
//   Energy2 sp = (momenta[0]+momenta[2]).m2();
//   Energy2 s0 = (momenta[0]+momenta[1]).m2();
//   int irho=-1;
//   if(ichan>=0) {
//     irho = ichan%3;
//   }
//   // isospin zero part of the current
//   complex<InvEnergy3> F_I0(ZERO);
//   if(Itotal==IsoSpin::IUnknown || Itotal==IsoSpin::IOne) {
//     // compute H rho
//     Complex Hrho = HChannel(irho,rhoMasses_[0],rhoWidths_[0],sp,sm,s0,mpip_,mpi0_);
//     // terms in the current
//     if((!resonance || resonance->id() == 223) && ichan<=2) {
//       F_I0 += Hrho*coup_I0_[0]*Resonance::BreitWignerFW(q2,omegaMasses_[0],omegaWidths_[0]);
//     }
//     if((!resonance || resonance->id() == 333) && (ichan<0 || (ichan>=9&&ichan<=11))) {
//       F_I0 += Hrho*coup_I0_[1]*Resonance::BreitWignerFW(q2,phiMass_       ,phiWidth_      );
//     }
//     if((!resonance || resonance->id() == 100223) && (ichan<0 || (ichan>=3&&ichan<=5))) {
//       F_I0 += Hrho*coup_I0_[2]*Resonance::BreitWignerFW(q2,omegaMasses_[1],omegaWidths_[1]);
//     }
//     if((!resonance || resonance->id() == 30223) && (ichan<0 || (ichan>=6&&ichan<=8))) {
//       F_I0 +=  Hrho*coup_I0_[3]*Resonance::BreitWignerFW(q2,omegaMasses_[2],omegaWidths_[2]);
//     }
//     if((!resonance || resonance->id() == 333) && (ichan<0 || (ichan>=12&&ichan<=14))) {
//       F_I0 += coup_I0_[4]*HChannel(irho,rhoMasses_[1],rhoWidths_[1],sp,sm,s0,mpip_,mpi0_)*
// 	Resonance::BreitWignerFW(q2,phiMass_,phiWidth_);
//     }
//     if((!resonance || resonance->id() == 100223) && (ichan<0 || (ichan>=15&&ichan<=17))) {
//       F_I0 += coup_I0_[5]*HChannel(irho,rhoMasses_[2],rhoWidths_[2],sp,sm,s0,mpip_,mpi0_)*
//     Resonance::BreitWignerFW(q2,omegaMasses_[2],omegaWidths_[2]);
//     }
//   }
//   // isospin = 1
//   complex<InvEnergy3> F_I1(ZERO);
//   if((Itotal==IsoSpin::IUnknown || Itotal==IsoSpin::IZero) && !resonance && ichan<0) {
//     F_I1 = GW_*
//       Resonance::BreitWignerFW(q2,omegaMass_I1_,omegaWidth_I1_)/sqr(omegaMass_I1_)*
//       (Resonance::BreitWignerPWave(s0,rhoMasses_I1_[0],
// 				   rhoWidths_I1_[0],mpip_,mpip_)/sqr(rhoMasses_I1_[0])+
//        sigma_*Resonance::BreitWignerPWave(s0,rhoMasses_I1_[1],
// 					  rhoWidths_I1_[1],mpip_,mpip_)/sqr(rhoMasses_I1_[0]));
//   }
//   // the current
//   LorentzPolarizationVector vect = (F_I0+F_I1)*
//     Helicity::epsilon(momenta[0],
//   		      momenta[1],
//   		      momenta[2]);
//   // factor to get dimensions correct
//   return vector<LorentzPolarizationVectorE>(1,q.mass()*vect);
  assert(false);
}
   
bool KKPiCurrent::accept(vector<int> id) {
  if(id.size()!=3){return false;}
  unsigned int npiplus(0),npi0(0),npiminus(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus) ++npiplus;
    else if(id[ix]==ParticleID::piminus) ++npiplus;
    else if(id[ix]==ParticleID::pi0)    ++npi0;
  }
  return (npiplus==1&&npiminus==1&&npi0==1);
}

// the decay mode
unsigned int KKPiCurrent::decayMode(vector<int> ) {
  assert(false);
  return 0;
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
