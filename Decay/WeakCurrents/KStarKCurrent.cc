// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the KStarKCurrent class.
//

#include "KStarKCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Decay/ResonanceHelpers.h"
#include "Herwig/Utilities/Kinematics.h"

using namespace Herwig;
using Kinematics::pstarTwoBodyDecay;

KStarKCurrent::KStarKCurrent() {
  using Constants::pi;
  // modes handled
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  setInitialModes(8);
  // masses for the isoscalar component
  isoScalarMasses_ = {782.65*MeV,1019.461*MeV,1425*MeV,1709*MeV,1625*MeV,2188*MeV};
  isoScalarWidths_ = {  8.49*MeV,   4.249*MeV, 215*MeV, 322*MeV, 315*MeV,  83*MeV};
  // masses for the isovector component
  isoVectorMasses_ = {775.26*MeV,1505*MeV,1720*MeV};
  isoVectorWidths_ = {149.1 *MeV, 418*MeV, 250*MeV};
  // iso scalar amplitudes
  isoScalarKStarAmp_    = {0./GeV,0./GeV,0./GeV,0./GeV,0./GeV,0./GeV};
  // isoScalarKStarAmp_    = {0./GeV,0.605/GeV,0./GeV,0.161/GeV,0./GeV,0./GeV};
  isoScalarKStarPhase_  = {    0.,      pi ,    0.,       0.,    0.,    0.};
  // iso vector amplitudes
  isoVectorKStarAmp_    = {0./GeV,0.2785/GeV,0./GeV};
  //isoVectorKStarAmp_    = {1.368/GeV,0.4464/GeV,0./GeV};
  isoVectorKStarPhase_  = {       pi,        0.,    0.};
  br4pi_                = {        0.,     0.65,    0.};
  // branching ratios
  brKK_  = 0.466;
  brPhi_ = 0.174;
}

void KStarKCurrent::doinit() {
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
  // phi mass
  mPhi_ = getParticleData(333)->mass();
  // eta mas
  mEta_ = getParticleData(221)->mass();
  // pion mass
  mpi_ = getParticleData(ParticleID::piplus)->mass();
}

IBPtr KStarKCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr KStarKCurrent::fullclone() const {
  return new_ptr(*this);
}

void KStarKCurrent::persistentOutput(PersistentOStream & os) const {
  os << ounit(isoScalarMasses_,GeV) << ounit(isoScalarWidths_,GeV)
     << ounit(isoVectorMasses_,GeV) << ounit(isoVectorWidths_,GeV)
     << ounit(isoScalarKStarAmp_,1./GeV) << ounit(isoVectorKStarAmp_,1./GeV)
     << isoScalarKStarPhase_ << isoVectorKStarPhase_
     << ounit(isoScalarKStarCoup_,1./GeV) << ounit(isoVectorKStarCoup_,1./GeV)
     << brKK_ << brPhi_ << br4pi_
     << ounit(mpi_,GeV) << ounit(mPhi_,GeV) << ounit(mEta_,GeV);
}

void KStarKCurrent::persistentInput(PersistentIStream & is, int) {
  is >> iunit(isoScalarMasses_,GeV) >> iunit(isoScalarWidths_,GeV)
     >> iunit(isoVectorMasses_,GeV) >> iunit(isoVectorWidths_,GeV)
     >> iunit(isoScalarKStarAmp_,1./GeV) >> iunit(isoVectorKStarAmp_,1./GeV)
     >> isoScalarKStarPhase_ >> isoVectorKStarPhase_
     >> iunit(isoScalarKStarCoup_,1./GeV) >> iunit(isoVectorKStarCoup_,1./GeV)
     >> brKK_ >> brPhi_ >> br4pi_
     >> iunit(mpi_,GeV) >> iunit(mPhi_,GeV) >> iunit(mEta_,GeV);
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<KStarKCurrent,WeakCurrent>
describeHerwigKStarKCurrent("Herwig::KStarKCurrent", "HwWeakCurrents.so");

void KStarKCurrent::Init() {

  static ClassDocumentation<KStarKCurrent> documentation
    ("There is no documentation for the KStarKCurrent class");

}


// complete the construction of the decay mode for integration
bool KStarKCurrent::createMode(int icharge, tcPDPtr resonance,
			     FlavourInfo flavour,
			     unsigned int imode,PhaseSpaceModePtr mode,
			     unsigned int iloc,int ires,
			     PhaseSpaceChannel phase, Energy upp ) {
  if(icharge!=0) return false;
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I==IsoSpin::IZero) {
      if(flavour.I3!=IsoSpin::I3Unknown || flavour.I3!=IsoSpin::I3Zero) return false;
    }
    else if(flavour.I==IsoSpin::IOne) {
      if(flavour.I3!=IsoSpin::I3Unknown&&flavour.I3!=IsoSpin::I3Zero) return false;
     }
    else
      return false;
  }
  // check the kinematics
  int iq(0),ia(0);
  tPDVector out=particles(icharge,imode,iq,ia);
  if(out[0]->mass()+out[1]->massMin()>upp) return false;
  // resonances we need
  tPDPtr omega[6] = {getParticleData(    223),getParticleData(    333),
		     getParticleData( 100223),getParticleData( 100333),
		     getParticleData(  30223),getParticleData( 100333)};//getParticleData(  30333)};
  tPDPtr rho0[3]  = {getParticleData( 113),getParticleData( 100113),getParticleData( 30113)};
  // I=0 channels
  if(flavour.I==IsoSpin::IUnknown || flavour.I==IsoSpin::IZero) {
    for(unsigned int ix=0;ix<6;++ix) {
      if(resonance && resonance != omega[ix]) continue;
      mode->addChannel((PhaseSpaceChannel(phase),ires,omega[ix],ires+1,iloc+1,ires+1,iloc+2));
    }
  }
  // I=1 channels
  if(flavour.I==IsoSpin::IUnknown || flavour.I==IsoSpin::IOne) {
    for(unsigned int ix=0;ix<3;++ix) {
      if(resonance && resonance != rho0[ix]) continue;
      mode->addChannel((PhaseSpaceChannel(phase),ires,rho0[ix],ires+1,iloc+1,ires+1,iloc+2));
    }
  }
  return true;
}

// the particles produced by the current
tPDVector KStarKCurrent::particles(int icharge, unsigned int imode,
				   int,int) {
  assert(icharge==0);
  if(imode==0 || imode==1) {
    return {getParticleData(ParticleID::Kminus),getParticleData(ParticleID::Kstarplus)};
  }
  else if(imode==2 || imode==3) {
    return {getParticleData(ParticleID::Kbar0),getParticleData(ParticleID::Kstar0)};
  }
  else if(imode==4 || imode==5) {
    return {getParticleData(ParticleID::Kplus),getParticleData(ParticleID::Kstarminus)};
  }
  else if(imode==6 || imode==7) {
    return {getParticleData(ParticleID::K0),getParticleData(ParticleID::Kstarbar0)};
  }
  else
    assert(false);
}


void KStarKCurrent::constructSpinInfo(ParticleVector decay) const {
  vector<LorentzPolarizationVector> temp(3);
  for(unsigned int ix=0;ix<3;++ix) {
    temp[ix] = HelicityFunctions::polarizationVector(-decay[1]->momentum()
						     ,ix,Helicity::outgoing);
  }
  ScalarWaveFunction::constructSpinInfo(decay[0],outgoing,true);
  VectorWaveFunction::constructSpinInfo(temp,decay[1],
					outgoing,true,true);
}

// hadronic current   
vector<LorentzPolarizationVectorE> 
KStarKCurrent::current(tcPDPtr resonance,
		     FlavourInfo flavour,
		     const int imode, const int ichan, Energy & scale, 
		     const tPDVector & outgoing,
		     const vector<Lorentz5Momentum> & momenta,
		     DecayIntegrator::MEOption) const {
  // check the total isospin
  if(flavour.I==IsoSpin::IHalf)
    return vector<LorentzPolarizationVectorE>();
  // check I3
  if(flavour.I3!=IsoSpin::I3Unknown&&flavour.I3!=IsoSpin::I3Zero)
    return vector<LorentzPolarizationVectorE>();
  // using this current
  useMe();
  // polarization vectors for the K*
  vector<LorentzPolarizationVector> temp(3);
  for(unsigned int ix=0;ix<3;++ix)
    temp[ix] = HelicityFunctions::polarizationVector(-momenta[1],ix,Helicity::outgoing);
  // calculate q2
  Lorentz5Momentum q = momenta[0]+momenta[1];
  q.rescaleMass();
  scale=q.mass();
  Energy2 q2=q.mass2();
  complex<InvEnergy> pre(ZERO);
  if((flavour.I==IsoSpin::IUnknown || flavour.I==IsoSpin::IOne) && ichan<6) {
    unsigned int imin=0, imax = 6;
    if(resonance) {
      if(ichan>0) {
	imin = ichan;
	imax = ichan+1;
      }
      switch(resonance->id()/1000) {
      case 0:
	imin = 0;
	break;
      case 100:
	imin = 1;
	break;
      case 30 :
	imin = 2;
	break;
      default:
	assert(false);
      }
      if(resonance->id()%1000==223)
	imin=2*imin;
      else
	imin=2*imin+1;
      imax=imin+1;
    }
    for(unsigned int ix=imin;ix<imax;++ix) {
      if(ix!=3) {
	pre += isoScalarKStarCoup_[ix]*Resonance::BreitWignerFW(q2,isoScalarMasses_[ix],isoScalarWidths_[ix]);
      }
      else {
	Energy m1 = outgoing[0]->mass(), m2 = outgoing[1]->mass();
	double r1(0.);
	if(m1+m2<q.mass() && m1+m2<isoScalarMasses_[ix])
	  r1 = pstarTwoBodyDecay(            q.mass(),m1,m2)/
	       pstarTwoBodyDecay(isoScalarMasses_[ix],m1,m2);
	double r2=0.;
	if(mEta_+mPhi_<q.mass() && mEta_+mPhi_<isoScalarMasses_[ix])
	  r2 = pstarTwoBodyDecay(            q.mass(),mEta_,mPhi_)/
	       pstarTwoBodyDecay(isoScalarMasses_[ix],mEta_,mPhi_);
	Energy gam = isoScalarWidths_[ix]*
	  (brKK_*pow(r1,3)+brPhi_*pow(r2,3)+1.-brKK_-brPhi_);
	Energy2 mR2 = sqr(isoScalarMasses_[ix]);
	pre += isoScalarKStarCoup_[ix]*mR2/(mR2-q2-Complex(0.,1.)*gam*q.mass());
      }
    }
  }
  double isoSign(1.);
  if(imode==2||imode==3||imode==6||imode==7)
    isoSign=-1.;
  if((flavour.I==IsoSpin::IUnknown || flavour.I==IsoSpin::IZero) &&
     (ichan<0 || ichan>=6)) {
    unsigned int imin=0, imax = 3;
    if(ichan>0) {
      imin = ichan-6;
      imax = imin+1;
    }
    if(resonance) {
      switch(resonance->id()/1000) {
      case 0:
	imin = 0;
	break;
      case 100:
	imin = 1;
	break;
      case 30 :
	imin = 2;
	break;
      default:
	assert(false);
      }
      imax=imin+1;
    }
    for(unsigned int ix=imin;ix<imax;++ix) {
      Energy2 mR2 = sqr(isoVectorMasses_[ix]);
      Energy wid = isoVectorWidths_[ix]*
	(1.-br4pi_[ix]+ br4pi_[ix]*mR2/q2*pow((q2-16.*sqr(mpi_))/(mR2-16.*sqr(mpi_)),1.5));
      pre += isoSign*isoVectorKStarCoup_[ix]*mR2/(mR2-q2-Complex(0.,1.)*q.mass()*wid);


      cerr << "testing in current " << ix << " " << isoVectorKStarCoup_[ix]*GeV << "\n";
      cerr << "testing width " << q.mass()/GeV << " " << wid/GeV << "\n";
    }
  }
  // calculate the current
  vector<LorentzPolarizationVectorE> ret(3);
  for(unsigned int ix=0;ix<3;++ix) {
    ret[ix] = pre*Helicity::epsilon(q,temp[ix],momenta[1]);
  }
  cerr << momenta[0]/GeV << " " << momenta[0].mass()/GeV << "\n";
  cerr << momenta[1]/GeV << " " << momenta[1].mass()/GeV << "\n";
  return ret;
}
   
bool KStarKCurrent::accept(vector<int> id) {
  if(id.size()!=2) return false;
  if(( id[0] == ParticleID::Kminus && id[1] == ParticleID::Kstarplus ) ||
     ( id[1] == ParticleID::Kminus && id[0] == ParticleID::Kstarplus ))
    return true;
  else if(( id[0] == ParticleID::Kbar0 && id[1] == ParticleID::Kstar0 ) ||
	  ( id[1] == ParticleID::Kbar0 && id[0] == ParticleID::Kstar0 ))
    return true;
  else if(( id[0] == ParticleID::Kplus && id[1] == ParticleID::Kstarminus ) ||
	  ( id[1] == ParticleID::Kplus && id[0] == ParticleID::Kstarminus ))
    return true;
  else if(( id[0] == ParticleID::K0 && id[1] == ParticleID::Kstarbar0 ) ||
	  ( id[1] == ParticleID::K0 && id[0] == ParticleID::Kstarbar0 ))
    return true;
  else
    return false;
}

// the decay mode
unsigned int KStarKCurrent::decayMode(vector<int> id) {
  assert(id.size()==2);
  if(( id[0] == ParticleID::Kminus && id[1] == ParticleID::Kstarplus ) ||
     ( id[1] == ParticleID::Kminus && id[0] == ParticleID::Kstarplus ))
    return 0;
  else if(( id[0] == ParticleID::Kbar0 && id[1] == ParticleID::Kstar0 ) ||
	  ( id[1] == ParticleID::Kbar0 && id[0] == ParticleID::Kstar0 ))
    return 2;
  else if(( id[0] == ParticleID::Kplus && id[1] == ParticleID::Kstarminus ) ||
	  ( id[1] == ParticleID::Kplus && id[0] == ParticleID::Kstarminus ))
    return 4;
  else if(( id[0] == ParticleID::K0 && id[1] == ParticleID::Kstarbar0 ) ||
	  ( id[1] == ParticleID::K0 && id[0] == ParticleID::Kstarbar0 ))
    return 6;
  else
    assert(false);
}

// output the information for the database
void KStarKCurrent::dataBaseOutput(ofstream & output,bool header,
					bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::KStarKCurrent " 
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
