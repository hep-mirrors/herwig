// -*- C++ -*-
//
// OneKaonTwoPionCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OneKaonTwoPionCurrent class.
//

#include "OneKaonTwoPionCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Helicity/epsilon.h"

using namespace Herwig;

DescribeClass<OneKaonTwoPionCurrent,WeakCurrent>
describeHerwigOneKaonTwoPionCurrent("Herwig::OneKaonTwoPionCurrent",
				    "HwWeakCurrents.so");
HERWIG_INTERPOLATOR_CLASSDESC(OneKaonTwoPionCurrent,Energy,Energy2)

IBPtr OneKaonTwoPionCurrent::clone() const {
  return new_ptr(*this);
}

IBPtr OneKaonTwoPionCurrent::fullclone() const {
  return new_ptr(*this);
}

OneKaonTwoPionCurrent::OneKaonTwoPionCurrent() {
  // the quarks for the different modes
  addDecayMode(2,-3);
  addDecayMode(2,-3);
  addDecayMode(2,-3);
  setInitialModes(3);
  // rho parameters
  // rho parameters for axial-vector pieces
  _rho1wgts  = {1.0,-0.145,0.};
  _rho1mass  = {0.773*GeV,1.370*GeV,1.750*GeV};
  _rho1width = {0.145*GeV,0.510*GeV,0.120*GeV};
  // K* parameters for the axial-vector pieces
  _kstar1wgts  = {1.0,-0.135,0.};
  _kstar1mass  = {0.892*GeV,1.412*GeV,1.714*GeV};
  _kstar1width = {0.050*GeV,0.227*GeV,0.323*GeV};
  // K* parameters for vector pieces
  _kstar2wgts  = {1.0,-0.25 ,-0.038};
  _kstar2mass  = {0.892*GeV,1.412*GeV,1.714*GeV};
  _kstar2width = {0.050*GeV,0.227*GeV,0.323*GeV};
  // K_1 parameters
  _k1mass  = {1.270*GeV,1.402*GeV};
  _k1width = {0.090*GeV,0.174*GeV};
  _k1wgta = {0.33,1.};
  _k1wgtb = {1.00,0.};
  // the pion decay constant
  _fpi = 130.7*MeV/sqrt(2.);
  _mpi = ZERO;
  _mK  = ZERO;
}


void OneKaonTwoPionCurrent::persistentOutput(PersistentOStream & os) const {
  os << _rho1wgts << ounit(_rho1mass,GeV) << ounit(_rho1width,GeV) 
     << _kstar1wgts << ounit(_kstar1mass,GeV) << ounit(_kstar1width,GeV) 
     << _kstar2wgts << ounit(_kstar2mass,GeV) << ounit(_kstar2width,GeV) 
     << ounit(_k1mass,GeV) << ounit(_k1width,GeV) << _k1wgta << _k1wgtb 
     << ounit(_fpi,GeV) << ounit(_mpi,GeV) << ounit(_mK,GeV);
}

void OneKaonTwoPionCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _rho1wgts >> iunit(_rho1mass,GeV) >> iunit(_rho1width,GeV)
     >> _kstar1wgts >> iunit(_kstar1mass,GeV) >> iunit(_kstar1width,GeV) 
     >> _kstar2wgts >> iunit(_kstar2mass,GeV) >> iunit(_kstar2width,GeV) 
     >> iunit(_k1mass,GeV) >> iunit(_k1width,GeV) >> _k1wgta >> _k1wgtb 
     >> iunit(_fpi,GeV) >> iunit(_mpi,GeV) >> iunit(_mK,GeV);
}


void OneKaonTwoPionCurrent::Init() {

  static ClassDocumentation<OneKaonTwoPionCurrent> documentation
    ("The OneKaonTwoPionCurrent class implements the model of "
     "Z. Phys.  C 69 (1996) 243 [arXiv:hep-ph/9503474]"
     " for the weak current with three "
     "mesons, at least one of which is a kaon",
     "The OneKaonTwoPionCurrent class implements the model of "
     "\\cite{Finkemeier:1995sr} for the weak current with three "
     "mesons, at least one of which is a kaon.",
     "\\bibitem{Finkemeier:1995sr}\n"
     "M.~Finkemeier and E.~Mirkes,\n"
     "Z.\\ Phys.\\  C {\\bf 69} (1996) 243 [arXiv:hep-ph/9503474].\n"
     " %%CITATION = ZEPYA,C69,243;%%\n");

  static Parameter<OneKaonTwoPionCurrent,Energy> interfaceFPi
    ("FPi",
     "The pion decay constant",
     &OneKaonTwoPionCurrent::_fpi, MeV, 92.4*MeV, ZERO, 200.0*MeV,
     false, false, true);

  static ParVector<OneKaonTwoPionCurrent,Energy> interfaceRhoAxialMasses
    ("RhoAxialMasses",
     "The masses for the rho resonances if used local values",
     &OneKaonTwoPionCurrent::_rho1mass, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<OneKaonTwoPionCurrent,Energy> interfaceRhoAxialWidths
    ("RhoAxialWidths",
     "The widths for the rho resonances if used local values",
     &OneKaonTwoPionCurrent::_rho1width, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<OneKaonTwoPionCurrent,Energy> interfaceKstarAxialMasses
    ("KstarAxialMasses",
     "The masses for the Kstar resonances if used local values",
     &OneKaonTwoPionCurrent::_kstar1mass, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<OneKaonTwoPionCurrent,Energy> interfaceKstarAxialWidths
    ("KstarAxialWidths",
     "The widths for the Kstar resonances if used local values",
     &OneKaonTwoPionCurrent::_kstar1width, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static ParVector<OneKaonTwoPionCurrent,Energy> interfaceKstarVectorMasses
    ("KstarVectorMasses",
     "The masses for the Kstar resonances if used local values",
     &OneKaonTwoPionCurrent::_kstar2mass, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<OneKaonTwoPionCurrent,Energy> interfaceKstarVectorWidths
    ("KstarVectorWidths",
     "The widths for the Kstar resonances if used local values",
     &OneKaonTwoPionCurrent::_kstar2width, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static ParVector<OneKaonTwoPionCurrent,double> interfaceAxialRhoWeight
    ("AxialRhoWeight",
     "The weights of the different rho resonances in the F1,2,3 form factor",
     &OneKaonTwoPionCurrent::_rho1wgts,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static ParVector<OneKaonTwoPionCurrent,double> interfaceAxialKStarWeight
    ("AxialKStarWeight",
     "The weights of the different Kstar resonances in the F1,2,3 form factor",
     &OneKaonTwoPionCurrent::_kstar1wgts,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static ParVector<OneKaonTwoPionCurrent,double> interfaceVectorKStarWeight
    ("VectorKStarWeight",
     "The weights of the different Kstar resonances in the F1,2,3 form factor",
     &OneKaonTwoPionCurrent::_kstar2wgts,
     0, 0, 0, -1000, 1000, false, false, true);

  static ParVector<OneKaonTwoPionCurrent,Energy> interfaceK1Masses
    ("K1Masses",
     "Masses of the K_1 mesons",
     &OneKaonTwoPionCurrent::_k1mass, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<OneKaonTwoPionCurrent,Energy> interfaceK1Widths
    ("K1Widths",
     "Widths of the K_1 mesons",
     &OneKaonTwoPionCurrent::_k1width, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<OneKaonTwoPionCurrent,double> interfaceK1WeightKStarPi
    ("K1WeightKStarPi",
     "The relative weights for the K_1 resonances in the K* pi final-state",
     &OneKaonTwoPionCurrent::_k1wgta, -1, 1.0, 0, 10.0,
     false, false, Interface::limited);

  static ParVector<OneKaonTwoPionCurrent,double> interfaceK1WeightRhoK
    ("K1WeightRhoK",
     "The relative weights for the K_1 resonances in the rho K final-state",
     &OneKaonTwoPionCurrent::_k1wgtb, -1, 1.0, 0, 10.0,
     false, false, Interface::limited);
}

// complete the construction of the decay mode for integration
bool OneKaonTwoPionCurrent::createMode(int icharge, tcPDPtr resonance,
				       IsoSpin::IsoSpin Itotal, IsoSpin::I3 i3,
				       unsigned int imode,PhaseSpaceModePtr mode,
				       unsigned int iloc,int ires,
				       PhaseSpaceChannel phase, Energy upp ) {
  // check the charge
  if(abs(icharge)!=3) return false; 
  // check the total isospin
  if(Itotal!=IsoSpin::IUnknown) {
    if(Itotal!=IsoSpin::IHalf) return false;
  }
  // check I_3
  if(i3!=IsoSpin::I3Unknown) {
    switch(i3) {
    case IsoSpin::I3Half:
      if(icharge ==-3) return false;
      break;
    case IsoSpin::I3MinusHalf:
      if(icharge == 3) return false;
      break;
    default:
      return false;
    }
  }
  // get the external particles and check the mass
  int iq(0),ia(0);
  tPDVector extpart(particles(1,imode,iq,ia));
  Energy min(ZERO);
  for(unsigned int ix=0;ix<extpart.size();++ix) min+=extpart[ix]->massMin();
  if(min>upp) return false;
  // the particles we will use a lot
  tPDPtr a1    = getParticleData(ParticleID::a_1minus);
  tPDPtr k1[2] = {getParticleData(ParticleID::K_1minus),
		  getParticleData(ParticleID::Kstar_1minus)};
  // the rho0 resonances
  tPDPtr rho0[3]  ={getParticleData( 113),getParticleData( 100113),
		    getParticleData( 30113)};
  // the charged rho resonances
  tPDPtr rhoc[3]  ={getParticleData(-213),getParticleData(-100213),
		    getParticleData(-30213)};
  // the K*0 resonances
  tPDPtr Kstar0[3]={getParticleData( 313),getParticleData( 100313),
		    getParticleData( 30313)};
  // the charged K* resonances
  tPDPtr Kstarc[3]={getParticleData(-323),getParticleData(-100323),
		    getParticleData(-30323)};
  if(icharge==3) {
    a1    = a1->CC();
    k1[0] = k1[0]->CC();
    k1[1] = k1[1]->CC();
    for(unsigned int ix=0;ix<3;++ix) {
      if(rhoc[ix]) rhoc[ix]=rhoc[ix]->CC();
      if(Kstar0[ix]) Kstar0[ix]=Kstar0[ix]->CC();
      if(Kstarc[ix]) Kstarc[ix]=Kstarc[ix]->CC();
    }
  }
  if(imode==0) {  
    // channels for pi0 pi0 K-
    for(unsigned int ix=0;ix<3;++ix) {
      for(unsigned int ik=0;ik<2;++ik) {
	if(resonance && resonance != k1[ik]) continue;
	mode->addChannel((PhaseSpaceChannel(phase),ires,k1[ik],ires+1,Kstarc[ix],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,k1[ik],ires+1,Kstarc[ix],ires+1,iloc+2,
			  ires+2,iloc+1,ires+2,iloc+3));
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(resonance && resonance !=Kstarc[ix]) continue;
	mode->addChannel((PhaseSpaceChannel(phase),ires,Kstarc[ix],ires+1,Kstarc[iy],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,Kstarc[ix],ires+1,Kstarc[iy],ires+1,iloc+2,
			  ires+2,iloc+1,ires+2,iloc+3));
      }
    }
  }
  else if(imode==1) {
    // channels for K- pi- pi+
    for(unsigned int ix=0;ix<3;++ix) {
      for(unsigned int ik=0;ik<2;++ik) {
	if(resonance && resonance != k1[ik]) continue;
	mode->addChannel((PhaseSpaceChannel(phase),ires,k1[ik],ires+1,rho0[ix],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,k1[ik],ires+1,Kstar0[ix],ires+1,iloc+2,
			  ires+2,iloc+1,ires+2,iloc+3));
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(resonance && resonance !=Kstarc[ix]) continue;
	mode->addChannel((PhaseSpaceChannel(phase),ires,Kstarc[ix],ires+1,rho0[ix],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,Kstarc[ix],ires+1,Kstar0[ix],ires+1,iloc+2,
			  ires+2,iloc+1,ires+2,iloc+3));
      }
    }
  }
  else if(imode==2) {
    // channels for pi- kbar0 pi0
    for(unsigned int ix=0;ix<3;++ix) {
      for(unsigned int ik=0;ik<2;++ik) {
	if(resonance && resonance != k1[ik]) continue;
	mode->addChannel((PhaseSpaceChannel(phase),ires,k1[ik],ires+1,rhoc[ix],ires+1,iloc+2,
			  ires+2,iloc+1,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,k1[ik],ires+1,Kstar0[ix],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,k1[ik],ires+1,Kstarc[ix],ires+1,iloc+3,
			  ires+2,iloc+1,ires+2,iloc+2));
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(resonance && resonance !=Kstarc[ix]) continue;
	mode->addChannel((PhaseSpaceChannel(phase),ires,Kstarc[ix],ires+1,Kstar0[iy],ires+1,iloc+1,
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,Kstarc[ix],ires+1,  rhoc[iy],ires+1,iloc+2,
			  ires+2,iloc+1,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,Kstarc[ix],ires+1,Kstarc[iy],ires+1,iloc+3,
			  ires+2,iloc+1,ires+2,iloc+2));
      }
    }
  }
  for(unsigned int ix=0;ix<_rho1mass.size();++ix) {
    mode->resetIntermediate(rhoc[ix],_rho1mass[ix],
			    _rho1width[ix]);
    mode->resetIntermediate(rho0[ix],_rho1mass[ix],
			    _rho1width[ix]);
  }
  // K star parameters in the base class
  for(unsigned int ix=0;ix<_kstar1mass.size();++ix) {
    mode->resetIntermediate(Kstarc[ix],_kstar1mass[ix],
			    _kstar1width[ix]);
    mode->resetIntermediate(Kstar0[ix],_kstar1mass[ix],
			    _kstar1width[ix]);
  }
  return true;
}

void OneKaonTwoPionCurrent::dataBaseOutput(ofstream & os,
					   bool header,bool create) const {
  if(header) os << "update decayers set parameters=\"";
  if(create) os << "create Herwig::OneKaonTwoPionCurrent " 
		<< name() << " HwWeakCurrents.so\n";
  for(unsigned int ix=0;ix<_rho1wgts.size();++ix) {
    if(ix<3) {
      os << "newdef " << name() << ":AxialRhoWeight " << ix 
	 << " " << _rho1wgts[ix] << "\n";
    }
    else {
      os << "insert " << name() << ":AxialRhoWeight " << ix 
	 << " " << _rho1wgts[ix] << "\n";
    }
  }
  for(unsigned int ix=0;ix<_kstar1wgts.size();++ix) {
    if(ix<3) {
      os << "newdef " << name() << ":AxialKStarWeight " << ix 
	 << " " << _kstar1wgts[ix] << "\n";}
    else {
      os << "insert " << name() << ":AxialKStarWeight " << ix 
	 << " " << _kstar1wgts[ix] << "\n";
    }
  }
  for(unsigned int ix=0;ix<_kstar2wgts.size();++ix) {
    if(ix<3) {
      os << "newdef " << name() << ":VectorKStarWeight " << ix 
	 << " " << _kstar2wgts[ix] << "\n";}
    else {
      os << "insert " << name() << ":VectorKStarWeight " << ix 
	 << " " << _kstar2wgts[ix] << "\n";
    }
  }
  os << "newdef " << name() << ":FPi " << _fpi/MeV << "\n";
  for(unsigned int ix=0;ix<_k1mass.size();++ix) {
    if(ix<2) {
      os << "newdef " << name() << ":K1Masses " << ix 
	 << " " << _k1mass[ix]/GeV << "\n";
    }
    else {
      os << "insert " << name() << ":K1Masses " << ix 
	 << " " << _k1mass[ix]/GeV << "\n";
    }
  }
  for(unsigned int ix=0;ix<_k1width.size();++ix) {
    if(ix<2) {
      os << "newdef " << name() << ":K1Widths " << ix 
	 << " " << _k1width[ix]/GeV << "\n";
    }
    else {
      os << "insert " << name() << ":K1Widths " << ix 
	 << " " << _k1width[ix]/GeV << "\n";
    }
  }
  for(unsigned int ix=0;ix<_rho1mass.size();++ix) {
    if(ix<3) os << "newdef " << name() << ":RhoAxialMasses " << ix 
		<< " " << _rho1mass[ix]/GeV << "\n";
    else     os << "insert " << name() << ": RhoAxialMasses" << ix 
		<< " " << _rho1mass[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_rho1width.size();++ix) {
    if(ix<3) os << "newdef " << name() << ":RhoAxialWidths " << ix 
		    << " " << _rho1width[ix]/GeV << "\n";
    else     os << "insert " << name() << ":RhoAxialWidths " << ix 
		    << " " << _rho1width[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_kstar1mass.size();++ix) {
    if(ix<3) os << "newdef " << name() << ":KstarAxialMasses " << ix 
		<< " " << _kstar1mass[ix]/GeV << "\n";
    else     os << "insert " << name() << ": KstarAxialMasses" << ix 
		<< " " << _kstar1mass[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_kstar1width.size();++ix) {
    if(ix<3) os << "newdef " << name() << ":KstarAxialWidths " << ix 
		    << " " << _kstar1width[ix]/GeV << "\n";
    else     os << "insert " << name() << ":KstarAxialWidths " << ix 
		    << " " << _kstar1width[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_kstar2mass.size();++ix) {
    if(ix<3) os << "newdef " << name() << ":KstarVectorMasses " << ix 
		<< " " << _kstar2mass[ix]/GeV << "\n";
    else     os << "insert " << name() << ": KstarVectorMasses" << ix 
		<< " " << _kstar2mass[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_kstar2width.size();++ix) {
    if(ix<3) os << "newdef " << name() << ":KstarVectorWidths " << ix 
		    << " " << _kstar2width[ix]/GeV << "\n";
    else     os << "insert " << name() << ":KstarVectorWidths " << ix 
		    << " " << _kstar2width[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_k1wgta.size();++ix) {
    if(ix<2) os << "newdef " << name() << ":K1WeightKStarPi " << ix
		<< " " << _k1wgta[ix] << "\n";
    else     os << "insert " << name() << ":K1WeightKStarPi " << ix
		<< " " << _k1wgta[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_k1wgtb.size();++ix) {
    if(ix<2) os << "newdef " << name() << ":K1WeightRhoK " << ix
		<< " " << _k1wgtb[ix] << "\n";
    else     os << "insert " << name() << ":K1WeightRhoK " << ix
		<< " " << _k1wgtb[ix] << "\n";
  }
  WeakCurrent::dataBaseOutput(os,false,false);
  if(header) os << "\n\" where BINARY ThePEGName=\"" 
		<< fullName() << "\";" << endl;
}  

void OneKaonTwoPionCurrent::doinit() {
  WeakCurrent::doinit();
  // masses for the running widths
  _mpi = getParticleData(ParticleID::piplus)->mass();
  _mK  = getParticleData(ParticleID::K0    )->mass();
}
  
Complex OneKaonTwoPionCurrent::TK1(Energy2 q2,unsigned int iopt,int ires) const {
  double denom(0);
  Complex num(0.);
  if(iopt==0) {
    if(ires>=int(_k1wgta.size())) return 0.;
    denom = std::accumulate(_k1wgta.begin(),_k1wgta.end(),0.0);
    unsigned int imin=0,imax=_k1wgta.size();
    if(ires>0) {
      imin=ires;
      imax=imin+1;
    }
    for(unsigned int ix=imin;ix<imax;++ix)
      num+=_k1wgta[ix]*Resonance::BreitWignerFW_GN(q2,_k1mass[ix],_k1width[ix]);
  }
  else if(iopt==1) {
    if(ires>=int(_k1wgtb.size())) return 0.;
    denom = std::accumulate(_k1wgtb.begin(),_k1wgtb.end(),0.0);
    unsigned int imin=0,imax=_k1wgtb.size();
    if(ires>0) {
      imin=ires;
      imax=imin+1;
    }
    for(unsigned int ix=imin;ix<imax;++ix)
      num+=_k1wgtb[ix]*Resonance::BreitWignerFW_GN(q2,_k1mass[ix],_k1width[ix]);
  }
  else assert(false);
  return num/denom;
}

// the hadronic currents    
vector<LorentzPolarizationVectorE> 
OneKaonTwoPionCurrent::current(tcPDPtr resonance,
			      IsoSpin::IsoSpin Itotal, IsoSpin::I3 i3,
			      const int imode, const int ichan, Energy & scale, 
			      const tPDVector & outgoing,
			      const vector<Lorentz5Momentum> & momenta,
			      DecayIntegrator::MEOption) const {
  // check the isospin
  if(Itotal!=IsoSpin::IUnknown && Itotal!=IsoSpin::IHalf)
    return vector<LorentzPolarizationVectorE>();
  int icharge = outgoing[0]->iCharge()+outgoing[1]->iCharge()+outgoing[2]->iCharge();
  // check I_3
  if(i3!=IsoSpin::I3Unknown) {
    switch(i3) {
    case IsoSpin::I3Half:
      if(icharge ==-3) return vector<LorentzPolarizationVectorE>();
      break;
    case IsoSpin::I3MinusHalf:
      if(icharge ==3) return vector<LorentzPolarizationVectorE>();
      break;
    default:
      return vector<LorentzPolarizationVectorE>();
    }
  }
  useMe();
  // check the resonance
  int ires1=-1;
  if(resonance) {
    switch(abs(resonance->id())/1000) {
    case 0:
      ires1=0; break;
    case 100:
      ires1=1; break;
    case  30:
      ires1=2; break;
    case  10:
      ires1=3; break;
    case  20:
      ires1=4; break;
    default:
      assert(false);
    }
  }
  // calculate q2,s1,s2,s3
  Lorentz5Momentum q;
  for(unsigned int ix=0;ix<momenta.size();++ix)
    q+=momenta[ix];
  q.rescaleMass();
  scale=q.mass();
  Energy2 q2=q.mass2();
  Energy2 s1 = (momenta[1]+momenta[2]).m2();
  Energy2 s2 = (momenta[0]+momenta[2]).m2();
  Energy2 s3 = (momenta[0]+momenta[1]).m2();
  Complex F1(0.), F2(0.), F5(0.);
  // calculate the pi0 pi0 K-
  if(imode==0) {
    if(ichan<0) {
      Complex K1fact;
      if(ires1<0)
	K1fact = TK1(q2,0,-1);
      else if(ires1<3)
	K1fact = 0.;
      else
	K1fact = TK1(q2,0,ires1-3);
      K1fact /= 6.;
      F1 = K1fact*TKstar1(s1,-1);
      F2 =-K1fact*TKstar1(s2,-1);
      if(ires1<0||ires1>2)
	F5 =-0.25*TKstar2(q2,   -1)*(TKstar1(s1,-1)-TKstar1(s2,-1));
      else
	F5 =-0.25*TKstar2(q2,ires1)*(TKstar1(s1,-1)-TKstar1(s2,-1));
    }
    else if(ichan%10==0) F1=  TK1(q2,0,0)/6.*TKstar1(s1,ichan/10);
    else if(ichan%10==1) F2= -TK1(q2,0,0)/6.*TKstar1(s2,ichan/10);
    else if(ichan%10==2) F1=  TK1(q2,0,1)/6.*TKstar1(s1,ichan/10);
    else if(ichan%10==3) F2= -TK1(q2,0,1)/6.*TKstar1(s2,ichan/10);
    else if(ichan%10<7 ) F5 =-sqrt(2.)/4*TKstar2(q2,ichan/10)*TKstar1(s1,(ichan-4)%10);
    else                 F5 = sqrt(2.)/4*TKstar2(q2,ichan/10)*TKstar1(s2,(ichan-7)%10);
  }
  // calculate the K- pi- pi+
  else if(imode==1) {
    double fact=sqrt(2.)/3.;
    if(ichan<0) {
      Complex K1facta(0.),K1factb(0.);
      if(ires1<0) {
	K1facta = TK1(q2,0,-1);
	K1factb = TK1(q2,1,-1);
      }
      else if(ires1<3) {
	K1facta = 0.;
      }
      else {
	K1facta = TK1(q2,0,ires1-3);
	K1factb = TK1(q2,1,ires1-3);
      }
      F1 = -fact*K1factb*Trho1(s1,-1);
      F2 =  fact*K1facta*TKstar1(s2,-1);
      if(ires1<0||ires1>2)
	F5 = -sqrt(0.5)*TKstar2(q2,-1)*(Trho1(s1,-1)+TKstar1(s2,-1));
      else
	F5 = -sqrt(0.5)*TKstar2(q2,ires1)*(Trho1(s1,-1)+TKstar1(s2,-1));
    }
    else if(ichan%10==0) F1 = -fact*TK1(q2,1,0)*Trho1  (s1,ichan/10);
    else if(ichan%10==1) F2 =  fact*TK1(q2,0,0)*TKstar1(s2,ichan/10);
    else if(ichan%10==2) F1 = -fact*TK1(q2,1,1)*Trho1(  s1,ichan/10);
    else if(ichan%10==3) F2 =  fact*TK1(q2,0,1)*TKstar1(s2,ichan/10);
    else if(ichan%10<7)  F5 = -sqrt(0.5)*TKstar2(q2,ichan/10)*Trho1(  s1,(ichan-4)%10);
    else                 F5 = -sqrt(0.5)*TKstar2(q2,ichan/10)*TKstar1(s2,(ichan-7)%10);
  }
  // calculate the pi- K0bar pi0
  else if(imode==2) {
    if(ichan<0) {
      Complex K1facta(0.),K1factb(0.);
      if(ires1<0) {
	K1facta = TK1(q2,0,-1);
	K1factb = TK1(q2,1,-1);
      }
      else if(ires1<3) {
	K1facta = 0.;
      }
      else {
	K1facta = TK1(q2,0,ires1-3);
	K1factb = TK1(q2,1,ires1-3);
      }
      F1 = K1facta*(TKstar1(s1,-1)-TKstar1(s3,-1))/3.;
      F2 =-(2.*K1factb*Trho1(s2,-1)+K1facta*TKstar1(s3,-1))/3.;
      if(ires1<0||ires1>2)
	F5 = -0.5*TKstar2(q2,-1)*(2.*Trho1(s2,-1)+TKstar1(s1,-1)+TKstar1(s3,-1));
      else
	F5 = -0.5*TKstar2(q2,ires1)*(2.*Trho1(s2,-1)+TKstar1(s1,-1)+TKstar1(s3,-1));
    }
    else if(ichan%15==0) F2 =-2.*TK1(q2,0,0)*Trho1  (s2,ichan/15)/3.;
    else if(ichan%15==1) F1 =    TK1(q2,1,0)*TKstar1(s1,ichan/15)/3.;
    else if(ichan%15==2) {
      F1 =-TK1(q2,1,0)*TKstar1(s3,ichan/15)/3.;
      F2 =-TK1(q2,1,0)*TKstar1(s3,ichan/15)/3.;
    }
    else if(ichan%15==3) F2 =-2.*TK1(q2,0,1)*Trho1  (s2,ichan/15)/3.;
    else if(ichan%15==4) F1 =    TK1(q2,1,1)*TKstar1(s1,ichan/15)/3.;
    else if(ichan%15==5) {
      F1 =-TK1(q2,1,1)*TKstar1(s3,ichan/15)/3.;
      F2 =-TK1(q2,1,1)*TKstar1(s3,ichan/15)/3.;
    }
    else if(ichan%15<9 ) F5 = -0.5*TKstar2(q2,ichan/15)*TKstar1(s1,(ichan- 6)%15);
    else if(ichan%15<12) F5 = -    TKstar2(q2,ichan/15)*Trho1  (s2,(ichan- 9)%15);
    else                 F5 = -0.5*TKstar2(q2,ichan/15)*TKstar1(s3,(ichan-12)%15);
  }
  // the first three form-factors
  LorentzPolarizationVectorE vect = (F2-F1)*momenta[2] + F1*momenta[1] - F2*momenta[0];
  // multiply by the transverse projection operator
  Complex dot=(vect*q)/q2;
  // scalar and parity violating terms
  vect -= dot*q;
  if(F5!=0.) 
    vect -= Complex(0.,1.)*F5/ sqr(Constants::twopi) / sqr(_fpi)*
      Helicity::epsilon(momenta[0],momenta[1],momenta[2]);
  // factor to get dimensions correct
  return vector<LorentzPolarizationVectorE>(1,q.mass()/_fpi*vect);
}

bool OneKaonTwoPionCurrent::accept(vector<int> id) {
  if(id.size()!=3) return false;
  int npip(0),npim(0),nkp(0),nkm(0);
  int npi0(0),nk0(0),nk0bar(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)       ++npip;
    else if(id[ix]==ParticleID::piminus) ++npim;
    else if(id[ix]==ParticleID::Kplus)   ++nkp;
    else if(id[ix]==ParticleID::Kminus)  ++nkm;
    else if(id[ix]==ParticleID::pi0)     ++npi0;
    else if(id[ix]==ParticleID::K0)      ++nk0;
    else if(id[ix]==ParticleID::Kbar0)   ++nk0bar;
  }
  if     ( (nkp==1&&npi0==2) || (npi0==2&&nkm==1) )   return true;
  else if( (npip==1&&npim==1&&nkp==1) ||
	   (nkm==1&&npim==1&&npip==1) )               return true;
  else if( (nk0==1&&npip==1&&npi0==1)  ||
	   (npim==1&&nk0bar==1&&npi0==1))             return true;
  else return false;
}

unsigned int OneKaonTwoPionCurrent::decayMode(vector<int> id) {
  assert(id.size()==3);
  int npip(0),npim(0),nkp(0),nkm(0);
  int npi0(0),nk0(0),nk0bar(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)       ++npip;
    else if(id[ix]==ParticleID::piminus) ++npim;
    else if(id[ix]==ParticleID::Kplus)   ++nkp;
    else if(id[ix]==ParticleID::Kminus)  ++nkm;
    else if(id[ix]==ParticleID::pi0)     ++npi0;
    else if(id[ix]==ParticleID::K0)      ++nk0;
    else if(id[ix]==ParticleID::Kbar0)   ++nk0bar;
  }
  if     ( (nkp==1&&npi0==2) || (npi0==2&&nkm==1) ) return 0;
  else if( (npip==1&&npim==1&&nkp==1) ||
	   (nkm==1&&npim==1&&npip==1) )             return 1;
  else if( (nk0==1&&npip==1&&npi0==1)  ||
	   (npim==1&&nk0bar==1&&npi0==1))           return 2;
  assert(false);
}


tPDVector OneKaonTwoPionCurrent::particles(int icharge, unsigned int imode,int,int) {
  tPDVector extpart(3);
  if(imode==0) {
    extpart[0]=getParticleData(ParticleID::pi0);
    extpart[1]=getParticleData(ParticleID::pi0);
    extpart[2]=getParticleData(ParticleID::Kminus);
  }
  else if(imode==1) {
    extpart[0]=getParticleData(ParticleID::Kminus);
    extpart[1]=getParticleData(ParticleID::piminus);
    extpart[2]=getParticleData(ParticleID::piplus);
  }
  else if(imode==2) {
    extpart[0]=getParticleData(ParticleID::piminus);
    extpart[1]=getParticleData(ParticleID::Kbar0);
    extpart[2]=getParticleData(ParticleID::pi0);
  }
  else
    assert(false);
  // conjugate the particles if needed
  if(icharge==3) {
    for(unsigned int ix=0;ix<3;++ix) {
      if(extpart[ix]->CC()) extpart[ix]=extpart[ix]->CC();
    }
  }
  // return the answer
  return extpart;
}
