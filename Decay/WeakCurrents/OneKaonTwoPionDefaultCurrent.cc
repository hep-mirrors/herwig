// -*- C++ -*-
//
// OneKaonTwoPionDefaultCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OneKaonTwoPionDefaultCurrent class.
//

#include "OneKaonTwoPionDefaultCurrent.h"
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
using namespace ThePEG;

DescribeClass<OneKaonTwoPionDefaultCurrent,WeakCurrent>
describeHerwigOneKaonTwoPionDefaultCurrent("Herwig::OneKaonTwoPionDefaultCurrent",
				       "HwWeakCurrents.so");
HERWIG_INTERPOLATOR_CLASSDESC(OneKaonTwoPionDefaultCurrent,Energy,Energy2)


OneKaonTwoPionDefaultCurrent::OneKaonTwoPionDefaultCurrent() {
  // the quarks for the different modes
  addDecayMode(2,-3);
  addDecayMode(2,-3);
  addDecayMode(2,-3);
  setInitialModes(3);
  // the pion decay constant
  _fpi=130.7*MeV/sqrt(2.);
  _mpi = ZERO;
  _mK  = ZERO;
  // set the initial weights for the resonances
  // the rho weights
  _rhoF123wgts = {1.0,-0.145,0.};
  // the Kstar weights
  _kstarF123wgts = {1.};
  _kstarF5wgts   = {1.};
  // relative rho/Kstar weights
  _rhoKstarwgt=-0.2;
  // local values of the K_1 parameters
  _k1mass  = 1.402*GeV;
  _k1width = 0.174*GeV;
  // local values of the rho parameters
  _rhoF123masses = {0.773*GeV,1.370*GeV,1.750*GeV};
  _rhoF123widths = {0.145*GeV,0.510*GeV,0.120*GeV};
  // local values for the Kstar parameters
  _kstarF123masses = {0.8921*GeV};
  _kstarF123widths = {0.0513*GeV};
  _kstarF5masses   = {0.8921*GeV};
  _kstarF5widths   = {0.0513*GeV};
}

void OneKaonTwoPionDefaultCurrent::doinit() {
  WeakCurrent::doinit();
  _mpi = getParticleData(ParticleID::piplus)->mass();
  _mK  = getParticleData(ParticleID::Kminus)->mass();

}

void OneKaonTwoPionDefaultCurrent::persistentOutput(PersistentOStream & os) const {
  os << _rhoF123wgts << _kstarF123wgts << _kstarF5wgts
     << _rhoKstarwgt << ounit(_k1mass,GeV)
     << ounit(_k1width,GeV) << ounit(_fpi,GeV) << ounit(_mpi,GeV) << ounit(_mK,GeV)
     << ounit(_rhoF123masses,GeV) << ounit(_rhoF123widths,GeV) 
     << ounit(_kstarF123masses,GeV) << ounit(_kstarF5masses,GeV)
     << ounit(_kstarF123widths,GeV) << ounit(_kstarF5widths,GeV);
}

void OneKaonTwoPionDefaultCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _rhoF123wgts >> _kstarF123wgts >> _kstarF5wgts
     >> _rhoKstarwgt >> iunit(_k1mass,GeV) 
     >> iunit(_k1width,GeV) >> iunit(_fpi,GeV) >> iunit(_mpi,GeV) >> iunit(_mK,GeV)
     >> iunit(_rhoF123masses,GeV) >> iunit(_rhoF123widths,GeV) 
     >> iunit(_kstarF123masses,GeV) >> iunit(_kstarF5masses,GeV)
     >> iunit(_kstarF123widths,GeV) >> iunit(_kstarF5widths,GeV);
}

void OneKaonTwoPionDefaultCurrent::Init() {
        
  static ClassDocumentation<OneKaonTwoPionDefaultCurrent> documentation
    ("The OneKaonTwoPionDefaultCurrent class is designed to implement "
     "the three meson decays of the tau, ie pi- pi- pi+, pi0 pi0 pi-, " 
     "K- pi- K+, K0 pi- Kbar0, K- pi0 K0,pi0 pi0 K-, K- pi- pi+, "
     "pi- Kbar0 pi0, pi- pi0 eta. It uses the same currents as those in TAUOLA.",
     "The three meson decays of the tau, ie pi- pi- pi+, pi0 pi0 pi-, "
     "K- pi- K+, K0 pi- Kbar0, K- pi0 K0,pi0 pi0 K-, K- pi- pi+, "
     "and pi- Kbar0 pi0, pi- pi0 eta "
     "use the same currents as \\cite{Jadach:1993hs,Kuhn:1990ad,Decker:1992kj}.",
     "%\\cite{Jadach:1993hs}\n"
     "\\bibitem{Jadach:1993hs}\n"
     "  S.~Jadach, Z.~Was, R.~Decker and J.~H.~Kuhn,\n"
     "  %``The Tau Decay Library Tauola: Version 2.4,''\n"
     "  Comput.\\ Phys.\\ Commun.\\  {\\bf 76}, 361 (1993).\n"
     "  %%CITATION = CPHCB,76,361;%%\n"
     "%\\cite{Kuhn:1990ad}\n"
     "\\bibitem{Kuhn:1990ad}\n"
     "  J.~H.~Kuhn and A.~Santamaria,\n"
     "  %``Tau decays to pions,''\n"
     "  Z.\\ Phys.\\  C {\\bf 48}, 445 (1990).\n"
     "  %%CITATION = ZEPYA,C48,445;%%\n"
     "%\\cite{Decker:1992kj}\n"
     "\\bibitem{Decker:1992kj}\n"
     "  R.~Decker, E.~Mirkes, R.~Sauer and Z.~Was,\n"
     "  %``Tau decays into three pseudoscalar mesons,''\n"
     "  Z.\\ Phys.\\  C {\\bf 58}, 445 (1993).\n"
     "  %%CITATION = ZEPYA,C58,445;%%\n"
     );

  
  static ParVector<OneKaonTwoPionDefaultCurrent,double> interfaceF123RhoWgt
    ("F123RhoWeight",
     "The weights of the different rho resonances in the F1,2,3 form factor",
     &OneKaonTwoPionDefaultCurrent::_rhoF123wgts,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static ParVector<OneKaonTwoPionDefaultCurrent,double> interfaceF123KstarWgt
    ("F123KstarWeight",
     "The weights of the different Kstar resonances in the F1,2,3 form factor",
     &OneKaonTwoPionDefaultCurrent::_kstarF123wgts,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static ParVector<OneKaonTwoPionDefaultCurrent,double> interfaceF5KstarWgt
    ("F5KstarWeight",
     "The weights of the different Kstar resonances in the F1,2,3 form factor",
     &OneKaonTwoPionDefaultCurrent::_kstarF5wgts,
     0, 0, 0, -1000, 1000, false, false, true);
  
  static Parameter<OneKaonTwoPionDefaultCurrent,double> interfaceRhoKstarWgt
    ("RhoKstarWgt",
     "The relative weights of the rho and K* in the F5 form factor",
     &OneKaonTwoPionDefaultCurrent::_rhoKstarwgt, -0.2, -10., 10.,
     false, false, false);
  
  static Parameter<OneKaonTwoPionDefaultCurrent,Energy> interfaceK1Width
    ("K1Width",
     "The K_1 width if using local values.",
     &OneKaonTwoPionDefaultCurrent::_k1width, GeV, 0.174*GeV, ZERO, 10.0*GeV,
     false, false, false);
  
  static Parameter<OneKaonTwoPionDefaultCurrent,Energy> interfaceK1Mass
    ("K1Mass",
     "The K_1 mass if using local values.",
     &OneKaonTwoPionDefaultCurrent::_k1mass, GeV, 1.402*GeV, ZERO, 10.0*GeV,
     false, false, false);
  
  static ParVector<OneKaonTwoPionDefaultCurrent,Energy> interfacerhoF123masses
    ("rhoF123masses",
     "The masses for the rho resonances if used local values",
     &OneKaonTwoPionDefaultCurrent::_rhoF123masses, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<OneKaonTwoPionDefaultCurrent,Energy> interfacerhoF123widths
    ("rhoF123widths",
     "The widths for the rho resonances if used local values",
     &OneKaonTwoPionDefaultCurrent::_rhoF123widths, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static ParVector<OneKaonTwoPionDefaultCurrent,Energy> interfaceKstarF123masses
    ("KstarF123masses",
     "The masses for the Kstar resonances if used local values",
     &OneKaonTwoPionDefaultCurrent::_kstarF123masses, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<OneKaonTwoPionDefaultCurrent,Energy> interfaceKstarF123widths
    ("KstarF123widths",
     "The widths for the Kstar resonances if used local values",
     &OneKaonTwoPionDefaultCurrent::_kstarF123widths, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<OneKaonTwoPionDefaultCurrent,Energy> interfaceKstarF5masses
    ("KstarF5masses",
     "The masses for the Kstar resonances if used local values",
     &OneKaonTwoPionDefaultCurrent::_kstarF5masses, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);
  
  static ParVector<OneKaonTwoPionDefaultCurrent,Energy> interfaceKstarF5widths
    ("KstarF5widths",
     "The widths for the Kstar resonances if used local values",
     &OneKaonTwoPionDefaultCurrent::_kstarF5widths, GeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<OneKaonTwoPionDefaultCurrent,Energy> interfaceFPi
    ("FPi",
     "The pion decay constant",
     &OneKaonTwoPionDefaultCurrent::_fpi, MeV, 92.4*MeV, ZERO, 200.0*MeV,
     false, false, true);
}

// complete the construction of the decay mode for integration
bool OneKaonTwoPionDefaultCurrent::createMode(int icharge, tcPDPtr resonance,
					      IsoSpin::IsoSpin Itotal, IsoSpin::I3 i3, Strangeness::Strange S,
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
  // get external particles and check mass
  int iq(0),ia(0);
  tPDVector extpart(particles(1,imode,iq,ia));
  Energy min(ZERO);
  for(unsigned int ix=0;ix<extpart.size();++ix) min+=extpart[ix]->massMin();
  if(min>upp) return false;
  // the particles we will use a lot
  tPDPtr k1 = getParticleData(ParticleID::Kstar_1minus);
  if(icharge==3) k1 = k1->CC();
  // the rho0 resonances
  tPDPtr rho0[3]   = { getParticleData(113), getParticleData(100113), getParticleData(30113)};
  tPDPtr rhoc[3]   = {getParticleData(-213),getParticleData(-100213),getParticleData(-30213)};
  tPDPtr Kstar0[3] = { getParticleData(313), getParticleData(100313), getParticleData(30313)};
  tPDPtr Kstarc[3] = {getParticleData(-323),getParticleData(-100323),getParticleData(-30323)};
  if(icharge==3) {
    for(unsigned int ix=0;ix<3;++ix) {
      rhoc  [ix] =   rhoc[ix]->CC();
      Kstar0[ix] = Kstar0[ix]->CC();
      Kstarc[ix] = Kstarc[ix]->CC();
    }
  }
  if(imode==0) {
    if(resonance && resonance != k1) return false;
    // channels for pi0 pi0 K-
    for(unsigned int ix=0;ix<3;++ix) {
      mode->addChannel((PhaseSpaceChannel(phase),ires,k1,ires+1,iloc+1,ires+1,Kstarc[ix],
			ires+2,iloc+2,ires+2,iloc+3));
      mode->addChannel((PhaseSpaceChannel(phase),ires,k1,ires+1,iloc+2,ires+1,Kstarc[ix],
			ires+2,iloc+1,ires+2,iloc+3));
    }
  }
  else if(imode==1) {
    // channels for K- pi- pi+
    for(unsigned int ix=0;ix<3;++ix) {
      if(!resonance || resonance==k1) {
	mode->addChannel((PhaseSpaceChannel(phase),ires,k1,ires+1,iloc+1,ires+1,rho0[ix],
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,k1,ires+1,iloc+2,ires+1,Kstar0[ix],
			  ires+2,iloc+1,ires+2,iloc+3));
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(resonance && resonance !=Kstarc[ix]) continue;
	mode->addChannel((PhaseSpaceChannel(phase),ires,Kstarc[ix],ires+1,iloc+1,ires+1,rho0[iy],
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,Kstarc[ix],ires+1,iloc+2,ires+1,Kstar0[iy],
			  ires+2,iloc+1,ires+2,iloc+3));
      }
    }
  }
  else if(imode==2) {
    // channels for pi- kbar0 pi0
    for(unsigned int ix=0;ix<3;++ix) {
      if(!resonance || resonance==k1) {
	mode->addChannel((PhaseSpaceChannel(phase),ires,k1,ires+1,iloc+2,ires+1,rhoc[ix],
			  ires+2,iloc+1,ires+2,iloc+3));
      }
      for(unsigned int iy=0;iy<3;++iy) {
	if(resonance && resonance !=Kstarc[ix]) continue;
	mode->addChannel((PhaseSpaceChannel(phase),ires,Kstarc[ix],ires+1,iloc+1,ires+1,Kstar0[iy],
			  ires+2,iloc+2,ires+2,iloc+3));
	mode->addChannel((PhaseSpaceChannel(phase),ires,Kstarc[ix],ires+1,iloc+2,ires+1,rhoc[iy],
			  ires+2,iloc+1,ires+2,iloc+3));
      }
    }
  }
  for(unsigned int ix=0;ix<_rhoF123masses.size();++ix) {
    mode->resetIntermediate(rhoc[ix],_rhoF123masses[ix],_rhoF123widths[ix]);
    mode->resetIntermediate(rho0[ix],_rhoF123masses[ix],_rhoF123widths[ix]);
  }
  for(unsigned int ix=0;ix<_kstarF123masses.size();++ix) {
    mode->resetIntermediate(Kstarc[ix],_kstarF123masses[ix],_kstarF123widths[ix]);
    mode->resetIntermediate(Kstar0[ix],_kstarF123masses[ix],_kstarF123widths[ix]);
  }
  return true;
}

void OneKaonTwoPionDefaultCurrent::dataBaseOutput(ofstream & output,bool header,
					      bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::OneKaonTwoPionDefaultCurrent " 
		    << name() << " HwWeakCurrents.so\n";
  for(unsigned int ix=0;ix<_rhoF123wgts.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":F123RhoWeight " << ix << " " << _rhoF123wgts[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_kstarF123wgts.size();++ix) {
    if(ix<1) output << "newdef ";
    else     output << "insert ";
    output << name() << ":F123KstarWeight " << ix << " " 
	   << _kstarF123wgts[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_kstarF5wgts.size();++ix) {
    if(ix<1) output << "newdef ";
    else     output << "insert ";
    output << name() << ":F5KstarWeight " << ix << " " << _kstarF5wgts[ix] << "\n";
  }
  output << "newdef " << name() << ":RhoKstarWgt "     << _rhoKstarwgt     << "\n";
  output << "newdef " << name() << ":K1Width " << _k1width/GeV << "\n";
  output << "newdef " << name() << ":K1Mass "  << _k1mass/GeV  << "\n";
  output << "newdef " << name() << ":FPi "     << _fpi/MeV     << "\n";
  for(unsigned int ix=0;ix<_rhoF123masses.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":rhoF123masses " << ix 
	   << " " << _rhoF123masses[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_rhoF123widths.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":rhoF123widths " << ix << " " 
	   << _rhoF123widths[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_kstarF123masses.size();++ix) {
    if(ix<1) output << "newdef ";
    else     output << "insert ";
    output << name() << ":KstarF123masses " << ix << " " 
	   << _kstarF123masses[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_kstarF123widths.size();++ix) {
    if(ix<1) output << "newdef ";
    else     output << "insert ";
    output << name() << ":KstarF123widths " << ix << " " 
	   << _kstarF123widths[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_kstarF5masses.size();++ix) {
    if(ix<1) output << "newdef ";
    else     output << "insert ";
    output << name() << ":KstarF5masses " << ix << " " 
	   << _kstarF5masses[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_kstarF5widths.size();++ix) {
    if(ix<1) output << "newdef ";
    else     output << "insert ";
    output << name() << ":KstarF5widths " << ix << " " 
	   << _kstarF5widths[ix]/GeV << "\n";
  }
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

// the hadronic currents    
vector<LorentzPolarizationVectorE> 
OneKaonTwoPionDefaultCurrent::current(tcPDPtr resonance,
			      IsoSpin::IsoSpin Itotal, IsoSpin::I3 i3, Strangeness::Strange S,
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
    case  20:
      ires1=3; break;
    default:
      assert(false);
    }
  }
  useMe();
  // calculate q2,s1,s2,s3
  Lorentz5Momentum q;
  for(unsigned int ix=0;ix<momenta.size();++ix)
    q+=momenta[ix];
  q.rescaleMass();
  scale=q.mass();
  Energy2 q2=q.mass2();
  Energy2 s1 = (momenta[1]+momenta[2]).m2();
  Energy2 s2 = (momenta[0]+momenta[2]).m2();
  // calculatebthe form factors
  Complex F1(0.), F2(0.), F3(0.), F5(0.);
  // calculate the K- pi0 k0
  // calculate the pi0 pi0 K-
  Complex K1fact = ires1<0 || ires1==3 ? Resonance::BreitWignerFW_GN(q2,_k1mass,_k1width) : 0.;
  if(imode==0) {
    K1fact /=6.;
    if(ichan<0) {
      F1 = K1fact*BKstarF123(s1,-1);
      F2 =-K1fact*BKstarF123(s2,-1);
    }
    else if(ichan%2==0) F1 = K1fact*BKstarF123(s1,ichan/2);
    else                F2 =-K1fact*BKstarF123(s2,(ichan-1)/2);
  }
  // calculate the K- pi- pi+
  else if(imode==1) {
    K1fact *= sqrt(2.)/3.;
    if(ichan<0) {
      F1 =-K1fact*  BrhoF123(s1,-1);
      F2 = K1fact*BKstarF123(s2,-1);
      if(ires1<0)
	F5 = -BKstarF123(q2,   -1)*FKrho(s2,s1,-1)*sqrt(2.);
      else if(ires1<3)
	F5 = -BKstarF123(q2,ires1)*FKrho(s2,s1,-1)*sqrt(2.);
      else
	F5 = 0.;
    }
    else if(ichan%8==0) F1 =-K1fact*BrhoF123(s1,ichan/8);
    else if(ichan%8==1) F2 = K1fact*BKstarF123(s2,(ichan-1)/8);
    else                F5 = -BKstarF123(q2,ichan/8)*FKrho(s2,s1,(ichan-2)%8)*sqrt(2.);
  }
  // calculate the pi- K0bar pi0
  else if(imode==2) {
    if(ichan<0) {
      F2 =-K1fact*BrhoF123(s2,-1);
      if(ires1<0)
	F5 =-2.*BKstarF123(q2,   -1)*FKrho(s1,s2,-1);
      else if(ires1<3)
	F5 =-2.*BKstarF123(q2,ires1)*FKrho(s1,s2,-1);
      else
	F5  =0.;
    }
    else if(ichan%7==0) F2 =-K1fact*BrhoF123(s2,ichan/7);
    else                F5 =-2.*BKstarF123(q2,ichan/7)*FKrho(s1,s2,(ichan-1)%7);
  }
  // the first three form-factors
  LorentzPolarizationVectorE vect;
  vect = (F2-F1)*momenta[2]
        +(F1-F3)*momenta[1]
        +(F3-F2)*momenta[0];
  // multiply by the transverse projection operator
  Complex dot=(vect*q)/q2;
  // scalar and parity violating terms
  vect -= dot*q;
  if(F5!=0.) {
    using Constants::twopi;
    vect -= Complex(0.,1.)*F5/sqr(twopi)/sqr(_fpi)*
      Helicity::epsilon(momenta[0],momenta[1],momenta[2]);
  }
  // factor to get dimensions correct
  return vector<LorentzPolarizationVectorE>(1,q.mass()/_fpi*vect);
}

bool OneKaonTwoPionDefaultCurrent::accept(vector<int> id) {
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
  if( (nkp==1&&npi0==2) || (npi0==2&&nkm==1) )   return 0;
  else if( (npip==1&&npim==1&&nkp==1) ||
	   (nkm==1&&npim==1&&npip==1) )          return 1;
  else if( (nk0==1&&npip==1&&npi0==1)  ||
	   (npim==1&&nk0bar==1&&npi0==1))        return 2;
  return -1;
}

unsigned int OneKaonTwoPionDefaultCurrent::decayMode(vector<int> id) {
  int npip(0),npim(0),nkp(0),nkm(0);
  int  npi0(0),nk0(0),nk0bar(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(id[ix]==ParticleID::piplus)       ++npip;
    else if(id[ix]==ParticleID::piminus) ++npim;
    else if(id[ix]==ParticleID::Kplus)   ++nkp;
    else if(id[ix]==ParticleID::Kminus)  ++nkm;
    else if(id[ix]==ParticleID::pi0)     ++npi0;
    else if(id[ix]==ParticleID::K0)      ++nk0;
    else if(id[ix]==ParticleID::Kbar0)   ++nk0bar;
  }
  if( (nkp==1&&npi0==2) || (npi0==2&&nkm==1) )   return 0;
  else if( (npip==1&&npim==1&&nkp==1) ||
	   (nkm==1&&npim==1&&npip==1) )          return 1;
  else if( (nk0==1&&npip==1&&npi0==1)  ||
	   (npim==1&&nk0bar==1&&npi0==1))        return 2;
  assert(false);
}

tPDVector OneKaonTwoPionDefaultCurrent::particles(int icharge, unsigned int imode,int,int) {
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
  // conjugate the particles if needed
  if(icharge==3) {
    for(unsigned int ix=0;ix<3;++ix) {
      if(extpart[ix]->CC()) extpart[ix]=extpart[ix]->CC();
    }
  }
  // return the answer
  return extpart;
}
