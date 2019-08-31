// -*- C++ -*-
//
// KPiKStarCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the KPiKStarCurrent class.
//
//  Author: Peter Richardson
//

#include "KPiKStarCurrent.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

KPiKStarCurrent::KPiKStarCurrent() {
  // set up for the modes in the base class
  addDecayMode(2,-3);
  addDecayMode(2,-3);
  addDecayMode(2,-3);
  setInitialModes(3);
  // the weights of the different resonances in the matrix elements
  _kmag    = {1.0,0.038,0.0};
  _kphase  = {0.0,180  ,0.0};
  // model to use
  _kmodel  = 0;
  // parameter for the masses (use the parameters freom the CLEO fit 
  // rather than the PDG masses etc)
  _kstarparameters=true;
  _kstarmasses = {0.8921*GeV,1.700*GeV};
  _kstarwidths = {0.0513*GeV,0.235*GeV};
}

void KPiKStarCurrent::doinit() {
  WeakCurrent::doinit();
  // check consistency of parametrers
  if(_kstarmasses.size()!=_kstarwidths.size()) {
    throw InitException() << "Inconsistent parameters in KPiKStarCurrent"
			  << "::doinit()" << Exception::abortnow;
  }
  // the resonances
  tPDPtr res[3]={getParticleData(-323   ),
		 getParticleData(-100323),
		 getParticleData(-30323 )};
  // reset the masses in the form-factors if needed
  if(_kstarparameters&&_kstarmasses.size()<3) {
    for(unsigned int ix=0;ix<3;++ix) {
      if(res[ix+3]) _kstarmasses.push_back(res[ix]->mass());
      if(res[ix+3]) _kstarwidths.push_back(res[ix]->width());
    }
  }
  else if(!_kstarparameters) {
    _kstarmasses.clear();_kstarwidths.clear();
    for(unsigned int ix=0;ix<3;++ix) {
      if(res[ix+3]) _kstarmasses.push_back(res[ix]->mass());
      if(res[ix+3]) _kstarwidths.push_back(res[ix]->width());
    }
  }
  // set up for the Breit Wigners
  Energy mpiplus(getParticleData(ParticleID::piplus)->mass());
  Energy mk0(    getParticleData(ParticleID::K0    )->mass());
  // Kstar resonances
  for(unsigned int ix=0;ix<3;++ix) {
    _mass.push_back(_kstarmasses[ix]);
    _width.push_back(_kstarwidths[ix]);
    _massa.push_back(mk0);
    _massb.push_back(mpiplus);
    _hres.push_back(Resonance::Hhat(sqr(_mass.back()),_mass.back(),_width.back(),_massa.back(),_massb.back()));
    _dh.push_back(Resonance::dHhatds(_mass.back(),_width.back(),_massa.back(),_massb.back()));
    _h0.push_back(Resonance::H(ZERO,_mass.back(),_width.back(),_massa.back(),_massb.back(),_dh.back(),_hres.back()));
  }
  // weights for the K* channels
  if(_kmag.size()!=_kphase.size()) 
    throw InitException() << "The vectors containing the weights and phase for the"
			  << " K* channel must be the same size in "
			  << "KPiKStarCurrent::doinit()" << Exception::runerror;
  _kwgt.resize(_kmag.size());
  for(unsigned int ix=0;ix<_kmag.size();++ix) {
    double angle = _kphase[ix]/180.*Constants::pi;
    _kwgt[ix] = _kmag[ix]*(cos(angle)+Complex(0.,1.)*sin(angle));
  }
}

void KPiKStarCurrent::persistentOutput(PersistentOStream & os) const {
  os << _kmodel << _kwgt << _kmag 
     << _kphase << _kstarparameters
     << ounit(_kstarmasses,GeV) << ounit(_kstarwidths,GeV) 
     << ounit(_mass,GeV) << ounit(_width,GeV)
     << ounit(_massa,GeV) <<ounit(_massb,GeV)
     << _dh << ounit(_hres,GeV2) << ounit(_h0,GeV2);
}

void KPiKStarCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _kmodel >> _kwgt >> _kmag 
     >> _kphase >> _kstarparameters
     >> iunit(_kstarmasses,GeV) >> iunit(_kstarwidths,GeV) 
     >> iunit(_mass,GeV) >> iunit(_width,GeV)
     >> iunit(_massa,GeV) >> iunit(_massb,GeV)
     >> _dh >> iunit(_hres,GeV2) >> iunit(_h0,GeV2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<KPiKStarCurrent,WeakCurrent>
describeHerwigKPiKStarCurrent("Herwig::KPiKStarCurrent", "HwWeakCurrents.so");

void KPiKStarCurrent::Init() {
  
  static ParVector<KPiKStarCurrent,Energy> interfaceKstarMasses
    ("KstarMasses",
     "The masses of the different K* resonances for the pi pi channel",
     &KPiKStarCurrent::_kstarmasses, MeV, -1, 891.66*MeV, ZERO, 10000.*MeV,
     false, false, true);

  static ParVector<KPiKStarCurrent,Energy> interfaceKstarWidths
    ("KstarWidths",
     "The widths of the different K* resonances for the pi pi channel",
     &KPiKStarCurrent::_kstarwidths, MeV, -1, 50.8*MeV, ZERO, 1000.*MeV,
     false, false, true);

  static Switch<KPiKStarCurrent,bool> interfaceKstarParameters
    ("KstarParameters",
     "Use local values for the Kstar meson masses and widths",
     &KPiKStarCurrent::_kstarparameters, true, false, false);
  static SwitchOption interfaceKstarParameterstrue
    (interfaceKstarParameters,
     "Local",
     "Use local values",
     true);
  static SwitchOption interfaceKstarParametersParticleData
    (interfaceKstarParameters,
     "ParticleData",
     "Use the value from the particle data objects",
     false);

  static ParVector<KPiKStarCurrent,double> interfaceKMagnitude
    ("KMagnitude",
     "Magnitude of the weight of the different resonances for the K pi channel",
     &KPiKStarCurrent::_kmag, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static ParVector<KPiKStarCurrent,double> interfaceKPhase
    ("KPhase",
     "Phase of the weight of the different resonances for the K pi channel",
     &KPiKStarCurrent::_kphase, -1, 0., 0, 0,
     false, false, Interface::nolimits);
  
  static Switch<KPiKStarCurrent,int> interfaceKModel
    ("KModel",
     "The model to use for the propagator for the kaon modes.",
     &KPiKStarCurrent::_kmodel, 0, false, false);
  static SwitchOption interfaceKModelKuhn
    (interfaceKModel,
     "Kuhn",
     "The model of Kuhn and Santamaria",
     0);
  static SwitchOption interfaceKModelGounaris
    (interfaceKModel,
     "Gounaris",
     "The model of Gounaris and Sakurai.",
     1);

  static ClassDocumentation<KPiKStarCurrent> documentation
    ("The KPiKStarCurrent class is designed to implement weak"
     "decay to two scalar mesons using the models of either Kuhn and "
     "Santamaria (Z. Phys. C48, 445 (1990)) or Gounaris and Sakurai Phys. Rev. "
     "Lett. 21, 244 (1968).  The mixing parameters are taken from "
     "Phys. Rev. D61:112002,2000 (CLEO), although the PDG values for the "
     "masses and widths are used, for the decay pi+/- pi0."
     " The decay K pi is assumed to  be dominated by the lowest lying K* resonance.",
     "The weak "
     "decay current to two scalar mesons is implemented "
     "using the models of either Kuhn and "
     "Santamaria \\cite{Kuhn:1990ad} or Gounaris and Sakurai \\cite{Gounaris:1968mw}. "
     "The mixing parameters are taken from "
     "\\cite{Asner:1999kj}, although the PDG values for the "
     "masses and widths are used, for the decay pi+/- pi0."
     " The decay K pi is assumed to  be dominated by the lowest lying K* resonance.",
     "%\\cite{Kuhn:1990ad}\n"
     "\\bibitem{Kuhn:1990ad}\n"
     "  J.~H.~Kuhn and A.~Santamaria,\n"
     "  %``Tau decays to pions,''\n"
     "  Z.\\ Phys.\\  C {\\bf 48}, 445 (1990).\n"
     "  %%CITATION = ZEPYA,C48,445;%%\n"
     "%\\cite{Gounaris:1968mw}\n"
     "\\bibitem{Gounaris:1968mw}\n"
     "  G.~J.~Gounaris and J.~J.~Sakurai,\n"
     "   ``Finite width corrections to the vector meson dominance prediction for rho\n"
     "  %$\\to$ e+ e-,''\n"
     "  Phys.\\ Rev.\\ Lett.\\  {\\bf 21}, 244 (1968).\n"
     "  %%CITATION = PRLTA,21,244;%%\n"
     "%\\cite{Asner:1999kj}\n"
     "\\bibitem{Asner:1999kj}\n"
     "  D.~M.~Asner {\\it et al.}  [CLEO Collaboration],\n"
     "   ``Hadronic structure in the decay tau- --> nu/tau pi- pi0 pi0 and the  sign\n"
     "  %of the tau neutrino helicity,''\n"
     "  Phys.\\ Rev.\\  D {\\bf 61}, 012002 (2000)\n"
     "  [arXiv:hep-ex/9902022].\n"
     "  %%CITATION = PHRVA,D61,012002;%%\n"
     );

}

// complete the construction of the decay mode for integration
bool KPiKStarCurrent::createMode(int icharge, tcPDPtr resonance,
				 FlavourInfo flavour,
				 unsigned int imode,PhaseSpaceModePtr mode,
				 unsigned int iloc,int ires,
				 PhaseSpaceChannel phase, Energy upp ) {
  // check the charge
  if(abs(icharge)!=3) return false;
  // check the total isospin
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I!=IsoSpin::IHalf) return false;
  } 
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
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
  // make sure that the decays are kinematically allowed
  tPDPtr part[2];
  if(imode==0) {
    part[0]=getParticleData(ParticleID::Kplus);
    part[1]=getParticleData(ParticleID::pi0);
  }
  else if(imode==1) {
    part[0]=getParticleData(ParticleID::K0);
    part[1]=getParticleData(ParticleID::piplus);
  }
  else if(imode==2) {
    part[0]=getParticleData(ParticleID::eta);
    part[1]=getParticleData(ParticleID::Kplus);
  }
  Energy min(part[0]->massMin()+part[1]->massMin());
  if(min>upp) return false;
  // set the resonances
  // K+ pi0 or K0 pi+ or K eta decay
  tPDPtr res[3]={getParticleData(323),getParticleData(100323),getParticleData(30323)};
  if(icharge==-3) {
    for(unsigned int ix=0;ix<3;++ix) {
      if(res[ix]&&res[ix]->CC()) res[ix]=res[ix]->CC();
    }
  }
  // create the channels
  for(unsigned int ix=0;ix<3;++ix) {
    if(!res[ix]) continue;
    if(resonance && resonance != res[ix]) continue;
    mode->addChannel((PhaseSpaceChannel(phase),ires,res[ix],
		      ires+1,iloc+1,ires+1,iloc+2));
  }
  // reset the masses in the intergrators if needed
  if(_kstarparameters) {
    for(unsigned int ix=0;ix<3;++ix) {
      if(ix<_kstarmasses.size()&&res[ix]) {
	mode->resetIntermediate(res[ix],_kstarmasses[ix],_kstarwidths[ix]);
      }
    }
  }
  // return if successful
  return true;
}

// the particles produced by the current
tPDVector KPiKStarCurrent::particles(int icharge, unsigned int imode,
				     int,int) {
  tPDVector output(2);
  if(imode==0) {
    output[0]=getParticleData(ParticleID::Kplus);
    output[1]=getParticleData(ParticleID::pi0);
  }
  else if(imode==1) {
    output[0]=getParticleData(ParticleID::K0);
    output[1]=getParticleData(ParticleID::piplus);
  }
  else if(imode==2) {
    output[0]=getParticleData(ParticleID::eta);
    output[1]=getParticleData(ParticleID::Kplus);
  }
  if(icharge==-3) {
    for(unsigned int ix=0;ix<output.size();++ix) {
      if(output[ix]->CC()) output[ix]=output[ix]->CC();
    }
  }
  return output;
}

// hadronic current   
vector<LorentzPolarizationVectorE> 
KPiKStarCurrent::current(tcPDPtr resonance,
		    FlavourInfo flavour,
		    const int imode, const int ichan,Energy & scale, 
		    const tPDVector & outgoing,
		    const vector<Lorentz5Momentum> & momenta,
		    DecayIntegrator::MEOption) const {
  useMe();
  // check the isospin
  if(flavour.I!=IsoSpin::IUnknown && flavour.I!=IsoSpin::IHalf)
    return vector<LorentzPolarizationVectorE>();
  int icharge = outgoing[0]->iCharge()+outgoing[1]->iCharge();
  // check I_3
  if(flavour.I3!=IsoSpin::I3Unknown) {
    switch(flavour.I3) {
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
  // momentum difference and sum of the mesons
  Lorentz5Momentum pdiff(momenta[0]-momenta[1]);
  Lorentz5Momentum psum (momenta[0]+momenta[1]);
  psum.rescaleMass();
  scale=psum.mass();
  // mass2 of vector intermediate state
  Energy2 q2(psum.m2());
  double dot(psum*pdiff/q2);
  psum *=dot;
  LorentzPolarizationVector vect;
  // calculate the current
  unsigned int imin=0, imax=_kwgt.size();
  if(ichan>0) {
    imin = ichan;
    imax = ichan+1;
  }
  if(resonance) {
    switch(abs(resonance->id())/1000) {
    case 0:
      imin=0; break;
    case 100:
      imin=1; break;
    case  30:
      imin=2; break;
    default:
      assert(false);
    }
    imax = imin+1;
  }
  Complex denom=std::accumulate(_kwgt.begin(),_kwgt.end(),Complex(0.));
  Complex FK(0.);
  for(unsigned int ix=imin;ix<imax;++ix) {
    FK+=_kwgt[ix]*BreitWigner(q2,_kmodel,ix);
  }
  // additional prefactors
  FK/=denom;
  // single kaon/pion modes
  if     (imode==0)      FK *= sqrt(0.5);
  else if(imode==1)      FK *= 1.       ;
  // the kaon eta mode
  else if(imode==2)      FK *=sqrt(1.5);
  // compute the current
  pdiff-=psum;
  return vector<LorentzPolarizationVectorE>(1,FK*pdiff);
}
   
bool KPiKStarCurrent::accept(vector<int> id) {
  // check there are only two particles
  if(id.size()!=2) return false;
  // single charged kaon
  if((abs(id[0])==ParticleID::Kplus  &&     id[1] ==ParticleID::pi0  ) ||
     (    id[0] ==ParticleID::pi0    && abs(id[1])==ParticleID::Kplus))
    return true;
  // single neutral kaon
  else if((id[0]==ParticleID::piminus && id[1]==ParticleID::Kbar0)   ||
	  (id[0]==ParticleID::Kbar0   && id[1]==ParticleID::piminus) ||
	  (id[0]==ParticleID::piplus  && id[1]==ParticleID::K0)      ||
	  (id[0]==ParticleID::K0      && id[1]==ParticleID::piplus))
    return true;
  // charged kaon and eta
  else if((id[0]==ParticleID::Kminus && id[1]==ParticleID::eta)    ||
	  (id[0]==ParticleID::eta    && id[1]==ParticleID::Kminus) ||
	  (id[0]==ParticleID::Kplus  && id[1]==ParticleID::eta)    ||
	  (id[0]==ParticleID::eta    && id[1]==ParticleID::Kplus))
    return true;
  else
    return false;
}

// the decay mode
unsigned int KPiKStarCurrent::decayMode(vector<int> idout) {
  unsigned int imode(0),nkaon(0);
  for(unsigned int ix=0;ix<idout.size();++ix) {
    if(abs(idout[ix])==ParticleID::K0) {
      imode=1;
      ++nkaon;
    }
    else if (abs(idout[ix])==ParticleID::Kplus) {
      imode=0;
      ++nkaon;
    }
    else if (idout[ix]==ParticleID::eta) {
      imode=2;
      break;
    }
  }
  return imode;
}

// output the information for the database
void KPiKStarCurrent::dataBaseOutput(ofstream & output,bool header,
					     bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::KPiKStarCurrent " 
		    << name() << " HwWeakCurrents.so\n";
  unsigned int ix;
  for(ix=0;ix<_kstarmasses.size();++ix) {
    if(ix<2) output << "newdef ";
    else     output << "insert ";
    output << name() << ":KstarMasses " << ix << " " << _kstarmasses[ix]/MeV << "\n";
  }
  for(ix=0;ix<_kstarwidths.size();++ix) {
    if(ix<2) output << "newdef ";
    else     output << "insert ";
    output << name() << ":KstarWidths " << ix << " " << _kstarwidths[ix]/MeV << "\n";
  }
  output << "newdef " << name() << ":KstarParameters " << _kstarparameters << "\n";
  for(ix=0;ix<_kwgt.size();++ix) {
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":KMagnitude " << ix << " " << _kmag[ix]   << "\n";
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":KPhase "     << ix << " " << _kphase[ix] << "\n";
  }
  output << "newdef " << name() << ":KModel  " << _kmodel  << "\n";
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
