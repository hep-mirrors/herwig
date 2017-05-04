// -*- C++ -*-
//
// KPiCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the KPiCurrent class.
//

#include "KPiCurrent.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"

using namespace Herwig;
using namespace ThePEG::Helicity;
using ThePEG::Helicity::outgoing;
using ThePEG::Helicity::ScalarWaveFunction;

KPiCurrent::KPiCurrent() :
  _localparameters(true),_transverse(false), _cV(1.),_cS(0.2),
  _mpi(ZERO), _mK(ZERO) {
  // set up for the modes in the base class
  addDecayMode(2,-3);
  addDecayMode(2,-3);
  setInitialModes(2);
  // parameters for the vector resonances
  _vecmag  .push_back(1.);_vecmag  .push_back(-0.135);
  _vecphase.push_back(0.);_vecphase.push_back(180.  );
  _vecmass .push_back(891.6*MeV);_vecmass .push_back(1412.*MeV);
  _vecwidth.push_back( 50. *MeV);_vecwidth.push_back( 227.*MeV);
  // parameters for the scalar resonances
  _scamag  .push_back(0.);_scamag  .push_back(1.);
  _scaphase.push_back(0.);_scaphase.push_back(0.);
  _scamass .push_back(841.*MeV);_scamass .push_back(1429.*MeV);
  _scawidth.push_back(618.*MeV);_scawidth.push_back( 287.*MeV);
}

void KPiCurrent::persistentOutput(PersistentOStream & os) const {
  os << _cV << _cS << _localparameters 
     << ounit(_mpi,GeV) << ounit(_mK,GeV) 
     << _resmap
     << _vecmag << _vecphase << _vecwgt 
     << ounit(_vecmass,MeV) << ounit(_vecwidth,MeV)
     << _scamag << _scaphase << _scawgt 
     << ounit(_scamass,MeV) << ounit(_scawidth,MeV)
     << _transverse;
}

void KPiCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _cV >> _cS >> _localparameters 
     >> iunit(_mpi,GeV) >> iunit(_mK,GeV) 
     >> _resmap
     >> _vecmag >> _vecphase >> _vecwgt 
     >> iunit(_vecmass,MeV) >> iunit(_vecwidth,MeV) 
     >> _scamag >> _scaphase >> _scawgt 
     >> iunit(_scamass,MeV) >> iunit(_scawidth,MeV)
     >> _transverse;
}

void KPiCurrent::doinit() {
  WeakDecayCurrent::doinit();
  // check consistency of parametrers
  if(_vecmass.size()!=_vecwidth.size()||
     _scamass.size()!=_scawidth.size()) {
    throw InitException() << "Inconsistent parameters in KPiCurrent"
			  << "doinit()" << Exception::abortnow;
  }
  // the resonances
  tPDPtr vec[3]={getParticleData(-323   ),getParticleData(-100323),
		 getParticleData(-30323 )};
  tPDPtr sca[3]={getParticleData(-9000321),getParticleData(-10321)};
  // reset the masses in the form-factors if needed
  if(_localparameters) {
    if(_vecmass.size()<3) {
      for(unsigned int ix=_vecmass.size();ix<3;++ix) {
	if(vec[ix]) {
	  _vecmass.push_back( vec[ix]->mass() );
	  _vecwidth.push_back(vec[ix]->width());
	}
      }
    }
    if(_scamass.size()<2) {
      for(unsigned int ix=_scamass.size();ix<2;++ix) {
	if(sca[ix]) {
	  _scamass.push_back( sca[ix]->mass() );
	  _scawidth.push_back(sca[ix]->width());
	}
      }
    }
  }
  else {
    _vecmass.clear();_vecwidth.clear();
    for(unsigned int ix=0;ix<3;++ix) {
      if(vec[ix]) {
	_vecmass .push_back(vec[ix]->mass() );
	_vecwidth.push_back(vec[ix]->width());
      }
    }
    _scamass.clear();_scawidth.clear();
    for(unsigned int ix=0;ix<2;++ix) {
      if(sca[ix]) {
	_scamass .push_back(sca[ix]->mass() );
	_scawidth.push_back(sca[ix]->width());
      }
    }
  }
  _mpi=getParticleData(ParticleID::piplus)->mass();
  _mK =getParticleData(ParticleID::K0    )->mass();
  // weight for the vector channels
  if(_vecmag.size()!=_vecphase.size())
    throw InitException() << "The vectors containing the weights and phase for the"
			  << "vector channel must be the same size in"
			  << "KPiCurrent::doinit()" 
			  << Exception::runerror;
  _vecwgt.resize(_vecmag.size());
  for(unsigned int ix=0;ix<_vecwgt.size();++ix) {
    double angle = _vecphase[ix]/180.*Constants::pi;
    _vecwgt[ix] = _vecmag[ix]*(cos(angle)+Complex(0.,1.)*sin(angle));
  }
  // weight for the scalar channels
  if(_scamag.size()!=_scaphase.size())
    throw InitException() << "The vectors containing the weights and phase for the"
			  << "scalar channel must be the same size in"
			  << "KPiCurrent::doinit()" 
			  << Exception::runerror;
  _scawgt.resize(_scamag.size());
  for(unsigned int ix=0;ix<_scawgt.size();++ix) {
    double angle = _scaphase[ix]/180.*Constants::pi;
    _scawgt[ix] = _scamag[ix]*(cos(angle)+Complex(0.,1.)*sin(angle));
  }
  // mapping for the resonaces
  int ires(-1);
  for(unsigned int ix=0;ix<3;++ix) {
    if(vec[ix]) ++ires;
    if(ires<int(_vecwgt.size())) _resmap.push_back(ires);
  }
  if(_resmap.size()<_vecwgt.size()) {
    for(unsigned int ix=_resmap.size();ix<_vecwgt.size();++ix) {
      _resmap.push_back(-1);
    }
  }
  ires=-1;
  for(unsigned int ix=0;ix<2;++ix) {
    if(sca[ix]) ++ires;
    if(ires<int(_scawgt.size())) _resmap.push_back(ires);
  }
  if(_resmap.size()<_vecwgt.size()+_scawgt.size()) {
    for(unsigned int ix=_resmap.size()-_scawgt.size();
	ix<_scawgt.size();++ix) {
      _resmap.push_back(-1);
    }
  }
}
ClassDescription<KPiCurrent> KPiCurrent::initKPiCurrent;
// Definition of the static class description member.

void KPiCurrent::Init() {

  static ClassDocumentation<KPiCurrent> documentation
    ("The KPiCurrent class",
     "The K pi weak current has the form of \\cite{Finkemeier:1996dh}.",
     "%\\cite{Finkemeier:1996dh}\n"
     "\\bibitem{Finkemeier:1996dh}\n"
     "  M.~Finkemeier and E.~Mirkes,\n"
     "  %``The scalar contribution to tau --> K pi nu/tau,''\n"
     "  Z.\\ Phys.\\  C {\\bf 72}, 619 (1996)\n"
     "  [arXiv:hep-ph/9601275].\n"
     "  %%CITATION = ZEPYA,C72,619;%%\n"
     );

  static Parameter<KPiCurrent,double> interfacecV
    ("cV",
     "The weight for the vector contribution",
     &KPiCurrent::_cV, 1., 0., 10.0,
     false, false, Interface::limited);

  static Parameter<KPiCurrent,double> interfacecS
    ("cS",
     "The weight for the scalar contribution",
     &KPiCurrent::_cS, 0.2, -10.0, 10.0,
     false, false, Interface::limited);
  
  static ParVector<KPiCurrent,double> interfaceVectorMagnitude
    ("VectorMagnitude",
     "Magnitude of the weight for the different vector resonances",
     &KPiCurrent::_vecmag, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static ParVector<KPiCurrent,double> interfaceVectorPhase
    ("VectorPhase",
     "Phase of the weight of the different vector resonances",
     &KPiCurrent::_vecphase, -1, 0., 0, 0,
     false, false, Interface::nolimits);
  
  static ParVector<KPiCurrent,double> interfaceScalarMagnitude
    ("ScalarMagnitude",
     "Magnitude of the weight for the different scalar resonances",
     &KPiCurrent::_scamag, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static ParVector<KPiCurrent,double> interfaceScalarPhase
    ("ScalarPhase",
     "Phase of the weight of the different scalar resonances",
     &KPiCurrent::_scaphase, -1, 0., 0, 0,
     false, false, Interface::nolimits);

  static Switch<KPiCurrent,bool> interfaceLocalParameters
    ("LocalParameters",
     "Use local values for the masses and widths of the resonances or those"
     " from the ParticleData objects",
     &KPiCurrent::_localparameters, true, false, false);
  static SwitchOption interfaceLocalParametersLocal
    (interfaceLocalParameters,
     "Local",
     "Use local values",
     true);
  static SwitchOption interfaceLocalParametersParticleData
    (interfaceLocalParameters,
     "ParticleData",
     "Use the values from the particle data objects",
     false);

  static Switch<KPiCurrent,bool> interfaceTransverse
    ("Transverse",
     "Form of the vector projection operator.",
     &KPiCurrent::_transverse, false, false, false);
  static SwitchOption interfaceTransverseTransverse
    (interfaceTransverse,
     "Transverse",
     "Use 1/Q^2 in the projection operator to force it to be transverse",
     true);
  static SwitchOption interfaceTransverseMass
    (interfaceTransverse,
     "Mass",
     "Use the on-shell mass in the projection operator",
     false);

  static ParVector<KPiCurrent,Energy> interfaceVectorMass
    ("VectorMass",
     "Masses of the vector resonances",
     &KPiCurrent::_vecmass, MeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<KPiCurrent,Energy> interfaceVectorWidth
    ("VectorWidth",
     "Widths of the vector resonances",
     &KPiCurrent::_vecwidth, MeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<KPiCurrent,Energy> interfaceScalarMass
    ("ScalarMass",
     "Masses of the scalar resonances",
     &KPiCurrent::_scamass, MeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<KPiCurrent,Energy> interfaceScalarWidth
    ("ScalarWidth",
     "Widths of the scalar resonances",
     &KPiCurrent::_scawidth, MeV, -1, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, Interface::limited);
}

bool KPiCurrent::accept(vector<int> id) {
  bool allowed(false);
  // check there are only two particles
  if(id.size()!=2){return false;}
  if      ((id[0]==ParticleID::Kminus && id[1]==ParticleID::pi0)    ||
	  (id[0]==ParticleID::pi0    && id[1]==ParticleID::Kminus) ||
	  (id[0]==ParticleID::Kplus  && id[1]==ParticleID::pi0)    ||
	  (id[0]==ParticleID::pi0    && id[1]==ParticleID::Kplus)) allowed=true;
  // single neutral kaon
  else if((id[0]==ParticleID::piminus && id[1]==ParticleID::Kbar0)   ||
	  (id[0]==ParticleID::Kbar0   && id[1]==ParticleID::piminus) ||
	  (id[0]==ParticleID::piplus  && id[1]==ParticleID::K0)      ||
	  (id[0]==ParticleID::K0      && id[1]==ParticleID::piplus)) allowed=true;
  return allowed;
}

tPDVector KPiCurrent::particles(int icharge, unsigned int imode, int,int) {
  if(abs(icharge)!=3) return tPDVector();
  tPDVector output(2);
  if(imode==0) {
    output[0]=getParticleData(ParticleID::Kplus);
    output[1]=getParticleData(ParticleID::pi0);
  }
  else if(imode==1) {
    output[0]=getParticleData(ParticleID::K0);
    output[1]=getParticleData(ParticleID::piplus);
  }
  if(icharge==-3) {
    for(unsigned int ix=0;ix<output.size();++ix) {
      if(output[ix]->CC()) output[ix]=output[ix]->CC();
    }
  }
  return output;
}

unsigned int KPiCurrent::decayMode(vector<int> id) {
  unsigned int imode(0);
  for(unsigned int ix=0;ix<id.size();++ix) {
    if(abs(id[ix])==ParticleID::K0) imode=1;
  }
  return imode;
}

bool KPiCurrent::createMode(int icharge,unsigned int imode,
			    DecayPhaseSpaceModePtr mode,
			    unsigned int iloc,unsigned int,
			    DecayPhaseSpaceChannelPtr phase,Energy upp) {
  if(abs(icharge)!=3) return false; 
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
  else {
    return false;
  }
  Energy min(part[0]->massMin()+part[1]->massMin());
  if(min>upp) return false;
  DecayPhaseSpaceChannelPtr newchannel;
  // possible resonances
  tPDPtr res[5]={getParticleData(-323   ),getParticleData(-100323),
		 getParticleData(-30323 ),getParticleData(-9000321),
		 getParticleData(-10321)};
  // create the channels
  for(unsigned int ix=0;ix<5;++ix) {
    if(res[ix]) {
      newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
      newchannel->addIntermediate(res[ix],0,0.0,iloc,iloc+1);
      mode->addChannel(newchannel);
    }
  }
  // reset the masses in the intergrators if needed
  if(_localparameters) {
    // for the vectors
    for(unsigned int ix=0;ix<3;++ix) {
      if(ix<_vecmass.size()&&res[ix]) {
	mode->resetIntermediate(res[ix],_vecmass[ix],_vecwidth[ix]);
      }
    }
    // for the scalars
    for(unsigned int ix=3;ix<5;++ix) {
      if(ix-3<_scamass.size()&&res[ix]) {
	mode->resetIntermediate(res[ix],_scamass[ix-3],_scawidth[ix-3]);
      }
    }
  }
  return true;
}

void KPiCurrent::dataBaseOutput(ofstream & output,bool header,
				bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::KPiCurrent " << name() 
		    << " HeWeakCurrents.so\n";
  output << "newdef " << name() << ":LocalParameters " << _localparameters << "\n";
  output << "newdef " << name() << ":Transverse "      << _transverse << "\n";
  output << "newdef " << name() << ":cV " << _cV << "\n";
  output << "newdef " << name() << ":cS " << _cS << "\n";
  for(unsigned int ix=0;ix<_vecmag.size();++ix) {
    if(ix<2) output << "newdef ";
    else     output << "insert ";
    output << name() << ":VectorMagnitude " << ix << " " << _vecmag[ix]   << "\n";
    if(ix<3) output << "newdef ";
    else     output << "insert ";
    output << name() << ":VectorPhase "     << ix << " " << _vecphase[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_scamag.size();++ix) {
    if(ix<2) output << "newdef ";
    else     output << "insert ";
    output << name() << ":ScalarMagnitude " << ix << " " << _scamag[ix]  << "\n";
    if(ix<2) output << "newdef ";
    else     output << "insert ";
    output << name() << ":ScalarPhase "     << ix << " " << _scaphase[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_vecmass.size();++ix) {
    if(ix<2) output << "newdef ";
    else     output << "insert ";
    output << name() << ":VectorMass "  << ix << " " << _vecmass[ix]/MeV  << "\n";
    if(ix<2) output << "newdef ";
    else     output << "insert ";
    output << name() << ":VectorWidth " << ix << " " << _vecwidth[ix]/MeV << "\n";
  }
  for(unsigned int ix=0;ix<_scamass.size();++ix) {
    if(ix<2) output << "newdef ";
    else     output << "insert ";
    output << name() << ":ScalarMass "  << ix << " " << _scamass[ix]/MeV  << "\n";
    if(ix<2) output << "newdef ";
    else     output << "insert ";
    output << name() << ":ScalarWidth " << ix << " " << _scawidth[ix]/MeV << "\n";
  }
  WeakDecayCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

vector<LorentzPolarizationVectorE> 
KPiCurrent::current(const int imode, const int ichan, Energy & scale,
		    const ParticleVector & decay,
		    DecayIntegrator::MEOption meopt) const {
  useMe();
  if(meopt==DecayIntegrator::Terminate) {
    for(unsigned int ix=0;ix<2;++ix)
      ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
    return vector<LorentzPolarizationVectorE>(1,LorentzPolarizationVectorE());
  }
  // momentum difference and sum of the mesons
  Lorentz5Momentum pdiff(decay[0]->momentum()-decay[1]->momentum());
  Lorentz5Momentum psum (decay[0]->momentum()+decay[1]->momentum());
  psum.rescaleMass();
  scale=psum.mass();
  // mass2 of intermediate state
  Energy2 q2 (psum.m2());
  Energy2 dot(psum*pdiff);
  // contribution of the vector resonances
  Complex vnorm(0.),gterm(0.),sterm(0.),snorm(0.);
  complex<InvEnergy2> qterm(ZERO);
  for(unsigned int ix=0;ix<_vecwgt.size();++ix) {
    vnorm += _vecwgt[ix];
    if(ichan<0||_resmap[ix]==ichan) {
      Complex bw=_vecwgt[ix]*pWaveBreitWigner(q2,ix);
      gterm +=bw;
      qterm += _transverse ? bw/sqr(scale) : bw/sqr(_vecmass[ix]);
    }
  }
  // contribution of the scalar resonances
  for(unsigned int ix=0;ix<_scawgt.size();++ix) {
    snorm += _scawgt[ix];
    if(ichan<0||_resmap[ix+_vecwgt.size()]==ichan) {
      sterm+=_scawgt[ix]*sWaveBreitWigner(q2,ix);
    }
  }
  // compute the current
  gterm *=_cV/vnorm;
  Complex qtermnew = qterm*_cV*dot/vnorm;
  sterm *= _cS/snorm;
  LorentzPolarizationVectorE output=gterm*pdiff+(-qtermnew+sterm)*psum;
  // return the answer
  if(imode==0) output *= sqrt(0.5);
  return vector<LorentzPolarizationVectorE>(1,output);
}
