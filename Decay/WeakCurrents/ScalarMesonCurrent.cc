// -*- C++ -*-
//
// ScalarMesonCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ScalarMesonCurrent class.
//

#include "ScalarMesonCurrent.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void ScalarMesonCurrent::doinit() {
  unsigned int isize=numberOfModes();
  if(_id.size()!=isize||_decay_constant.size()!=isize)
    {throw InitException() << "Inconsistent parameters in ScalarMesonCurrent::doinit()"
			   << Exception::abortnow;}
  WeakCurrent::doinit();
}

ScalarMesonCurrent::ScalarMesonCurrent() {
  // the eta/eta' mixing angle
  _thetaeta=-0.194;
  // the decay constants for the different modes
  _id = {211,111,111,221,221,221,
	 331,331,331,
	 311,321,411,421,431,10431};
  _decay_constant = {130.7*MeV,130.7*MeV,130.7*MeV,130.7*MeV,130.7*MeV,130.7*MeV,
		     130.7*MeV,130.7*MeV,130.7*MeV,
		     159.8*MeV,159.8*MeV,200.0*MeV,200.0*MeV,241.0*MeV,73.7*MeV};
  addDecayMode(2,-1);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  addDecayMode(3,-3);
  addDecayMode(1,-1);
  addDecayMode(2,-2);
  addDecayMode(3,-3);
  addDecayMode(1,-3);
  addDecayMode(2,-3);
  addDecayMode(4,-1);
  addDecayMode(4,-2);
  addDecayMode(4,-3);
  addDecayMode(4,-3);
  // initial size of the arrays
  _initsize = _id.size();
  setInitialModes(_initsize);
}

void ScalarMesonCurrent::persistentOutput(PersistentOStream & os) const {
  os << _id << ounit(_decay_constant,GeV) << _thetaeta;
}

void ScalarMesonCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _id >> iunit(_decay_constant,GeV) >> _thetaeta;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ScalarMesonCurrent,WeakCurrent>
describeHerwigScalarMesonCurrent("Herwig::ScalarMesonCurrent", "HwWeakCurrents.so");

void ScalarMesonCurrent::Init() {

  static ClassDocumentation<ScalarMesonCurrent> documentation
    ("The ScalarMesonCurrent class implements the current"
     " for the decay of the weak current into a pseudoscalar meson.");

  static ParVector<ScalarMesonCurrent,long> interfaceID
    ("ID",
     "The PDG code for the outgoing meson.",
     &ScalarMesonCurrent::_id,
     0, 0, 0, -1000000, 1000000, false, false, true);

  static ParVector<ScalarMesonCurrent,Energy> interfaceDecay_Constant
    ("Decay_Constant",
     "The decay constant for the meson.",
     &ScalarMesonCurrent::_decay_constant, MeV, -1, 100.*MeV,-1000.0*MeV, 1000.0*MeV,
     false, false, true);

  static Parameter<ScalarMesonCurrent,double> interfaceThetaEtaEtaPrime
    ("ThetaEtaEtaPrime",
     "The eta-eta' mixing angle",
     &ScalarMesonCurrent::_thetaeta, -0.194, -Constants::pi, Constants::pi,
     false, false, true);

}

// create the decay phase space mode
bool ScalarMesonCurrent::createMode(int icharge, tcPDPtr resonance,
				    FlavourInfo flavour,
				    unsigned int imode,PhaseSpaceModePtr mode,
				    unsigned int iloc,int ires,
				    PhaseSpaceChannel phase, Energy upp ) {
  assert(!resonance);
  assert(flavour.I==IsoSpin::IUnknown && flavour.I3==IsoSpin::I3Unknown);
  // check the mode has the correct charge
  if(abs(icharge)!=abs(int(getParticleData(_id[imode])->iCharge()))) return false;
  // check if the particle is kinematically allowed
  tPDPtr part(getParticleData(_id[imode]));
  Energy min=part->massMin();
  if(min>upp) return false;
  // construct the mode
  mode->addChannel((PhaseSpaceChannel(phase),ires,iloc+1));
  return true;
}

// outgoing particles 
tPDVector ScalarMesonCurrent::particles(int icharge, unsigned int imode, int iq, int ia) {
  tPDPtr part(getParticleData(_id[imode]));
  tPDVector output;
  if(icharge==int(part->iCharge())) {
    if(icharge==0) {
      int iqb,iab; 
      decayModeInfo(imode,iqb,iab);
      if(iq==iqb&&ia==iab) output.push_back(part);
      else                 output.push_back(part->CC());
    }
    else {
      output.push_back(part);
    }
  }
  else if(icharge==-int(part->iCharge())) {
    output.push_back(part->CC());
  }
  return output;
}

vector<LorentzPolarizationVectorE> 
ScalarMesonCurrent::current(tcPDPtr resonance,
			    FlavourInfo flavour,
			    const int imode, const int , Energy & scale, 
			    const tPDVector & outgoing,
			    const vector<Lorentz5Momentum> & momenta,
			    DecayIntegrator::MEOption) const {
  assert(!resonance);
  assert(flavour.I==IsoSpin::IUnknown && flavour.I3==IsoSpin::I3Unknown);
  static const Complex ii(0.,1.);
  scale =momenta[0].mass();
  Complex pre(-ii*_decay_constant[imode]/scale);
  // quarks in the current
  int iq,ia;
  decayModeInfo(imode,iq,ia);
  if(abs(iq)==abs(ia)) {
    int id(outgoing[0]->id());
    if(id==ParticleID::eta) {
      if(abs(iq)==3) pre*=-2.*cos(_thetaeta)/sqrt(6.)-sin(_thetaeta)/sqrt(3.);
      else           pre*=cos(_thetaeta)/sqrt(6.)-sin(_thetaeta)/sqrt(3.);
    }
    else if(id==ParticleID::etaprime) {
      if(abs(iq)==3) pre*=-2.*sin(_thetaeta)/sqrt(6.)+cos(_thetaeta)/sqrt(3.);
      else           pre*=sin(_thetaeta)/sqrt(6.)+cos(_thetaeta)/sqrt(3.);
    }
    else if(id==ParticleID::pi0&&abs(iq)==1) {
      pre*=-sqrt(0.5);
    }
    else {
      pre*= sqrt(0.5);
    }
  }
  // return the answer
  return vector<LorentzPolarizationVectorE>(1,pre*momenta[0]);
}
  
bool ScalarMesonCurrent::accept(vector<int> id) {
  if(id.size()!=1){return false;}
  int idtemp(abs(id[0]));
  for(unsigned int ix=0;ix<_id.size();++ix) {
    if(abs(_id[ix])==idtemp) return true;
  }
  return false;
}

unsigned int ScalarMesonCurrent::decayMode(vector<int> idout) {
  int idtemp(abs(idout[0])); unsigned int ix(0);
  bool found(false);
  do {
    if(idtemp==abs(_id[ix])) found=true;
    else                     ++ix;
  }
  while(!found);
  return ix;
}

void ScalarMesonCurrent::dataBaseOutput(ofstream & output,
					bool header,bool create) const {
  if(header) {
    output << "update decayers set parameters=\"";
  }
  if(create) {
    output << "create Herwig::ScalarMesonCurrent " << name() 
	   << " HwWeakCurrents.so\n";
  }
  output << "newdef " << name() << ":ThetaEtaEtaPrime " << _thetaeta  << "\n";
  unsigned int ix;
  for(ix=0;ix<_id.size();++ix) {
    if(ix<_initsize) {
      output << "newdef " << name() << ":ID " << ix 
	     << " " << _id[ix] << "\n";
      output << "newdef " << name() << ":Decay_Constant " << ix 
	     << " " << _decay_constant[ix]/MeV << "\n";
    }
    else {
      output << "insert " << name() << ":ID " << ix 
	     << " " << _id[ix] << "\n";
      output << "insert " << name() << ":Decay_Constant " << ix 
	     << " " << _decay_constant[ix]/MeV << "\n";
    }
  }
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) {
    output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";\n";
  }
}
