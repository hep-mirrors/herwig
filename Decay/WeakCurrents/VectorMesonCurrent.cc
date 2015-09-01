// -*- C++ -*-
//
// VectorMesonCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMesonCurrent class.
//

#include "VectorMesonCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void VectorMesonCurrent::doinit() {
  unsigned int isize=numberOfModes();
  if(_id.size()!=isize||_decay_constant.size()!=isize)
    {throw InitException() << "Inconsistent parameters in VectorMesonCurrent::doinit()"
			   << Exception::abortnow;}
  WeakDecayCurrent::doinit();
}
VectorMesonCurrent::VectorMesonCurrent()  {
  _id.push_back(213);_decay_constant.push_back(0.1764*GeV2);
  addDecayMode(2,-1);
  _id.push_back(113);_decay_constant.push_back(0.1764*GeV2);
  addDecayMode(1,-1);
  _id.push_back(113);_decay_constant.push_back(0.1764*GeV2);
  addDecayMode(2,-2);
  _id.push_back(223);_decay_constant.push_back(0.1764*GeV2);
  addDecayMode(1,-1);
  _id.push_back(223);_decay_constant.push_back(0.1764*GeV2);
  addDecayMode(2,-2);
  _id.push_back(333);_decay_constant.push_back(0.2380*GeV2);
  addDecayMode(3,-3);
  _id.push_back(313);_decay_constant.push_back(0.2019*GeV2);
  addDecayMode(1,-3);
  _id.push_back(323);_decay_constant.push_back(0.2019*GeV2);
  addDecayMode(2,-3);
  _id.push_back(20213);_decay_constant.push_back(0.4626*GeV2);
  addDecayMode(2,-1);
  _id.push_back(20113);_decay_constant.push_back(0.4626*GeV2);
  addDecayMode(1,-1);
  _id.push_back(20113);_decay_constant.push_back(0.4626*GeV2);
  addDecayMode(2,-2);
  _id.push_back(413);_decay_constant.push_back(0.402*GeV2);
  addDecayMode(4,-1);
  _id.push_back(423);_decay_constant.push_back(0.402*GeV2);
  addDecayMode(4,-2);
  _id.push_back(433);_decay_constant.push_back(0.509*GeV2);
  addDecayMode(4,-3);
  _id.push_back(443);_decay_constant.push_back(1.223*GeV2);
  addDecayMode(4,-4);
  _id.push_back(100443);_decay_constant.push_back(1.08*GeV2);
  addDecayMode(4,-4);
  _id.push_back(10433);_decay_constant.push_back(0.397*GeV2);
  addDecayMode(4,-3);
  // initial size of the vectors
  _initsize=_id.size();
  setInitialModes(_initsize);
}

void VectorMesonCurrent::persistentOutput(PersistentOStream & os) const {
  os << _id << ounit(_decay_constant,GeV2);
}

void VectorMesonCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _id >> iunit(_decay_constant,GeV2);
}

ClassDescription<VectorMesonCurrent> VectorMesonCurrent::initVectorMesonCurrent;
// Definition of the static class description member.

void VectorMesonCurrent::Init() {

  static ClassDocumentation<VectorMesonCurrent> documentation
    ("The VectorMesonCurrent class implements the current"
     " for the decay of the weak current into a (pseudo)vector meson.");

  static ParVector<VectorMesonCurrent,int> interfaceID
    ("ID",
     "The PDG code for the outgoing meson.",
     &VectorMesonCurrent::_id,
     0, 0, 0, -1000000, 1000000, false, false, true);

  static ParVector<VectorMesonCurrent,Energy2> interfaceDecay_Constant
    ("Decay_Constant",
     "The decay constant for the meson.",
     &VectorMesonCurrent::_decay_constant, GeV2, -1, 1.0*GeV2,-10.0*GeV2, 10.0*GeV2,
     false, false, true);

}

// create the decay phase space mode
bool VectorMesonCurrent::createMode(int icharge, unsigned int imode,
				    DecayPhaseSpaceModePtr mode,
				    unsigned int iloc,unsigned int ires,
				    DecayPhaseSpaceChannelPtr phase,Energy upp) {
  tPDPtr part(getParticleData(_id[imode]));
  // check the mode has the correct charge
  if(abs(icharge)!=abs(int(getParticleData(_id[imode])->iCharge()))) return false;
  // check if the particle is kinematically allowed
  Energy min=part->massMin();
  if(min>upp) return false;
  // construct the mode
  DecayPhaseSpaceChannelPtr newchannel(new_ptr(DecayPhaseSpaceChannel(*phase)));
  newchannel->resetDaughter(-ires,iloc);
  mode->addChannel(newchannel);
  return true;
}

// outgoing particles 
tPDVector VectorMesonCurrent::particles(int icharge, unsigned int imode, int iq, int ia) {
  tPDPtr part(getParticleData(_id[imode]));
  tPDVector output;
  if(icharge==int(part->iCharge())) {
    if(icharge==0) {
      int iqb,iab;
      decayModeInfo(imode,iqb,iab);
      if(iq==iqb&&ia==iab) output.push_back(part);
      else                 output.push_back(part->CC());
    }
    else output.push_back(part);
  }
  else if(icharge==-int(part->iCharge())) output.push_back(part->CC());
  return output;
}

vector<LorentzPolarizationVectorE> 
VectorMesonCurrent::current(const int imode, const int, 
			    Energy & scale,const ParticleVector & decay,
			    DecayIntegrator::MEOption meopt) const {
  // set up the spin information for the particle and calculate the wavefunctions
  vector<LorentzPolarizationVector> temp;
  VectorWaveFunction::
    calculateWaveFunctions(temp,decay[0],outgoing,false);
  if(meopt==DecayIntegrator::Terminate) {
    VectorWaveFunction::constructSpinInfo(temp,decay[0],
					  outgoing,true,false);
    return vector<LorentzPolarizationVectorE>(1,LorentzPolarizationVectorE());
  }
  scale=decay[0]->mass();
  // polarization vector
  Energy fact(_decay_constant[imode]/scale);
  // quarks in the current
  int iq,ia;
  decayModeInfo(imode,iq,ia);
  if(abs(iq)==abs(ia)&&abs(iq)<3) {
    fact *= sqrt(0.5);
    if(decay[0]->id()==ParticleID::rho0&&abs(iq)==1) fact=-fact;
  }
  // normalise the current
  vector<LorentzPolarizationVectorE> returnval(3);
  for(unsigned int ix=0;ix<3;++ix) {
    returnval[ix] = temp[ix] * fact;
  }
  // return the answer
  return returnval;
}

bool VectorMesonCurrent::accept(vector<int> id) {
  if(id.size()!=1) return false;
  int idtemp(abs(id[0]));
  for(unsigned int ix=0;ix<_id.size();++ix) {
    if(abs(_id[ix])==idtemp) return true;
  }
  return false;
}

unsigned int VectorMesonCurrent::decayMode(vector<int> idout) {
  int idtemp(abs(idout[0])); unsigned int ix(0);
  bool found(false);
  do {
    if(idtemp==abs(_id[ix])) found=true;
    else                     ++ix;
  }
  while(!found);
  return ix;
}

void VectorMesonCurrent::dataBaseOutput(ofstream & output,
					bool header,bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::VectorMesonCurrent " << name() 
		    << " HwWeakCurrents.so\n";
  for(unsigned int ix=0;ix<_id.size();++ix) {
    if(ix<_initsize) {
      output << "newdef " << name() << ":ID " << ix 
	     << " " << _id[ix] << "\n";
      output << "newdef " << name() << ":Decay_Constant " << ix 
	     << " " << _decay_constant[ix]/GeV2 << "\n";
    }
    else {
      output << "insert " << name() << ":ID " << ix 
	     << " " << _id[ix] << "\n";
      output << "insert " << name() << ":Decay_Constant " << ix 
	     << " " << _decay_constant[ix]/GeV2 << "\n";
    }
  }
  WeakDecayCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
