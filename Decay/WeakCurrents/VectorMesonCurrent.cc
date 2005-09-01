// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMesonCurrent class.
//

#include "VectorMesonCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "VectorMesonCurrent.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"

namespace Herwig {
using namespace ThePEG;
using Helicity::VectorWaveFunction;
using Helicity::Direction;
using Helicity::incoming;
using Helicity::outgoing;

VectorMesonCurrent::VectorMesonCurrent() 
{
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

VectorMesonCurrent::~VectorMesonCurrent() {}

void VectorMesonCurrent::persistentOutput(PersistentOStream & os) const {
  os << _id << _decay_constant;
}

void VectorMesonCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _id >> _decay_constant;
}

ClassDescription<VectorMesonCurrent> VectorMesonCurrent::initVectorMesonCurrent;
// Definition of the static class description member.

void VectorMesonCurrent::Init() {

  static ClassDocumentation<VectorMesonCurrent> documentation
    ("The \\classname{VectorMesonCurrent} class implements the current"
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
				    DecayPhaseSpaceChannelPtr phase,Energy upp)
{
  tPDPtr part(getParticleData(_id[imode]));
  // check the mode has the correct charge
  if(abs(icharge)!=abs(int(getParticleData(_id[imode])->iCharge()))){return false;}
  // check if the particle is kinematically allowed
  bool kineallowed(true);
  Energy min=part->mass()-part->widthLoCut();
  if(min>upp){kineallowed=false;return false;}
  // construct the mode
  DecayPhaseSpaceChannelPtr newchannel(new_ptr(DecayPhaseSpaceChannel(*phase)));
  newchannel->resetDaughter(-ires,iloc);
  mode->addChannel(newchannel);
  return kineallowed;
}

// outgoing particles 
PDVector VectorMesonCurrent::particles(int icharge, unsigned int imode, int iq, int ia)
{
  tPDPtr part(getParticleData(_id[imode]));
  PDVector output;
  if(icharge==int(part->iCharge()))
    {
      if(icharge==0)
	{
	  int iqb,iab;
	  decayModeInfo(imode,iqb,iab);
	  if(iq==iqb&&ia==iab){output.push_back(part);}
	  else{output.push_back(part->CC());}
	}
      else{output.push_back(part);}
    }
  else if(icharge==-int(part->iCharge())){output.push_back(part->CC());}
  return output;
}

vector<LorentzPolarizationVector> 
VectorMesonCurrent::current(bool vertex, const int imode, const int ichan, 
			    Energy & scale,const ParticleVector & decay) const
{
  scale=decay[0]->mass();
  // polarization vector
  vector<LorentzPolarizationVector> temp;
  Energy fact(_decay_constant[imode]/scale);
  // quarks in the current
  int iq,ia;
  decayModeInfo(imode,iq,ia);
  if(abs(iq)==abs(ia)&&abs(iq)<3)
    {
      fact*=(sqrt(0.5));
      if(decay[0]->id()==ParticleID::rho0&&abs(iq)==1){fact=-fact;}
    }
  // set up the spin information for the particle and calculate the wavefunctions
  VectorWaveFunction(temp,decay[0],outgoing,true,false,vertex);
  // normalise the current
  for(unsigned int ix=0;ix<3;++ix){temp[ix]*=fact;}
  // return the answer
  return temp;
}

bool VectorMesonCurrent::accept(vector<int> id)
{
  if(id.size()!=1){return false;}
  int idtemp(abs(id[0]));
  for(unsigned int ix=0;ix<_id.size();++ix){if(abs(_id[ix])==idtemp){return true;}}
  return false;
}

unsigned int VectorMesonCurrent::decayMode(vector<int> idout)
{
  int idtemp(abs(idout[0])); unsigned int ix(0);
  bool found(false);
  do
    {
      if(idtemp==abs(_id[ix])){found=true;}
      else{++ix;}
    }
  while(!found);
  return ix;
}

void VectorMesonCurrent::dataBaseOutput(ofstream & output)
{
  output << "create /Herwig++/VectorMesonCurrent " << fullName() << " \n";
  for(unsigned int ix=0;ix<_id.size();++ix)
    {
      if(ix<_initsize)
	{
	  output << "set " << fullName() << ":ID " << ix 
		 << " " << _id[ix] << "\n";
	  output << "set " << fullName() << ":Decay_Constant " << ix 
		 << " " << _decay_constant[ix]/GeV2 << "\n";
	}
      else
	{
	  output << "insert " << fullName() << ":ID " << ix 
		 << " " << _id[ix] << "\n";
	  output << "insert " << fullName() << ":Decay_Constant " << ix 
		 << " " << _decay_constant[ix]/GeV2 << "\n";
	}
    }
}

}
