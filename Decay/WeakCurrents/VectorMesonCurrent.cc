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
#include "ThePEG/Helicity//VectorSpinInfo.h"
#include "Herwig++/Helicity/WaveFunction/VectorWaveFunction.h"

namespace Herwig {
using namespace ThePEG;
using Helicity::VectorWaveFunction;
using ThePEG::Helicity::VectorSpinInfo;
using ThePEG::Helicity::VectorSpinPtr;
using Helicity::Direction;
using Helicity::incoming;
using Helicity::outgoing;

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

  static ParVector<VectorMesonCurrent,Energy2> interfaceDecayConstant
    ("Decay_Constant",
     "The decay constant for the meson.",
     &VectorMesonCurrent::_decay_constant,
     0, 0, 0, 0, 100000000, false, false, true);
}

// create the decay phase space mode
bool VectorMesonCurrent::createMode(int icharge, unsigned int imode,
				    DecayPhaseSpaceModePtr mode,
				    unsigned int iloc,unsigned int ires,
				    DecayPhaseSpaceChannelPtr phase,Energy upp)
{
  // check the mode has the correct charge
  if(abs(icharge)!=abs(int(getParticleData(_id[imode])->iCharge()))){return false;}
  // check if the particle is kinematically allowed
  bool kineallowed(true);
  Energy min=getParticleData(_id[imode])->mass()
    -getParticleData(_id[imode])->widthLoCut();
  if(min>upp){kineallowed=false;}
  if(kineallowed==false){return false;}
  // construct the mode
  DecayPhaseSpaceChannelPtr newchannel(new_ptr(DecayPhaseSpaceChannel(*phase)));
  newchannel->resetDaughter(-ires,iloc);
  newchannel->init();
  mode->addChannel(newchannel);
  return kineallowed;
}

// outgoing particles 
PDVector VectorMesonCurrent::particles(int icharge, unsigned int imode, int iq, int ia)
{
  PDVector output;
  if(icharge==int(getParticleData(_id[imode])->iCharge()))
    {
      if(icharge==0)
	{
	  int iqb,iab; decayModeInfo(imode,iqb,iab);
	  if(iq==iqb&&ia==iab){output.push_back(getParticleData(_id[imode]));}
	  else{output.push_back(getParticleData(_id[imode])->CC());}
	}
      else
	{output.push_back(getParticleData(_id[imode]));}
    }
  else if(icharge==-int(getParticleData(_id[imode])->iCharge()))
    {output.push_back(getParticleData(_id[imode])->CC());}
  return output;
}

vector<LorentzPolarizationVector> 
VectorMesonCurrent::current(bool vertex, const int imode, const int ichan, 
			    Energy & scale,const ParticleVector & decay) const
{
  scale=decay[0]->mass();
  // polarization vector
  vector<LorentzPolarizationVector> temp;
  Energy fact(_decay_constant[imode]/decay[0]->mass());
  // set up the spin information for the particle
  VectorSpinPtr hwtemp;
  if(vertex)
    {
      SpinPtr vtemp(new_ptr(VectorSpinInfo(decay[0]->momentum(),true)));
      decay[0]->spinInfo(vtemp);
      hwtemp= dynamic_ptr_cast<VectorSpinPtr>(vtemp);
    }
  // calculate the vectors
  VectorWaveFunction wave(decay[0]->momentum(), decay[0]->dataPtr(), outgoing);
  for(int iy=-1;iy<2;++iy)
    {
      wave.reset(iy);
      temp.push_back(fact*wave.Wave());
      if(vertex){hwtemp->setBasisState(iy,wave.Wave());}
    }
  return temp;
}

bool VectorMesonCurrent::accept(vector<int> id)
{
  if(id.size()!=1){return false;}
  int idtemp=abs(id[0]);
  for(unsigned int ix=0;ix<_id.size();++ix)
    {if(abs(_id[ix])==idtemp){return true;}}
  return false;
}

unsigned int VectorMesonCurrent::decayMode(vector<int> idout)
{
  int idtemp(abs(idout[0])); unsigned int ix(0);
  bool found(false);
  do
    {if(idtemp==abs(_id[ix])){found=true;}++ix;}
  while(!found);
  return ix;
}

void VectorMesonCurrent::dataBaseOutput(ofstream & output)
{
  output << "create /Herwig++/VectorMesonCurrent " << fullName() << " \n";
  for(unsigned int ix=0;ix<_id.size();++ix)
    {output << "insert " << fullName() << ":ID " << ix 
	    << " " << _id[ix] << "\n";}
  for(unsigned int ix=0;ix<_decay_constant.size();++ix)
    {output << "insert " << fullName() << ":Decay_Constant " << ix 
	    << " " << _decay_constant[ix] << "\n";}
}

}
