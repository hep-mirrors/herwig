// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ScalarMesonCurrent class.
//

#include "ScalarMesonCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Interface/ParVector.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ScalarMesonCurrent.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity//ScalarSpinInfo.h"

namespace Herwig {
using namespace ThePEG;
using ThePEG::Helicity::ScalarSpinInfo;
using ThePEG::Helicity::tcSpinfoPtr;

ScalarMesonCurrent::~ScalarMesonCurrent() {}

void ScalarMesonCurrent::persistentOutput(PersistentOStream & os) const {
  os << _id << _decay_constant;
}

void ScalarMesonCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _id >> _decay_constant;
}

ClassDescription<ScalarMesonCurrent> ScalarMesonCurrent::initScalarMesonCurrent;
// Definition of the static class description member.

void ScalarMesonCurrent::Init() {

  static ClassDocumentation<ScalarMesonCurrent> documentation
    ("The \\classname{ScalarMesonCurrent} class implements the current"
     " for the decay of the weak current into a pseudoscalar meson.");

  static ParVector<ScalarMesonCurrent,int> interfaceID
    ("ID",
     "The PDG code for the outgoing meson.",
     &ScalarMesonCurrent::_id,
     0, 0, 0, -1000000, 1000000, false, false, true);

  static ParVector<ScalarMesonCurrent,Energy> interfaceDecayConstant
    ("Decay_Constant",
     "The decay constant for the meson.",
     &ScalarMesonCurrent::_decay_constant,
     0, 0, 0, 0, 10000000, false, false, true);

}

// create the decay phase space mode
bool ScalarMesonCurrent::createMode(int icharge,unsigned int imode,
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
  return true;
}

// outgoing particles 
  PDVector ScalarMesonCurrent::particles(int icharge, unsigned int imode, int iq, int ia)
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
ScalarMesonCurrent::current(bool vertex, const int imode, const int ichan, 
			    const Particle & inpart,
			    const ParticleVector & decay) const
{
  unsigned int imes=decay.size()-1;
  // polarization vector
  vector<LorentzPolarizationVector> temp;
  temp.push_back(_decay_constant[imode]*(decay[imes]->momentum())/inpart.mass());
  // set up the spin information for the particle
  if(vertex)
    {
      SpinPtr stemp = new_ptr(ScalarSpinInfo(decay[imes]->momentum(),true));
      decay[imes]->spinInfo(stemp);
    }
  return temp;
}

bool ScalarMesonCurrent::accept(vector<int> id)
{
  if(id.size()!=1){return false;}
  int idtemp=abs(id[0]);
  for(unsigned int ix=0;ix<_id.size();++ix)
    {if(abs(_id[ix])==idtemp){return true;}}
  return false;
}

unsigned int ScalarMesonCurrent::decayMode(vector<int> idout)
{
  int idtemp(abs(idout[0])); unsigned int ix(0);
  bool found(false);
  do
    {if(idtemp==abs(_id[ix])){found=true;}++ix;}
  while(!found);
  return ix;
}

void ScalarMesonCurrent::dataBaseOutput(ofstream & output)
{
  output << "create /Herwig++/ScalarMesonCurrent " << fullName() << " \n";
  for(unsigned int ix=0;ix<_id.size();++ix)
    {output << "insert " << fullName() << ":ID " << ix 
	    << " " << _id[ix] << "\n";}
  for(unsigned int ix=0;ix<_decay_constant.size();++ix)
    {output << "insert " << fullName() << ":Decay_Constant " << ix 
	    << " " << _decay_constant[ix] << "\n";}
}

}

/*

 */
