// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ScalarMesonCurrent class.
//

#include "ScalarMesonCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ScalarMesonCurrent.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig++/Helicity/WaveFunction/ScalarWaveFunction.h"

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;
using Helicity::ScalarWaveFunction;
using Helicity::outgoing;

ScalarMesonCurrent::ScalarMesonCurrent() 
{
  // the eta/eta' mixing angle
  _thetaeta=-0.194;
  // the decay constants for the different modes
  _id.push_back(211);_decay_constant.push_back(130.7*MeV);
  addDecayMode(2,-1);
  _id.push_back(111);_decay_constant.push_back(130.7*MeV);
  addDecayMode(1,-1);
  _id.push_back(111);_decay_constant.push_back(130.7*MeV);
  addDecayMode(2,-2);
  _id.push_back(221);_decay_constant.push_back(130.7*MeV);
  addDecayMode(1,-1);
  _id.push_back(221);_decay_constant.push_back(130.7*MeV);
  addDecayMode(2,-2);
  _id.push_back(221);_decay_constant.push_back(130.7*MeV);
  addDecayMode(3,-3);
  _id.push_back(331);_decay_constant.push_back(130.7*MeV);
  addDecayMode(1,-1);
  _id.push_back(331);_decay_constant.push_back(130.7*MeV);
  addDecayMode(2,-2);
  _id.push_back(331);_decay_constant.push_back(130.7*MeV);
  addDecayMode(3,-3);
  _id.push_back(311);_decay_constant.push_back(159.8*MeV);
  addDecayMode(1,-3);
  _id.push_back(321);_decay_constant.push_back(159.8*MeV);
  addDecayMode(2,-3);
  _id.push_back(411);_decay_constant.push_back(200.0*MeV);
  addDecayMode(4,-1);
  _id.push_back(421);_decay_constant.push_back(200.0*MeV);
  addDecayMode(4,-2);
  _id.push_back(431);_decay_constant.push_back(241.0*MeV);
  addDecayMode(4,-3);
  _id.push_back(10431);_decay_constant.push_back(73.7*MeV);
  addDecayMode(4,-3);
  // initial size of the arrays
  _initsize = _id.size();
  setInitialModes(_initsize);
}

void ScalarMesonCurrent::persistentOutput(PersistentOStream & os) const {
  os << _id << _decay_constant << _thetaeta;
}

void ScalarMesonCurrent::persistentInput(PersistentIStream & is, int) {
  is >> _id >> _decay_constant >> _thetaeta;
}

ClassDescription<ScalarMesonCurrent> ScalarMesonCurrent::initScalarMesonCurrent;
// Definition of the static class description member.

void ScalarMesonCurrent::Init() {

  static ClassDocumentation<ScalarMesonCurrent> documentation
    ("The ScalarMesonCurrent class implements the current"
     " for the decay of the weak current into a pseudoscalar meson.");

  static ParVector<ScalarMesonCurrent,int> interfaceID
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
     &ScalarMesonCurrent::_thetaeta, -0.194, -pi, pi,
     false, false, true);

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
  tPDPtr part(getParticleData(_id[imode]));
  Energy min=part->mass()-part->widthLoCut();
  if(min>upp){kineallowed=false; return false;}
  // construct the mode
  DecayPhaseSpaceChannelPtr newchannel(new_ptr(DecayPhaseSpaceChannel(*phase)));
  newchannel->resetDaughter(-ires,iloc);
  mode->addChannel(newchannel);
  return true;
}

// outgoing particles 
PDVector ScalarMesonCurrent::particles(int icharge, unsigned int imode, int iq, int ia)
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
      else
	{output.push_back(part);}
    }
  else if(icharge==-int(part->iCharge()))
    {output.push_back(part->CC());}
  return output;
}

vector<LorentzPolarizationVector> 
ScalarMesonCurrent::current(bool vertex, const int imode, const int, 
			    Energy & scale,const ParticleVector & decay) const
{
  static const Complex ii(0.,1.);
  scale = decay[0]->mass();
  // workaround for gcc 3.2.3 bug*
  // set up the spin information for the particle
  //ALB ScalarWaveFunction(decay[0],outgoing,true,vertex);
  PPtr mytemp = decay[0];
  ScalarWaveFunction(mytemp,outgoing,true,vertex);

  Complex pre(-ii*_decay_constant[imode]/scale);
  // quarks in the current
  int iq,ia;
  decayModeInfo(imode,iq,ia);
  if(abs(iq)==abs(ia))
    {
      int id(decay[0]->id());
      if(id==ParticleID::eta)
	{
	  if(abs(iq)==3){pre*=-2.*cos(_thetaeta)/sqrt(6.)-sin(_thetaeta)/sqrt(3.);}
	  else{pre*=cos(_thetaeta)/sqrt(6.)-sin(_thetaeta)/sqrt(3.);}
	}
      else if(id==ParticleID::etaprime)
	{
	  if(abs(iq)==3){pre*=-2.*sin(_thetaeta)/sqrt(6.)+cos(_thetaeta)/sqrt(3.);}
	  else{pre*=sin(_thetaeta)/sqrt(6.)+cos(_thetaeta)/sqrt(3.);}
	}
      else if(id==ParticleID::pi0&&abs(iq)==1){pre*=-sqrt(0.5);}
      else{pre*= sqrt(0.5);}
    }
  // return the answer
  return vector<LorentzPolarizationVector>(1,pre*decay[0]->momentum());
}

bool ScalarMesonCurrent::accept(vector<int> id)
{
  if(id.size()!=1){return false;}
  int idtemp(abs(id[0]));
  for(unsigned int ix=0;ix<_id.size();++ix)
    {if(abs(_id[ix])==idtemp){return true;}}
  return false;
}

unsigned int ScalarMesonCurrent::decayMode(vector<int> idout)
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

void ScalarMesonCurrent::dataBaseOutput(ofstream & output,bool header,bool create) const
{
  if(header){output << "update decayers set parameters=\"";}
  if(create)
    {output << "create Herwig++::ScalarMesonCurrent " << fullName() << " \n";}
  output << "set " << fullName() << ":ThetaEtaEtaPrime " << _thetaeta  << "\n";
  unsigned int ix;
  for(ix=0;ix<_id.size();++ix)
    {
      if(ix<_initsize)
	{
	  output << "set " << fullName() << ":ID " << ix 
		 << " " << _id[ix] << "\n";
	  output << "set " << fullName() << ":Decay_Constant " << ix 
	    << " " << _decay_constant[ix]/MeV << "\n";
	}
      else
	{
	  output << "insert " << fullName() << ":ID " << ix 
		 << " " << _id[ix] << "\n";
	  output << "insert " << fullName() << ":Decay_Constant " << ix 
	    << " " << _decay_constant[ix]/MeV << "\n";
	}
    }
  WeakDecayCurrent::dataBaseOutput(output,false,false);
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}

}
