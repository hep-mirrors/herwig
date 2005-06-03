// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LeptonNeutrinoCurrent class.
//

#include "LeptonNeutrinoCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "Herwig++/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "Herwig++/Helicity/WaveFunction/SpinorBarWaveFunction.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "LeptonNeutrinoCurrent.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Herwig {
using namespace ThePEG;
using Helicity::SpinorWaveFunction;
using Helicity::SpinorBarWaveFunction;
using ThePEG::Helicity::LorentzPolarizationVector;
using Helicity::Direction;
using Helicity::incoming;
using Helicity::outgoing;

LeptonNeutrinoCurrent::~LeptonNeutrinoCurrent() {}

void LeptonNeutrinoCurrent::persistentOutput(PersistentOStream & os) const {
  // no data
}

void LeptonNeutrinoCurrent::persistentInput(PersistentIStream & is, int) {
  // no data
}

ClassDescription<LeptonNeutrinoCurrent> LeptonNeutrinoCurrent::initLeptonNeutrinoCurrent;
// Definition of the static class description member.

void LeptonNeutrinoCurrent::Init() {

  static ClassDocumentation<LeptonNeutrinoCurrent> documentation
    ("The \\classname{LeptonNeutrinoCurrent} class is designed to handle the "
     "leptonic decay of the weak current.");
}


// complete the construction of the decay mode for integration
bool LeptonNeutrinoCurrent::createMode(int icharge, unsigned int imode,
				       DecayPhaseSpaceModePtr mode,
				       unsigned int iloc,unsigned int ires,
				       DecayPhaseSpaceChannelPtr phase,Energy upp)
{
  // make sure the the decays are kinematically allowed
  bool kineallowed(true);
  Energy min = getParticleData(11+2*imode)->mass()+getParticleData(12+2*imode)->mass();
  if(min>=upp){kineallowed=false; return false;}
  DecayPhaseSpaceChannelPtr newchannel;
  // set the resonances
  tPDPtr res;
  if(icharge==3){res=getParticleData(ParticleID::Wplus);}
  else if(icharge==-3){res=getParticleData(ParticleID::Wminus);}
  else{return false;}
  // create the channel
  newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
  newchannel->addIntermediate(res,0,0.0,iloc,iloc+1);
  mode->addChannel(newchannel);
  // return if successful
  return kineallowed;
}

// the particles produced by the current
PDVector LeptonNeutrinoCurrent::particles(int icharge, unsigned int imode,
					  int iq,int ia)
{
  PDVector output(2);
  if(icharge==3)
    {
      output[0]=getParticleData(-11-2*imode);
      output[1]=getParticleData(12+2*imode);
    }
  else if(icharge==-3)
    {
      output[0]=getParticleData(11+2*imode);
      output[1]=getParticleData(-12-2*imode);
    }
  return output;
}

// hadronic current   
vector<LorentzPolarizationVector> 
LeptonNeutrinoCurrent::current(bool vertex, const int imode, const int ichan,
			       Energy & scale,const ParticleVector & outpart) const
{
  Lorentz5Momentum q(outpart[0]->momentum()+outpart[1]->momentum());q.rescaleMass();
  scale=q.mass();
  // storage for the currents
  vector<LorentzPolarizationVector> temp(4);
  vector<LorentzSpinor> wave;
  vector<LorentzSpinorBar> wavebar;
  // construct the spin information objects for the  decay products and calculate
  // their wavefunctions
  if(outpart[0]->id()>0)
    {
      SpinorWaveFunction(   wave   ,outpart[1],outgoing,true,vertex);
      SpinorBarWaveFunction(wavebar,outpart[0],outgoing,true,vertex);
    }
  else
    {
      SpinorWaveFunction(   wave   ,outpart[0],outgoing,true,vertex);
      SpinorBarWaveFunction(wavebar,outpart[1],outgoing,true,vertex);
    }
  // now compute the currents
  int iloc(0);
  unsigned int ix,iy;
  for(ix=0;ix<2;++ix)
    {
      for(iy=0;iy<2;++iy)
	{
	  // location in the vector
	  if(outpart[0]->id()>0){iloc=2*iy+ix;}
	  else{iloc=2*ix+iy;}
	  // add it to the vector
	  temp[iloc]=2.*wave[ix].leftCurrent(wavebar[iy]);
	}
    }
  // return the answer
  return temp;
}

bool LeptonNeutrinoCurrent::accept(vector<int> id)
{
  bool allowed(false);
  if(id.size()!=2){return false;}
  if(abs(id[0])%2==0)
    {if((id[0]> 10&&id[0]< 18&&id[1]==-id[0]+1)||
	(id[0]<-10&&id[0]>-18&&id[1]==-id[0]-1)){allowed=true;}}
  else
    {if((id[1]> 10&&id[1]< 18&&id[0]==-id[1]+1)||
	(id[1]<-10&&id[1]>-18&&id[0]==-id[1]-1)){allowed=true;}}
  return allowed;
}

// the decay mode
unsigned int LeptonNeutrinoCurrent::decayMode(vector<int> idout)
{
  unsigned int imode=((abs(idout[0])+abs(idout[0])%2)-12)/2;
  return imode;
}

// output the information for the database
void LeptonNeutrinoCurrent::dataBaseOutput(ofstream & output)
{
  output << "create /Herwig++/LeptonNeutrinoCurrent " << fullName() << " \n";
  WeakDecayCurrent::dataBaseOutput(output);
}

}



