// -*- C++ -*-
//
// LeptonNeutrinoCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LeptonNeutrinoCurrent class.
//

#include "LeptonNeutrinoCurrent.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
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

NoPIOClassDescription<LeptonNeutrinoCurrent> 
LeptonNeutrinoCurrent::initLeptonNeutrinoCurrent;
// Definition of the static class description member.

void LeptonNeutrinoCurrent::Init() {

  static ClassDocumentation<LeptonNeutrinoCurrent> documentation
    ("The LeptonNeutrinoCurrent class is designed to handle the "
     "leptonic decay of the weak current.");
}


// complete the construction of the decay mode for integration
bool LeptonNeutrinoCurrent::createMode(int icharge, unsigned int imode_in,
				       DecayPhaseSpaceModePtr mode,
				       unsigned int iloc,unsigned int,
				       DecayPhaseSpaceChannelPtr phase,Energy upp) {
  int imode = imode_in;
  // make sure the the decays are kinematically allowed
  Energy min = getParticleData(11+2*imode)->mass()+getParticleData(12+2*imode)->mass();
  if(min>=upp) return false;
  DecayPhaseSpaceChannelPtr newchannel;
  // set the resonances
  tPDPtr res;
  if(icharge==3)       res=getParticleData(ParticleID::Wplus);
  else if(icharge==-3) res=getParticleData(ParticleID::Wminus);
  else                 return false;
  // create the channel
  newchannel=new_ptr(DecayPhaseSpaceChannel(*phase));
  newchannel->addIntermediate(res,0,0.0,iloc,iloc+1);
  mode->addChannel(newchannel);
  // return if successful
  return true;
}

// the particles produced by the current
tPDVector LeptonNeutrinoCurrent::particles(int icharge, unsigned int imode_in,
					  int,int)
{
  int imode = imode_in;
  tPDVector output(2);
  if(icharge==3)
    {
      int id = -11-2*imode;
      output[0]=getParticleData(id);
      output[1]=getParticleData(12+2*imode);
    }
  else if(icharge==-3)
    {
      output[0]=getParticleData(11+2*imode);
      int id = -12-2*imode;
      output[1]=getParticleData(id);
    }
  return output;
}

// hadronic current   
vector<LorentzPolarizationVectorE> 
LeptonNeutrinoCurrent::current(const int, const int,
			       Energy & scale,const ParticleVector & outpart,
			       DecayIntegrator::MEOption meopt) const {
  Lorentz5Momentum q(outpart[0]->momentum()+outpart[1]->momentum());q.rescaleMass();
  scale=q.mass();
  // storage for the currents
  vector<LorentzPolarizationVectorE> temp(4);
  vector<LorentzSpinor<SqrtEnergy> > wave;
  vector<LorentzSpinorBar<SqrtEnergy> > wavebar; 
  // their wavefunctions
  if(outpart[0]->id()>0) {
    SpinorWaveFunction   ::
      calculateWaveFunctions(wave   ,outpart[1],outgoing);
    SpinorBarWaveFunction::
      calculateWaveFunctions(wavebar,outpart[0],outgoing);
  }
  else {
    SpinorWaveFunction   ::
      calculateWaveFunctions(wave   ,outpart[0],outgoing);
    SpinorBarWaveFunction::
      calculateWaveFunctions(wavebar,outpart[1],outgoing);
  }
  if(meopt==DecayIntegrator::Terminate) {
    if(outpart[0]->id()>0) {
      SpinorWaveFunction   ::constructSpinInfo(wave   ,outpart[1],outgoing,true);
      SpinorBarWaveFunction::constructSpinInfo(wavebar,outpart[0],outgoing,true);
    }
    else {
      SpinorWaveFunction   ::constructSpinInfo(   wave,outpart[0],outgoing,true);
      SpinorBarWaveFunction::constructSpinInfo(wavebar,outpart[1],outgoing,true);
    }
  }
  // now compute the currents
  int iloc(0);
  unsigned int ix,iy;
  for(ix=0;ix<2;++ix) {
    for(iy=0;iy<2;++iy) {
      iloc = outpart[0]->id()>0 ? 2*iy+ix : 2*ix+iy;
      temp[iloc]=2.*wave[ix].leftCurrent(wavebar[iy]);
    }
  }
  // return the answer
  return temp;
}

bool LeptonNeutrinoCurrent::accept(vector<int> id) {
  bool allowed(false);
  if(id.size()!=2) return false;
  if(abs(id[0])%2==0) {
    if((id[0]> 10&&id[0]< 18&&id[1]==-id[0]+1)||
       (id[0]<-10&&id[0]>-18&&id[1]==-id[0]-1)) allowed=true;
  }
  else {
    if((id[1]> 10&&id[1]< 18&&id[0]==-id[1]+1)||
       (id[1]<-10&&id[1]>-18&&id[0]==-id[1]-1)) allowed=true;
  }
  return allowed;
}

// the decay mode
unsigned int LeptonNeutrinoCurrent::decayMode(vector<int> idout) {
  unsigned int imode=((abs(idout[0])+abs(idout[0])%2)-12)/2;
  return imode;
}

// output the information for the database
void LeptonNeutrinoCurrent::dataBaseOutput(ofstream & output,bool header,
					   bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::LeptonNeutrinoCurrent " << name() 
		    << "  HwWeakCurrents.so\n";
  WeakDecayCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}

}
