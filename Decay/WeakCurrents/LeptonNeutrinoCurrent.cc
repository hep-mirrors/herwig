// -*- C++ -*-
//
// LeptonNeutrinoCurrent.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LeptonNeutrinoCurrent class.
//

#include "LeptonNeutrinoCurrent.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/HelicityFunctions.h"

namespace Herwig {
using namespace ThePEG;
using Helicity::SpinorWaveFunction;
using Helicity::SpinorBarWaveFunction;
using ThePEG::Helicity::LorentzPolarizationVector;
using Helicity::Direction;
using Helicity::incoming;
using Helicity::outgoing;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<LeptonNeutrinoCurrent,WeakCurrent>
describeHerwigLeptonNeutrinoCurrent("Herwig::LeptonNeutrinoCurrent", "HwWeakCurrents.so");

void LeptonNeutrinoCurrent::Init() {

  static ClassDocumentation<LeptonNeutrinoCurrent> documentation
    ("The LeptonNeutrinoCurrent class is designed to handle the "
     "leptonic decay of the weak current.");
}


// complete the construction of the decay mode for integration
bool LeptonNeutrinoCurrent::createMode(int icharge, tcPDPtr ,
				       IsoSpin::IsoSpin Itotal, IsoSpin::I3 i3,
				       unsigned int imode,PhaseSpaceModePtr mode,
				       unsigned int iloc,int ires,
				       PhaseSpaceChannel phase, Energy upp ) {
  // no isospin here
  if(Itotal!=IsoSpin::IUnknown || i3 !=IsoSpin::I3Unknown) return false;
  // make sure the the decays are kinematically allowed
  Energy min =
    getParticleData(11+2*int(imode))->mass()+
    getParticleData(12+2*int(imode))->mass();
  if(min>=upp) return false;
  // set the resonances and check charge
  tPDPtr res;
  if(icharge==3)       res=getParticleData(ParticleID::Wplus);
  else if(icharge==-3) res=getParticleData(ParticleID::Wminus);
  else                 return false;
  // create the channel
  mode->addChannel((PhaseSpaceChannel(phase),ires,res,ires+1,iloc+1,ires+1,iloc+2));
  // return if successful
  return true;
}

// the particles produced by the current
tPDVector LeptonNeutrinoCurrent::particles(int icharge, unsigned int imode_in,
					  int,int) {
  int imode = imode_in;
  tPDVector output(2);
  if(icharge==3) {
    int id = -11-2*imode;
    output[0]=getParticleData(id);
    output[1]=getParticleData(12+2*imode);
  }
  else if(icharge==-3) {
    output[0]=getParticleData(11+2*imode);
    int id = -12-2*imode;
    output[1]=getParticleData(id);
  }
  return output;
}

void LeptonNeutrinoCurrent::constructSpinInfo(ParticleVector decay) const {
  if(decay[0]->id()>0) {
    SpinorWaveFunction   ::constructSpinInfo(wave_   ,decay[1],outgoing,true);
    SpinorBarWaveFunction::constructSpinInfo(wavebar_,decay[0],outgoing,true);
  }
  else {
    SpinorWaveFunction   ::constructSpinInfo(   wave_,decay[0],outgoing,true);
    SpinorBarWaveFunction::constructSpinInfo(wavebar_,decay[1],outgoing,true);
  }
}

// hadronic current   
vector<LorentzPolarizationVectorE> 
LeptonNeutrinoCurrent::current(tcPDPtr ,
			       IsoSpin::IsoSpin Itotal, IsoSpin::I3 i3,
			       const int imode, const int ichan,Energy & scale, 
			       const tPDVector & outgoing,
			       const vector<Lorentz5Momentum> & momenta,
			       DecayIntegrator2::MEOption) const {
  // no isospin here
  if(Itotal!=IsoSpin::IUnknown || i3 !=IsoSpin::I3Unknown) return vector<LorentzPolarizationVectorE>();
  useMe();
  Lorentz5Momentum q = momenta[0]+momenta[1];
  q.rescaleMass();
  scale=q.mass();
  wave_.resize(2);
  wavebar_.resize(2);
  // their wavefunctions
  if(outgoing[0]->id()>0) {
    for(unsigned int ihel=0;ihel<2;++ihel) {
      wavebar_[ihel] = HelicityFunctions::dimensionedSpinorBar(-momenta[0],ihel,Helicity::outgoing);
      wave_   [ihel] = HelicityFunctions::dimensionedSpinor   (-momenta[1],ihel,Helicity::outgoing);
    }
  }
  else {
    for(unsigned int ihel=0;ihel<2;++ihel) {
      wavebar_[ihel] = HelicityFunctions::dimensionedSpinorBar(-momenta[1],ihel,Helicity::outgoing);
      wave_   [ihel] = HelicityFunctions::dimensionedSpinor   (-momenta[0],ihel,Helicity::outgoing);
    }
  }
  // storage for the currents
  vector<LorentzPolarizationVectorE> temp(4);
  // now compute the currents
  int iloc(0);
  unsigned int ix,iy;
  for(ix=0;ix<2;++ix) {
    for(iy=0;iy<2;++iy) {
      iloc = outgoing[0]->id()>0 ? 2*iy+ix : 2*ix+iy;
      temp[iloc]=2.*wave_[ix].leftCurrent(wavebar_[iy]);
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
  WeakCurrent::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}

}
